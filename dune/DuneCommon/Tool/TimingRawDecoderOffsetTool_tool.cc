// TimingRawDecoderOffsetTool_tool.cc

#include "TimingRawDecoderOffsetTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

using Offset = TimingRawDecoderOffsetTool::Offset;

//**********************************************************************

TimingRawDecoderOffsetTool::TimingRawDecoderOffsetTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_TpcTickPhase(ps.get<Index>("TpcTickPhase")),
  m_Unit(ps.get<Name>("Unit")),
  m_RunDataTool(ps.get<Name>("RunDataTool")) {
  const Name myname = "TimingRawDecoderOffsetTool::ctor: ";
  if ( m_RunDataTool.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_pRunDataTool = ptm->getShared<RunDataTool>(m_RunDataTool);
    if ( m_pRunDataTool == nullptr ) {
      cout << myname << "WARNING: RunDataTool not found: " << m_RunDataTool << endl;
    }
  }
  if ( m_LogLevel ) {
    cout << myname << "Configuration:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "  TpcTickPhase: " << m_TpcTickPhase << endl;
    cout << myname << "          Unit: " << m_Unit << endl;
    cout << myname << "   RunDataTool: " << m_RunDataTool << " @ " << m_pRunDataTool << endl;
  }
}

//**********************************************************************

Offset TimingRawDecoderOffsetTool::offset(const Data& dat) const {
  const Name myname = "TimingRawDecoderOffsetTool::offset: ";
  Name ifname = "artdaqTimestamp-Run" + std::to_string(dat.run) + 
                "-Event" + std::to_string(dat.event) + ".dat";
  Offset res;
  unsigned long daqVal = dat.triggerClock;
  static Index checkCount = 0;
  if ( checkCount ) {
    ifstream fin(ifname.c_str());
    if ( ! fin ) {
      cout << myname << "Unable to find time offset file: " << ifname << endl;
      //return res.setStatus(1);
    }
    unsigned long daqValFile = 0;
    fin >> daqValFile;
    if ( daqVal == daqValFile ) {
      cout << myname << "Input clock matches file clock." << endl;
    } else {
      cout << myname << "Input clock does not match file: " << daqVal << " != " << daqValFile << endl;
    }
    --checkCount;
  }
  if ( m_Unit == "daq" ) {
    res.value = daqVal;
    res.rem = 0.0;
  } else if ( m_Unit == "ns" ) {
    res.value = 20*daqVal;
    res.rem = 0.0;
  } else if ( m_Unit == "tick" ) {
    Index runPhase = 0;
    if ( m_pRunDataTool != nullptr ) {
      Index icha = dat.channel;
      RunData rdat = m_pRunDataTool->runData(dat.run, dat.subrun);
      if ( rdat.havePhaseGroup() ) {
        bool haveGroup = true;
        Index igrp = 0;
        art::ServiceHandle<dune::PdspChannelMapService> pchanMap;
        bool useWib = false;
        bool useFemb = false;
        if ( rdat.phaseGroup() == "channel" ) igrp = icha;
        else if ( rdat.phaseGroup() == "all" ) igrp = 0;
        else if ( rdat.phaseGroup() == "wib" ) useWib = true;
        else if ( rdat.phaseGroup() == "femb" ) useFemb = true;
        else {
          haveGroup = false;
          cout << myname << "WARNING: Invalid phase group: " << rdat.phaseGroup() << endl;
        }
        if ( haveGroup ) {
          art::ServiceHandle<dune::PdspChannelMapService> pchanMap;
          if ( pchanMap.get() != nullptr ) {
            Index kapa = pchanMap->APAFromOfflineChannel(icha);
            Index jwib = pchanMap->WIBFromOfflineChannel(icha);
            Index jcon = pchanMap->FEMBFromOfflineChannel(icha);
            Index jfmb = 5*jcon + jwib;
            Index kfmb = 20*kapa + jfmb;
            Index kwib = 5*kapa + jfmb%5;
            if ( useWib  ) igrp = kwib;
            if ( useFemb ) igrp = kfmb;
          } else {
            haveGroup = false;
            cout << myname << "WARNING: Channel map service not found." << endl;
          }
        }
        if ( haveGroup ) {
          const RunData::IndexVector& phases = rdat.phases();
          if ( igrp + 1 > phases.size() ) {
            cout << myname << "WARNING: Phases is too short: " << igrp << "/" << phases.size() << endl;
          } else {
            runPhase = phases[igrp];
          }
        }
      }
      if ( m_LogLevel >= 3 ) cout << myname << "Run phase is " << runPhase << endl;
    } else {
      if ( m_LogLevel >= 3 ) cout << myname << "Run data tool not found." << endl;
    }
    long daqoff = daqVal + m_TpcTickPhase + runPhase;
    res.value = daqoff/25;
    res.rem = (daqoff % 25)/25.0;
  } else {
    cout << myname << "Invalid unit: " << m_Unit << ifname << endl;
    return res.setStatus(2);
  }
  res.unit = m_Unit;
  if ( m_LogLevel >= 2 ) cout << myname << "Offset is " << res.value
                              << " " << res.unit << endl;
  return res;
}

//**********************************************************************
