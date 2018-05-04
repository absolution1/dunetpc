// AdcDetectorPlotter_tool.cc

#include "AdcDetectorPlotter.h"
#include <iostream>
#include <sstream>
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneCommon/RootPalette.h"
#include "dune/DuneCommon/LineColors.h"
#include "larcore/Geometry/Geometry.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::ostringstream;

using Tick = AdcSignalVector::size_type;
using geo::GeometryCore;
using std::vector;

//**********************************************************************
// Local methods.
//**********************************************************************

namespace {

void initializeState(AdcDetectorPlotter::State& state, const AdcChannelData& acd) {
  state.reportCount = 0;
  state.channelCount = 0;
  state.run    = acd.run;
  state.subrun = acd.subRun;
  state.event  = acd.event;
  state.ppad.reset(nullptr);
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDetectorPlotter::AdcDetectorPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_WireAngle(ps.get<float>("WireAngle")),
  m_DataType(ps.get<int>("DataType")),
  m_Tick0(ps.get<float>("Tick0")),
  m_DriftSpeed(ps.get<float>("DriftSpeed")),
  m_XMin(ps.get<float>("XMin")),
  m_XMax(ps.get<float>("XMax")),
  m_ZMin(ps.get<float>("ZMin")),
  m_ZMax(ps.get<float>("ZMax")),
  m_SignalThreshold(ps.get<float>("SignalThreshold")),
  m_ShowWires(ps.get<bool>("ShowWires")),
  m_ShowCathode(ps.get<bool>("ShowCathode")),
  m_ShowGrid(ps.get<bool>("ShowGrid")),
  m_Title(ps.get<string>("Title")),
  m_FileName(ps.get<string>("FileName")),
  m_state(new State) {
  const string myname = "AdcDetectorPlotter::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "        WireAngle: " << m_WireAngle << endl;
    cout << myname << "         DataType: " << m_DataType << endl;
    cout << myname << "            Tick0: " << m_Tick0 << endl;
    cout << myname << "       DriftSpeed: " << m_DriftSpeed << " cm/tick" << endl;
    cout << myname << "             XMin: " << m_XMin << " cm" << endl;
    cout << myname << "             XMax: " << m_XMax << " cm" << endl;
    cout << myname << "             ZMin: " << m_ZMin << " cm" << endl;
    cout << myname << "             ZMax: " << m_ZMax << " cm" << endl;
    cout << myname << "  SignalThreshold: " << m_SignalThreshold << endl;
    cout << myname << "        ShowWires: " << m_ShowWires << endl;
    cout << myname << "      ShowCathode: " << m_ShowCathode << endl;
    cout << myname << "         ShowGrid: " << m_ShowGrid << endl;
    cout << myname << "            Title: " << m_Title << endl;
    cout << myname << "         FileName: " << m_FileName << endl;
  }
  WireSelector& sel = getState()->sel;
  sel.selectWireAngle(m_WireAngle);
  const WireSelector::WireInfoVector& wdat = sel.fillData();
  const WireSelector::WireInfoMap& wmap = sel.fillDataMap();
  const WireSelector::WireSummary& wsum = sel.fillWireSummary();
  const GeometryCore* pgeo = sel.geometry();
  string gname = pgeo == nullptr ? "NONE" : pgeo->DetectorName();
  if ( m_LogLevel ) {
    cout << myname << "Derived: " << endl;
    cout << myname << "           Geometry name: " << gname << endl;
    cout << myname << "       # selected planes: " << sel.planeIDs().size() << endl;
    cout << myname << "        # selected wires: " << wdat.size() << endl;
    cout << myname << "  # selected wire points: " << wsum.size() << endl;
    cout << myname << "     # selected channels: " << wmap.size() << endl;
  }
}

//**********************************************************************

void AdcDetectorPlotter::initialize() {
}

//**********************************************************************

DataMap AdcDetectorPlotter::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcDetectorPlotter::view: ";
  DataMap ret;
  State& state = *getState();
  if ( m_LogLevel >= 2 ) cout << myname << "Begin call " << state.reportCount
                              << "/" << state.jobCount << "." << endl;
  if ( acds.size() == 0 ) {
    cout << myname << "WARNING: Channel map is empty. No data extracted." << endl;
    return ret.setStatus(1);
  }
  const AdcChannelData& acdFirst = acds.begin()->second;
  const AdcChannelData& acdLast = acds.rbegin()->second;
  string hname = "hdet";
  int npadx = 1200;
  int npady = 1000;
  string sttlx = "Drift coordinate [cm]";
  string sttly = "Wire coordinate [cm]";
  if ( state.jobCount == 0 ||
       acdFirst.run != state.run || acdFirst.subRun != state.subrun || acdFirst.event != state.event ) {
    if ( m_LogLevel >= 2 ) cout << myname << "  Starting new event." << endl;
    initializeState(state, acdFirst);
    string sttl = m_Title;
    state.ofname = m_FileName;
    vector<string*> strs = {&sttl, &state.ofname};
    for ( string* pstr : strs ) {
      string& str = *pstr;
      StringManipulator sman(str);
      //sman.replace("%PAT%", fpat);
      if ( acdFirst.run != AdcChannelData::badIndex ) sman.replace("%RUN%", acdFirst.run);
      else sman.replace("%RUN%", "RunNotFound");
      if ( acdFirst.subRun != AdcChannelData::badIndex ) sman.replace("%SUBRUN%", acdFirst.subRun);
      else sman.replace("%SUBRUN%", "SubRunNotFound");
      if ( acdFirst.event != AdcChannelData::badIndex ) sman.replace("%EVENT%", acdFirst.event);
      else sman.replace("%EVENT%", "EventNotFound");
    }
    // Create graph.
    state.ppad.reset(new TPadManipulator(npadx, npady));
    TGraph* pg = new TGraph;
    pg->SetTitle(sttl.c_str());
    pg->GetXaxis()->SetTitle(sttlx.c_str());
    pg->GetYaxis()->SetTitle(sttly.c_str());
    state.ppad->add(pg, "P");
    state.ppad->graph()->Expand(512);  // Allocate 1024 points.
    state.ppad->setRangeX(m_XMin, m_XMax);
    state.ppad->setRangeY(m_ZMin, m_ZMax);
    if ( m_ShowGrid ) state.ppad->setGrid();
    LineColors cols;
    const WireSelector::WireInfoVector& wins = state.sel.fillData();
    if ( m_ShowWires ) {
      TGraph* pgw = new TGraph;
      pgw->SetMarkerColor(cols.red());
      pgw->Expand(wins.size());
      for ( Index iwin=0; iwin<wins.size(); ++iwin ) {
        const WireSelector::WireInfo& win = wins[iwin];
        pgw->SetPoint(iwin, win.x, win.z);
      }
      state.ppad->add(pgw, "P");
    }
    if ( m_ShowCathode ) {
      TGraph* pgc = new TGraph;
      pgc->SetMarkerColor(cols.green());
      pgc->Expand(wins.size());
      for ( Index iwin=0; iwin<wins.size(); ++iwin ) {
        const WireSelector::WireInfo& win = wins[iwin];
        pgc->SetPoint(iwin, win.x + win.driftMax, win.z);
      }
      state.ppad->add(pgc, "P");
    }
  }
  ++state.jobCount;
  ++state.reportCount;
  Tick maxtick = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( iacd.first == AdcChannelData::badChannel ) {
      cout << myname << "WARNING: Channel map has invalid channels. No plot is created." << endl;
    }
    Tick ntick = iacd.second.samples.size();
    if ( ntick > maxtick ) maxtick = ntick;
  }
  AdcIndex chanFirst = acdFirst.channel;
  AdcIndex chanLast = acdLast.channel;
  AdcIndex nchan = chanLast + 1 - chanFirst;
  if ( m_LogLevel >= 2 ) cout << myname << "  Input channel count is " << nchan << endl;
  // Fill graph.
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( m_LogLevel >= 3 ) cout << myname << "    Filling with channel " << iacd.first << endl;
    addChannel(iacd.second);
  }
  if ( m_LogLevel >= 2 ) cout << myname << "  Printing plot " << state.ofname << endl;
  state.ppad->print(state.ofname);
  return ret;
}

//**********************************************************************

int AdcDetectorPlotter::addChannel(const AdcChannelData& acd) const {
  const string myname = "AdcDetectorPlotter::addChannel: ";
  bool isRaw = m_DataType == 1;
  bool isPrep = m_DataType == 0;
  if ( !isRaw && !isPrep ) {
    cout << myname << "ERROR: Invalid data type: " << m_DataType << endl;
    return 2;
  }
  State& state = *getState();
  TGraph* pg = state.ppad->graph();
  if ( pg == nullptr ) return 1;
  AdcChannel icha = acd.channel;
  auto rng = state.sel.dataMap().equal_range(icha);
  for ( auto ient=rng.first; ient!=rng.second; ++ient) {
    const WireSelector::WireInfo& win = *(ient->second);
    float z = win.z;
    float driftVelocity = win.driftSign()*m_DriftSpeed;
    if ( isRaw ) {
    } else {
      for ( Index isam=0; isam<acd.samples.size(); ++isam ) {
        float sig = acd.samples[isam];
        if ( sig > m_SignalThreshold ) {
          float x = win.x + driftVelocity*(isam - m_Tick0);
          Index ipt = pg->GetN();
          pg->SetPoint(ipt, x, z);
          if ( m_LogLevel >= 4 ) {
            ostringstream sout;
            sout.precision(2);
            sout << myname << "Added point " << ipt << " (" << x << ", " << z << ")";
            cout << sout.str() << endl;
          }
        }
      }
    }
  }
/*
  AdcChannel chan = acd.channel;
  const AdcSignalVector& sams = acd.samples;
  const AdcCountVector& raw = acd.raw;
  AdcSignal ped = 0.0;
  bool isRawPed = false;
  if ( isRaw ) {
    ped = acd.pedestal;
    isRawPed = ped != AdcChannelData::badSignal;
  }
    unsigned int ibin = ph->GetBin(1, chan-chanFirst+1);
    for ( Tick isam=0; isam<sams.size(); ++isam, ++ibin ) {
      unsigned int isig = isam + m_FirstTick;
      if ( isPrep ) {
        if ( isig < sams.size() ) ph->SetBinContent(ibin, sams[isig]);
      } else if ( isRawPed ) {
        if ( isig < raw.size() ) ph->SetBinContent(ibin, raw[isig] - ped);
      } else {
        cout << myname << "Fill failed for bin " << ibin << endl;
      }
  }
*/
  return 0;
}

//**********************************************************************
