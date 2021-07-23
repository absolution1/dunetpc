// FloatArrayGainCalibration_tool.cc

#include "FloatArrayGainCalibration.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Data/RunData.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune/DuneInterface/Tool/FloatArrayTool.h"
#include "dune/DuneCommon/Utility/RootParFormula.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;

using Index = unsigned int;

//**********************************************************************

FloatArrayGainCalibration::FloatArrayGainCalibration(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Unit(ps.get<string>("Unit")),
  m_GainDefault(ps.get<float>("GainDefault")),
  m_AdcUnderflowDefault(ps.get<unsigned int>("AdcUnderflowDefault")),
  m_AdcOverflowDefault(ps.get<unsigned int>("AdcOverflowDefault")),
  m_GainTool(ps.get<string>("GainTool")),
  m_ScaleFactor(new RootParFormula("ScaleFactor", ps.get<Name>("ScaleFactor"))),
  m_prdtool(nullptr) {
  const string myname = "FloatArrayGainCalibration::ctor: ";
  DuneToolManager* pdtm = DuneToolManager::instance();
  if ( pdtm == nullptr ) {
    cout << myname << "ERROR: Unable to retrieve tool manager." << endl;
  } else {
    m_pgains = pdtm->getShared<FloatArrayTool>(m_GainTool);
    if ( ! m_pgains ) {
      cout << myname << "ERROR: Unable to retrieve gains tool " << m_GainTool << endl;
    }
    string stnam = "runDataTool";
    if ( m_ScaleFactor->npar() ) {
      m_prdtool = pdtm->getShared<RunDataTool>(stnam);
      if ( m_prdtool == nullptr ) {
        cout << myname << "ERROR: RunDataTool " << stnam
             << " not found. Scale factor formula will not be evaluated." << endl;
      } else {
        cout << myname << "RunDataTool retrieved." << endl;
      }
    }
  }
  m_ScaleFactor->setDefaultEval(1.0);
  if ( m_LogLevel >= 1 ) {
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "          Unit: " << m_Unit << endl;
    cout << myname << "   GainDefault: " << m_GainDefault << endl;
    cout << myname << "      GainTool: " << m_GainTool  << " (@" << m_pgains << ")" << endl;
    cout << myname << "   ScaleFactor: " << m_ScaleFactor->formulaString() << endl;
  }
}

//**********************************************************************

DataMap FloatArrayGainCalibration::view(const AdcChannelData& acd) const {
  DataMap result;
  AdcChannelData acdtmp(acd);
  acdtmp.raw = acd.raw;
  return update(acdtmp);
}

//**********************************************************************

DataMap FloatArrayGainCalibration::update(AdcChannelData& acd) const {
  const string myname = "FloatArrayGainCalibration::update: ";
  DataMap res;
  if ( ! m_pgains ) return res;
  AdcChannel icha = acd.channel();
  if ( icha == AdcChannelData::badChannel() ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Data does not have a channel ID." << endl;
    }
    return res.setStatus(2);
  }
  const FloatArrayTool& gains = *m_pgains;
  if ( m_LogLevel >= 2 ) {
    if ( ! gains.inRange(icha) ) {
      cout << myname << "Gain not found for channel " << icha << endl;
    }
  }
  if ( m_prdtool != nullptr ) {
    RunData rdat = m_prdtool->runData(acd.run());
    if ( ! rdat.isValid() ) cout << myname << "WARNING: RunData not found." << endl;
    else rdat.setFormulaPars(*m_ScaleFactor);
  }
  if ( ! m_ScaleFactor->ready() ) {
    cout << myname << "WARNING: Using default scale factor " << m_ScaleFactor->defaultEval() << endl;
  } 
  float gain = m_GainDefault >= 0 ? gains.value(icha, m_GainDefault) : gains.value(icha);
  gain *= m_ScaleFactor->eval();
  AdcCount adcudr = m_AdcUnderflowDefault;
  AdcCount adcovr = m_AdcOverflowDefault;
  acd.samples.resize(acd.raw.size(), 0.0);
  acd.flags.resize(acd.raw.size(), AdcGood);
  AdcCount adcmin = 0;
  AdcCount adcmax = 0;
  Index nunder = 0;
  Index nover = 0;
  Index nsam = acd.raw.size();
  Name unit = m_Unit;
  if ( m_LogLevel >= 3 ) cout << myname << "Processing " << nsam << " samples for channel " << icha << endl;
  for ( Index isam=0; isam<nsam; ++isam ) {
    AdcCount adcin = acd.raw[isam];
    if ( adcin <= adcudr ) {
      acd.flags[isam] = AdcUnderflow;
      ++nunder;
    } else if ( adcin >= adcovr ) {
      acd.flags[isam] = AdcOverflow;
      ++nover;
    }
    if ( isam ==0 || adcin < adcmin ) adcmin = adcin;
    if ( isam ==0 || adcin > adcmax ) adcmax = adcin;
    float sigout = gain*(adcin - acd.pedestal);
    acd.samples[isam] = sigout;
  }
  acd.sampleUnit = unit;
  res.setFloat("calibGain", gain);
  res.setInt("calibSampleCount", nsam);
  res.setInt("calibUnderflowCount", nunder);
  res.setInt("calibOverflowCount", nover);
  res.setInt("calibAdcMin", adcmin);
  res.setInt("calibAdcMax", adcmax);
  return res;
}

//**********************************************************************
