// AdcSampleScaler_tool.cc

#include "AdcSampleScaler.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcSampleScaler::AdcSampleScaler(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_ScaleFactor(ps.get<float>("ScaleFactor")),
  m_OutputUnit(ps.get<string>("OutputUnit")),
  m_InputUnit(ps.get<string>("InputUnit")) {
  const string myname = "AdcSampleScaler::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "      ScaleFactor: " << m_ScaleFactor << endl;
    cout << myname << "       OutputUnit: " << m_OutputUnit << endl;
    cout << myname << "        InputUnit: " << m_InputUnit << endl;
  }
}

//**********************************************************************

DataMap AdcSampleScaler::update(AdcChannelData& acd) const {
  const string myname = "AdcSampleScaler::update: ";
  if ( m_LogLevel >= 2 ) cout << "Processing run " << acd.run() << " event " << acd.event()
                              << " channel " << acd.channel() << endl;
  DataMap ret;

  // Check input data unit.
  int unitCheck = 0;
  if ( m_InputUnit.size() ) {
    if ( acd.sampleUnit.size() == 0 ) {
      cout << myname << "WARNING: Input data does not have a sample unit." << endl;
      unitCheck = 1;
    } else if ( acd.sampleUnit != m_InputUnit ) {
      cout << myname << "WARNING: Unexpected input data unit: " << acd.sampleUnit
           << " != " << m_InputUnit << endl;
      unitCheck = 2;
    }
  }

  // Scale samples.
  for ( float& sam : acd.samples ) sam *= m_ScaleFactor;
  if ( m_OutputUnit.size() ) acd.sampleUnit = m_OutputUnit;

  ret.setInt("acsUnitCheck", unitCheck);
  ret.setInt("acsSampleCount", acd.samples.size());

  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcSampleScaler)
