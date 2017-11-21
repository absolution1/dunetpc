// FembLinearCalibration_tool.cc

#include "FembLinearCalibration.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

using Index = unsigned int;

//**********************************************************************

FembLinearCalibration::FembLinearCalibration(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Gains(ps.get<AdcSignalVector>("Gains")),
  m_AdcMins(ps.get<AdcCountVector>("AdcMins")) {
  const string myname = "FembLinearCalibration::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "      Log level: " << m_LogLevel << endl;
    cout << myname << "     Gain count: " << m_Gains.size() << endl;
    cout << myname << "  ADC min count: " << m_AdcMins.size() << endl;
  }
}

//**********************************************************************

DataMap FembLinearCalibration::view(const AdcChannelData& acd) const {
  DataMap result;
  AdcChannelData acdtmp(acd);
  return update(acdtmp);
}

//**********************************************************************

DataMap FembLinearCalibration::update(AdcChannelData& acd) const {
  const string myname = "FembLinearCalibration::update: ";
  DataMap res;
  AdcChannel icha = acd.channel;
  if ( icha == AdcChannelData::badChannel ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Data does not have a channel ID." << endl;
    }
    return res.setStatus(2);
  }
  if ( icha >= m_Gains.size() ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Gain not found for channel " << icha << endl;
    }
    return res.setStatus(4);
  }
  if ( m_AdcMins.size() && icha >= m_AdcMins.size() ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "ADC min not found for channel " << icha << endl;
    }
    return res.setStatus(5);
  }
  AdcSignal gain = m_Gains[icha];
  AdcCount adcmin = m_AdcMins[icha];
  acd.samples.resize(acd.raw.size(), 0.0);
  acd.flags.resize(acd.raw.size(), AdcGood);
  Index nunder = 0;
  Index nsam = acd.raw.size();
  for ( Index isam=0; isam<nsam; ++isam ) {
    acd.samples[isam] = gain*(acd.raw[isam] - acd.pedestal);
    if ( acd.raw[isam] <= adcmin ) {
      acd.flags[isam] = AdcUnderflow;
      ++nunder;
    }
  }
  res.setInt("nSample", nsam);
  res.setInt("nUnderflow", nunder);
  return res;
}

//**********************************************************************
