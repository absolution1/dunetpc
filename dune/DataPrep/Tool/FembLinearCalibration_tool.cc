// FembLinearCalibration_tool.cc

#include "FembLinearCalibration.h"
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

FembLinearCalibration::FembLinearCalibration(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Units(ps.get<string>("Units")),
  m_Gains(ps.get<AdcSignalVector>("Gains")),
  m_AdcMin(ps.get<AdcCount>("AdcMin")),
  m_AdcMins(ps.get<AdcCountVector>("AdcMins")),
  m_AdcMax(ps.get<AdcCount>("AdcMax")),
  m_AdcMaxs(ps.get<AdcCountVector>("AdcMaxs")) {
  const string myname = "FembLinearCalibration::ctor: ";
  int w = 28;
  if ( m_LogLevel >= 1 ) {
    cout << myname << setw(w) << "Log level: " << m_LogLevel << endl;
    ostringstream sslab;
    sslab << "Gains[" << m_Gains.size() << "]" << (m_Gains.size() ? ": " : "  ");;
    cout << myname << setw(w) << sslab.str();
    for ( Index icha=0; icha<m_Gains.size(); ++icha ) {
      if ( icha > 0 ) cout << ", ";
      if ( icha >= 10 ) { cout << "..."; break; }
      cout << m_Gains[icha];
    }
    cout << endl;
    cout << myname << setw(w) << "Default ADC min: " << m_AdcMin << endl;
    sslab.str("");
    sslab << "Channel ADC min[" << m_AdcMins.size() << "]" << (m_AdcMins.size() ? ": " : "  ");;
    cout << myname << setw(w) << sslab.str();
    for ( Index icha=0; icha<m_AdcMins.size(); ++icha ) {
      if ( icha > 0 ) cout << ", ";
      if ( icha >= 10 ) { cout << "..."; break; }
      cout << m_AdcMins[icha];
    }
    cout << endl;
    cout << myname << setw(w) << "Default ADC max: " << m_AdcMax << endl;
    sslab.str("");
    sslab << "Channel ADC max[" << m_AdcMaxs.size() << "]" << (m_AdcMaxs.size() ? ": " : "  ");
    cout << myname << setw(w) << sslab.str();
    for ( Index icha=0; icha<m_AdcMaxs.size(); ++icha ) {
      if ( icha > 0 ) cout << ", ";
      if ( icha >= 10 ) { cout << "..."; break; }
      cout << m_AdcMaxs[icha];
    }
    cout << endl;
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
  if ( m_AdcMaxs.size() && icha >= m_AdcMaxs.size() ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "ADC max not found for channel " << icha << endl;
    }
    return res.setStatus(6);
  }
  AdcSignal gain = m_Gains[icha];
  if ( gain > 0.0 ) {
    acd.samples.resize(acd.raw.size(), 0.0);
    acd.flags.resize(acd.raw.size(), AdcGood);
    AdcCount adcmin = m_AdcMins.size() ? m_AdcMins[icha] : m_AdcMin;
    AdcCount adcmax = m_AdcMaxs.size() ? m_AdcMaxs[icha] : m_AdcMax;
    Index nunder = 0;
    Index nover = 0;
    Index nsam = acd.raw.size();
    for ( Index isam=0; isam<nsam; ++isam ) {
      acd.samples[isam] = gain*(acd.raw[isam] - acd.pedestal);
      if ( acd.raw[isam] <= adcmin ) {
        acd.flags[isam] = AdcUnderflow;
        ++nunder;
      } else if ( acd.raw[isam] >= adcmax ) {
        acd.flags[isam] = AdcOverflow;
        ++nover;
      }
    }
    acd.sampleUnit = m_Units;
    res.setInt("calibSampleCount", nsam);
    res.setInt("calibUnderflowCount", nunder);
    res.setInt("calibOverflowCount", nover);
    res.setInt("calibAdcMin", adcmin);
    res.setInt("calibAdcMax", adcmax);
  } else {
    acd.samples.resize(0);
  }
  res.setFloat("calibGain", gain);
  return res;
}

//**********************************************************************
