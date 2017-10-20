// AdcThresholdSignalFinder_tool.cc

#include "AdcThresholdSignalFinder.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcThresholdSignalFinder::AdcThresholdSignalFinder(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Threshold(ps.get<float>("Threshold")),
  m_BinsBefore(ps.get<unsigned int>("BinsBefore")),
  m_BinsAfter(ps.get<unsigned int>("BinsAfter")),
  m_FlagPositive(ps.get<bool>("FlagPositive")),
  m_FlagNegative(ps.get<bool>("FlagNegative")) {
  const string myname = "AdcThresholdSignalFinder::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "     Threshold: " << m_Threshold<< endl;
    cout << myname << "    BinsBefore: " << m_BinsBefore << endl;
    cout << myname << "     BinsAfter: " << m_BinsAfter << endl;
    cout << myname << "  FlagPositive: " << m_FlagPositive << endl;
    cout << myname << "  FlagNegative: " << m_FlagNegative << endl;
  }
}

//**********************************************************************

DataMap AdcThresholdSignalFinder::update(AdcChannelData& acd) const {
  AdcIndex nsam = acd.samples.size();
  AdcIndex nsamlo = m_BinsBefore;
  AdcIndex nsamhi = m_BinsAfter;
  acd.signal.resize(nsam, false);
  AdcIndex nbinAbove = 0;
  for ( AdcIndex isam=0; isam<nsam; ++isam ) {
    bool keep = ( m_FlagPositive && acd.samples[isam] >  m_Threshold ) ||
                ( m_FlagNegative && acd.samples[isam] < -m_Threshold );
    if ( keep ) {
      ++nbinAbove;
      AdcIndex jsam1 = isam > nsamlo ? isam - nsamlo : 0;
      AdcIndex jsam2 = isam + nsamhi + 1;
      if ( jsam2 > nsam ) jsam2 = nsam;
      for ( AdcIndex jsam=jsam1; jsam<jsam2; ++jsam ) acd.signal[jsam] = true;
    }
  }
  acd.roisFromSignal();
  DataMap res(0);
  res.setInt("nThresholdBins", nbinAbove);
  return res;
}

//**********************************************************************

DataMap AdcThresholdSignalFinder::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************
