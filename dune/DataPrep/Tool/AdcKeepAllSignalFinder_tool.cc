// AdcKeepAllSignalFinder_tool.cc

#include "AdcKeepAllSignalFinder.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcKeepAllSignalFinder::AdcKeepAllSignalFinder(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "AdcKeepAllSignalFinder::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
  }
}

//**********************************************************************

DataMap AdcKeepAllSignalFinder::update(AdcChannelData& acd) const {
  const string myname = "AdcKeepAllSignalFinder::update: ";
  DataMap ret;
  AdcIndex nsam = acd.samples.size();
  acd.signal.clear();
  acd.rois.clear();
  if ( nsam == 0 ) {
    cout << myname << "ERROR: No samples found in channel " << acd.channel << endl;
    return ret.setStatus(1);
  }
  if ( m_LogLevel >= 2 ) cout << myname << "Building ROI for channel " << acd.channel << "." << endl;
  acd.signal.resize(nsam, true);
  acd.roisFromSignal();
  ret.setInt("nroi", acd.rois.size());
  return ret;
}

//**********************************************************************

DataMap AdcKeepAllSignalFinder::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************
