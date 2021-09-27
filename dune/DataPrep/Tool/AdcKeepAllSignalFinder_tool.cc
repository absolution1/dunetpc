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
: m_LogLevel(ps.get<int>("LogLevel")),
  m_DataPath(ps.get<Name>("DataPath")) {
  const string myname = "AdcKeepAllSignalFinder::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "      DataPath: " << m_DataPath << endl;
  }
}

//**********************************************************************

DataMap AdcKeepAllSignalFinder::update(AdcChannelData& acdtop) const {
  const string myname = "AdcKeepAllSignalFinder::update: ";
  DataMap ret;
  Index nvie = acdtop.viewSize(m_DataPath);
  Index nroi = 0;
  for ( Index ivie=0; ivie<nvie; ++ivie ) {
    AdcChannelData* pacd = acdtop.mutableViewEntry(m_DataPath, ivie);
    if ( pacd == nullptr ) {
      cout << myname << "WARNING: Skipping null data." << endl;
      continue;
    }
    AdcChannelData& acd = *pacd;
    AdcIndex nsam = acd.samples.size();
    acd.signal.clear();
    acd.rois.clear();
    if ( nsam == 0 ) {
      cout << myname << "ERROR: No samples found in channel " << acd.channel()
           << " view " << ivie << endl;
      return ret.setStatus(1);
    }
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Building ROI for channel " << acd.channel()
           << " view " << ivie << "/" << nvie << "." << endl;
    }
    acd.signal.resize(nsam, true);
    acd.roisFromSignal();
    nroi += acd.rois.size();
  }
  ret.setInt("nroi", nroi);
  return ret;
}

//**********************************************************************

DataMap AdcKeepAllSignalFinder::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcKeepAllSignalFinder)
