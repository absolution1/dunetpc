// ThresholdNoiseRemovalService_service.cc

#include "ThresholdNoiseRemovalService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;
using art::ServiceHandle;

//**********************************************************************

ThresholdNoiseRemovalService::
ThresholdNoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "ThresholdNoiseRemovalService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_Threshold         = pset.get<AdcSignal>("Threshold");
  print(cout, myname);
}

//**********************************************************************

int ThresholdNoiseRemovalService::update(AdcChannelData& data) const {
  const string myname = "ThresholdNoiseRemovalService:update: ";
  AdcSignalVector& sigs = data.samples;
  AdcFlagVector& flags = data.flags;
  if ( sigs.size() == 0 ) {
    cout << myname << "WARNING: No data found." << endl;
    return 1;
  }
  if ( flags.size() == 0 ) flags.resize(sigs.size(), AdcGood );
  if ( flags.size() != sigs.size() ) {
    cout << myname << "ERROR: samples and flags vectors have different sizes." << endl;
    return 2;
  }
  AdcIndex count = 0;
  for ( unsigned int isig=0; isig<sigs.size(); ++isig ) {
    AdcSignal& sig = sigs[isig];
    if ( fabs(sig) < m_Threshold ) {
      sig = 0.0;
      flags[isig] = AdcSetFixed;
      ++count;
    }
  }
  if ( m_LogLevel >= 2 ) cout << myname << "Suppressed " << count << " of "
                              << sigs.size() << " ticks." << endl;
  return 0;
}

//**********************************************************************

ostream& ThresholdNoiseRemovalService::
print(ostream& out, string prefix) const {
  out << prefix << "ThresholdNoiseRemovalService:" << endl;
  out << prefix << "   LogLevel: " << m_LogLevel   << endl;
  out << prefix << "  Threshold: " << m_Threshold  << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(ThresholdNoiseRemovalService, AdcChannelNoiseRemovalService)

//**********************************************************************
