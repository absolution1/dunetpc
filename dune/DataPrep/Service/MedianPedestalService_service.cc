// MedianPedestalService_service.cc

#include "MedianPedestalService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;
using art::ServiceHandle;

//**********************************************************************

MedianPedestalService::
MedianPedestalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1),
  m_SkipFlaggedSamples(false),
  m_SkipSignals(false) {
  const string myname = "MedianPedestalService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_UseMean = pset.get<bool>("UseMean");
  m_SkipFlaggedSamples = pset.get<bool>("SkipFlaggedSamples");
  m_SkipSignals = pset.get<bool>("SkipSignals");
  print(cout, myname);
}

//**********************************************************************

int MedianPedestalService::
evaluate(const AdcChannelData& data, AdcSignal* pped, AdcSignal* prms,
         AdcSignal* ppederr, AdcSignal* prmserr) const {
  const string myname = "MedianPedestalService::evaluate: ";
  AdcSignal ped = 0.0;
  AdcSignal rms = 0.0;
  AdcSignal pederr = 0.0;
  AdcSignal rmserr = 0.0;
  const AdcSignalVector sigsin = data.samples;
  const AdcFlagVector flags = data.flags;
  const AdcFilterVector signal = data.signal;
  if ( sigsin.size() == 0 ) return 0;
  AdcSignalVector sigs;
  unsigned int nsig = sigsin.size();
  double sigsum = 0.0;
  unsigned int nsigout = 0;
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    AdcSignal sig = sigsin[isig];
    if ( m_SkipFlaggedSamples && flags.size() ) {
      AdcFlag flag = flags[isig];
      bool skip = flag != AdcGood;
      if ( skip ) continue;
    }
    if ( m_SkipSignals && signal.size() ) {
      bool isSignal = signal[isig];
      if ( isSignal ) continue;
    }
    if ( m_UseMean ) sigsum += sig;
    else sigs.push_back(sig);
    ++nsigout;
  }
  if ( nsigout ) {
    if ( m_UseMean ) {
      ped = sigsum/nsigout;
    } else {
      sort(sigs.begin(), sigs.end());
      unsigned int isig = sigs.size()/2;
      bool isodd = sigs.size()%2;
      ped = isodd ? sigs[isig-1] : 0.5*(sigs[isig-1] + sigs[isig]);
    }
  }
  if ( pped != nullptr ) *pped = ped;
  if ( prms != nullptr ) *prms = rms;
  if ( ppederr != nullptr ) *ppederr = pederr;
  if ( prmserr != nullptr ) *prmserr = rmserr;
  return 0;
}

//**********************************************************************

ostream& MedianPedestalService::
print(ostream& out, string prefix) const {
  out << prefix << "MedianPedestalService:"                          << endl;
  out << prefix << "             LogLevel: " << m_LogLevel           << endl;
  out << prefix << "   SkipFlaggedSamples: " << m_SkipFlaggedSamples << endl;
  out << prefix << "          SkipSignals: " << m_SkipSignals        << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(MedianPedestalService, PedestalEvaluationService)

//**********************************************************************
