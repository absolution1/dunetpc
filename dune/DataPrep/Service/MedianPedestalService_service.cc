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
: m_LogLevel(1) {
  const string myname = "MedianPedestalService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  print(cout, myname);
}

//**********************************************************************

int MedianPedestalService::
evaluate(const AdcSignalVector& sigsin, AdcSignal* pped, AdcSignal* prms,
         AdcSignal* ppederr, AdcSignal* prmserr) const {
  AdcSignal ped = 0.0;
  AdcSignal rms = 0.0;
  AdcSignal pederr = 0.0;
  AdcSignal rmserr = 0.0;
  if ( sigsin.size() == 0 ) return 0;
  AdcSignalVector sigs = sigsin;
  sort(sigs.begin(), sigs.end());
  unsigned int isig = sigs.size()/2;
  bool isodd = sigs.size()%2;
  ped = isodd ? sigs[isig-1] : 0.5*(sigs[isig-1] + sigs[isig]);
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
  out << prefix << "               LogLevel: " << m_LogLevel              << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(MedianPedestalService, PedestalEvaluationService)

//**********************************************************************
