// AdcSuppressSignalFindingService_service.cc

#include "AdcSuppressSignalFindingService.h"
#include <iostream>
#include <sstream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/DuneInterface/AdcSuppressService.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using art::ServiceHandle;

//**********************************************************************

AdcSuppressSignalFindingService::
AdcSuppressSignalFindingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_psup(&*ServiceHandle<AdcSuppressService>()) { }

//**********************************************************************

int AdcSuppressSignalFindingService::find(AdcChannelData& data) const {
  const string myname = "AdcSuppressSignalFindingService:find: ";
  const AdcCountVector& raw = data.raw;
  if ( raw.size() == 0 && data.samples.size() != 0 ) {
    cout << "ERROR: Input data does not include raw counts." << endl;
    return 1;
  }
  if ( raw.size() != data.samples.size() ) {
    cout << "ERROR: Input prep and raw data have inconsistent sizes." << endl;
    return 2;
  }
  m_psup->filter(raw, data.channel, data.pedestal, data.signal);
  return 0;
}

//**********************************************************************

ostream& AdcSuppressSignalFindingService::
print(ostream& out, string prefix) const {
  out << prefix << "AdcSuppressSignalFindingService:" << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(AdcSuppressSignalFindingService, AdcSignalFindingService)

//**********************************************************************
