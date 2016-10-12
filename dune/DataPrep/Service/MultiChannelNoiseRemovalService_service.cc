// MultiChannelNoiseRemovalService_service.cc

#include "MultiChannelNoiseRemovalService.h"
#include "dune/DuneInterface/AdcChannelNoiseRemovalService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"

using std::cout;
using std::endl;
using std::ostream;
using std::string;
using art::ServiceHandle;

//**********************************************************************

MultiChannelNoiseRemovalService::
MultiChannelNoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "MultiChannelNoiseRemovalService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  // Get service.
  m_pAdcChannelNoiseRemovalService = &*ServiceHandle<AdcChannelNoiseRemovalService>();
  print(cout, myname);
}

//**********************************************************************

int MultiChannelNoiseRemovalService::update(AdcChannelDataMap& datamap) const {
  const string myname = "MultiChannelNoiseRemovalService:update: ";
  if ( m_pAdcChannelNoiseRemovalService == nullptr ) {
    cout << myname << "ERROR: AdcChannelNoiseRemovalService not found." << endl;
    return 1;
  }
  if ( m_LogLevel >= 2 ) cout << myname << "Updating " << datamap.size() << " channel"
                              << (datamap.size()!=1 ? "" : "s") << "." << endl;
  int rstat = 0;
  for ( AdcChannelDataMap::value_type& chdata : datamap ) {
    AdcChannelData& data = chdata.second;
    if ( m_pAdcChannelNoiseRemovalService->update(data) ) ++rstat;
  }
  if ( rstat != 0 ) cout << myname << "WARNING: Update failed for " << rstat << " of "
                         << datamap.size() << " channel" << (rstat!=1 ? "" : "s") << "." << endl;
  return 0;
}

//**********************************************************************

ostream& MultiChannelNoiseRemovalService::
print(ostream& out, string prefix) const {
  out << prefix << "MultiChannelNoiseRemovalService:"  << endl;
  out << prefix << "     LogLevel: " << m_LogLevel     << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(MultiChannelNoiseRemovalService, AdcNoiseRemovalService)

//**********************************************************************
