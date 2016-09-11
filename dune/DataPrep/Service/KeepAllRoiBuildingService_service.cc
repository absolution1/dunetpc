// KeepAllRoiBuildingService_service.cc

#include "KeepAllRoiBuildingService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::setw;
using art::ServiceHandle;

//**********************************************************************

KeepAllRoiBuildingService::
KeepAllRoiBuildingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "KeepAllRoiBuildingService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  if ( m_LogLevel > 0 ) print(cout, myname);
}


//**********************************************************************

int KeepAllRoiBuildingService::build(AdcChannelData& data) const {
  const string myname = "KeepAllRoiBuildingService:build: ";
  if ( m_LogLevel >= 2 ) cout << myname << "Building ROIs for channel "
                              << data.channel << "." << endl;
  data.signal.clear();
  data.signal.resize(data.samples.size(), true);
  data.roisFromSignal();
  return 0;
}

//**********************************************************************

ostream& KeepAllRoiBuildingService::
print(ostream& out, string prefix) const {
  out << prefix << "KeepAllRoiBuildingService:" << endl;
  out << prefix << "    LogLevel: " << m_LogLevel << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(KeepAllRoiBuildingService, AdcRoiBuildingService)

//**********************************************************************
