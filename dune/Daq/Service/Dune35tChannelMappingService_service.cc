// Dune35tChannelMappingService.cxx

#include "dune/Daq/Service/Dune35tChannelMappingService.h"
#include <limits>
#include <fstream>

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ifstream;

#undef UseSeedService

typedef ChannelMappingService::Channel Channel;

//**********************************************************************

Dune35tChannelMappingService::
Dune35tChannelMappingService(fhicl::ParameterSet const& pset)
: m_LogLevel(1) {
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  if ( m_LogLevel > 0 ) {
    print(cout, "Dune35tChannelMappingService::ctor: ");
  }
}

//**********************************************************************

Dune35tChannelMappingService::
Dune35tChannelMappingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: Dune35tChannelMappingService(pset) { }

//**********************************************************************

Channel Dune35tChannelMappingService::offline(Channel chin) const {
  return hlbnesvc->Offline(chin);
}

//**********************************************************************

Channel Dune35tChannelMappingService::online(Channel chin) const {
  return hlbnesvc->Online(chin);
}

//**********************************************************************

ostream& Dune35tChannelMappingService::print(ostream& out, string prefix) const {
  string myname = prefix + "  ";
  cout << myname << "   LogLevel: " << m_LogLevel << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(Dune35tChannelMappingService, ChannelMappingService)

//**********************************************************************
