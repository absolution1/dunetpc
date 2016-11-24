// FixedChannelGroupService.cxx

#include "dune/DuneCommon/Service/FixedChannelGroupService.h"
#include <sstream>
#include <iomanip>

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;
using std::setfill;
using std::setw;
using std::vector;
using std::move;

typedef ChannelGroupService::Index Index;
typedef ChannelGroupService::Name Name;
typedef ChannelGroupService::ChannelVector ChannelVector;

//**********************************************************************

FixedChannelGroupService::
FixedChannelGroupService(fhicl::ParameterSet const& pset) {
  vector<string> names = pset.get_names();
  for ( string name : names ) {
    if ( name == "service_provider" ) continue;
    if ( name == "service_type" ) continue;
    ChannelVector chans = pset.get<ChannelVector>(name);
    m_names.push_back(name);
    m_chanvecs.emplace_back(move(chans));
  }
}

//**********************************************************************

FixedChannelGroupService::
FixedChannelGroupService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: FixedChannelGroupService(pset) { }

//**********************************************************************

Index FixedChannelGroupService::size() const {
  return m_names.size();
}

//**********************************************************************

Name FixedChannelGroupService::name(Index igrp) const {
  if ( igrp < m_names.size() ) return m_names[igrp];
  return "NoSuchGroup";
}

//**********************************************************************

const ChannelVector& FixedChannelGroupService::channels(Index igrp) const {
  if ( igrp < m_chanvecs.size() ) return m_chanvecs[igrp];
  static ChannelVector empty;
  return empty;
}

//**********************************************************************

ostream& FixedChannelGroupService::print(ostream& out, string prefix) const {
  out << prefix << "FixedChannelGroupService: " << endl;
  out << prefix << "  # group: " << size() << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(FixedChannelGroupService, ChannelGroupService)

//**********************************************************************
