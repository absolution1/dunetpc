// GeoApaChannelGroupService.cxx

#include "dune/DuneCommon/Service/GeoApaChannelGroupService.h"
#include <sstream>
#include <iomanip>
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;
using std::setfill;
using std::setw;
using std::vector;

typedef ChannelGroupService::Index Index;
typedef ChannelGroupService::Name Name;
typedef ChannelGroupService::ChannelVector ChannelVector;

typedef readout::TPCsetID APAID;
typedef readout::ROPID ROPID;

//**********************************************************************

GeoApaChannelGroupService::
GeoApaChannelGroupService(fhicl::ParameterSet const& pset)
: m_LogLevel(pset.get<int>("LogLevel")),
  m_ApaRops(pset.get<vector<int>>("ApaRops")),
  m_size(0) {
  const string myname = "GeoApaChannelGroupService::ctor: ";
  // Find the total # APAs and assign width for APA index.
  Index ngrp = 0;
  for ( const geo::CryostatID& cryid : m_pgeo->IterateCryostatIDs() ) {
    ngrp += m_pgeo->NTPCsets(cryid);
  }
  unsigned int w = log10(ngrp-1) + 1;
  m_names.resize(ngrp);
  m_chanvecs.resize(ngrp);
  Index igrp = 0;
  Index ncry = m_pgeo->Ncryostats();
  // Loop over cryostats.
  for ( const geo::CryostatID& cryid : m_pgeo->IterateCryostatIDs() ) {
    Index napa = m_pgeo->NTPCsets(cryid);
    Index icry = cryid.Cryostat;
    // Loop over APAs.
    for ( Index iapa=0; iapa<napa; ++iapa, ++igrp ) {
      APAID apaid(cryid, iapa);
      Index nrop = m_pgeo->NROPs(apaid);
      // Build name.
      ostringstream ssname;
      if ( ncry > 1 ) ssname << "cry" << icry << "-";
      ssname << "apa" << setfill('0') << setw(w) << iapa;
      m_names[igrp] = ssname.str();
      // Loop over ROPs and fill channel vector.
      vector<bool> skipRop(nrop, false);
      Index naparop = m_ApaRops.size();
      if ( naparop > 0 ) {
        Index pat = m_ApaRops[iapa%naparop];
        for ( Index irop=0; irop<nrop; ++irop ) {
          Index mask = 1 << irop;
          bool keep = mask & pat;
          if ( ! keep ) skipRop[irop] = true;
        }
      }
      Index nropUsed = 0;
      for ( Index irop=0; irop<nrop; ++irop ) {
        if ( skipRop[irop] ) continue;
        ROPID ropid(apaid, irop);
        // Find channels.
        Index icha1 = m_pgeo->FirstChannelInROP(ropid); 
        Index icha2 = icha1 + m_pgeo->Nchannels(ropid);
        for ( unsigned int icha=icha1; icha<icha2; ++icha ) m_chanvecs[igrp].push_back(icha);
        ++nropUsed;
      }
      if ( m_LogLevel >= 1 ) {
        cout << myname << "Group " << igrp << " (" << m_names[igrp] << ") has "
             << m_chanvecs[igrp].size() << " channels from "
             << nropUsed << "/" << nrop << " readout planes." << endl;
      }
    }
  }
  m_size = ngrp;
}

//**********************************************************************

GeoApaChannelGroupService::
GeoApaChannelGroupService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: GeoApaChannelGroupService(pset) { }

//**********************************************************************

Index GeoApaChannelGroupService::size() const {
  return m_size;
}

//**********************************************************************

Name GeoApaChannelGroupService::name(Index igrp) const {
  if ( igrp < m_names.size() ) return m_names[igrp];
  return "NoSuchApa";
}

//**********************************************************************

const ChannelVector& GeoApaChannelGroupService::channels(Index igrp) const {
  if ( igrp < m_chanvecs.size() ) return m_chanvecs[igrp];
  static ChannelVector empty;
  return empty;
}

//**********************************************************************

ostream& GeoApaChannelGroupService::print(ostream& out, string prefix) const {
  out << prefix << "GeoApaChannelGroupService: " << endl;
  out << prefix << "  # ROP: " << size() << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(GeoApaChannelGroupService, ChannelGroupService)

//**********************************************************************
