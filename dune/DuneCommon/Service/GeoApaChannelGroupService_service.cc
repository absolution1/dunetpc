// GeoApaChannelGroupService.cxx

#include "dune/DuneCommon/Service/GeoApaChannelGroupService.h"
#include <sstream>
#include <iomanip>
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/ChannelMapAlg.h"

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
: m_size(0) {
  // Find the total # APAs and assign width for APA index.
  Index napatot = 0;
  for ( const geo::CryostatID& cryid : m_pgeo->IterateCryostatIDs() ) {
    napatot += m_pgeo->NTPCsets(cryid);
  }
  unsigned int w = log10(napatot-1) + 1;
  m_names.resize(napatot);
  m_chanvecs.resize(napatot);
  // Loop over cryostats.
  for ( const geo::CryostatID& cryid : m_pgeo->IterateCryostatIDs() ) {
    Index napa = m_pgeo->NTPCsets(cryid);
    // Loop over APAs.
    for ( Index iapa=0; iapa<napa; ++iapa ) {
      APAID apaid(cryid, iapa);
      Index nrop = m_pgeo->NROPs(apaid);
      // Build name.
      ostringstream ssname;
      ssname << "apa" << setfill('0') << setw(w) << iapa;
      m_names[iapa] = ssname.str();
      // Loop over ROPs and fill channel vector.
      for ( Index irop=0; irop<nrop; ++irop ) {
        ROPID ropid(apaid, irop);
        // Find channels.
        Index icha1 = m_pgeo->FirstChannelInROP(ropid); 
        Index icha2 = icha1 + m_pgeo->Nchannels(ropid);
        for ( unsigned int icha=icha1; icha<icha2; ++icha ) m_chanvecs[iapa].push_back(icha);
      }
    }
  }
  m_size = napatot;
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
