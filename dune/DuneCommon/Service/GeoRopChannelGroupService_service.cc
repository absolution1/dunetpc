// GeoRopChannelGroupService.cxx

#include "dune/DuneCommon/Service/GeoRopChannelGroupService.h"
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

GeoRopChannelGroupService::
GeoRopChannelGroupService(fhicl::ParameterSet const& pset)
: m_size(0) {
  // Find the total # APAs and assign width for APA index.
  Index napatot = 0;
  for ( const geo::CryostatID& cryid : m_pgeo->IterateCryostatIDs() ) {
    napatot += m_pgeo->NTPCsets(cryid);
  }
  unsigned int w = log10(napatot-1) + 1;
  // Loop over cryostats.
  Index krop = 0;
  for ( const geo::CryostatID& cryid : m_pgeo->IterateCryostatIDs() ) {
    Index napa = m_pgeo->NTPCsets(cryid);
    // Loop over APAs.
    for ( Index iapa=0; iapa<napa; ++iapa ) {
      APAID apaid(cryid, iapa);
      Index nrop = m_pgeo->NROPs(apaid);
      // Loop over ROPs.
      unsigned int nvu = 0;
      unsigned int nvv = 0;
      unsigned int nvz = 0;
      // First pass to get the number of each view.
      for ( Index irop=0; irop<nrop; ++irop ) {
        ROPID ropid(apaid, irop);
        geo::View_t view = m_pgeo->View(ropid);
        if      ( view == geo::kU ) ++nvu;
        else if ( view == geo::kV ) ++nvv;
        else if ( view == geo::kZ ) ++nvz;
      }
      unsigned int ivu = 0;
      unsigned int ivv = 0;
      unsigned int ivz = 0;
      // Loop over ROPs, assign name and fill channel vector.
      for ( Index irop=0; irop<nrop; ++irop ) {
        ROPID ropid(apaid, irop);
        m_names.resize(krop+1);
        m_chanvecs.resize(krop+1);
        // Build name.
        ostringstream ssname;
        ssname << "apa" << setfill('0') << setw(w) << iapa;
        geo::View_t view = m_pgeo->View(ropid);
        if        ( view == geo::kU ) {
          ssname << "u";
          ++ivu;
          if ( nvu > 1 ) ssname << ivu;
        } else if ( view == geo::kV ) {
          ssname << "v";
          ++ivv;
          if ( nvv > 1 ) ssname << ivv;
        } else if ( view == geo::kZ ) {
          ssname << "z";
          ++ivz;
          if ( nvz > 1 ) ssname << ivz;
        } else {
          ssname << "UnknownView";
        }
        m_names[krop] = ssname.str();
        // Find channels.
        Index icha1 = m_pgeo->FirstChannelInROP(ropid); 
        Index icha2 = icha1 + m_pgeo->Nchannels(ropid);
        for ( unsigned int icha=icha1; icha<icha2; ++icha ) m_chanvecs[krop].push_back(icha);
        ++krop;
      }
    }
  }
  m_size = krop;
}

//**********************************************************************

GeoRopChannelGroupService::
GeoRopChannelGroupService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: GeoRopChannelGroupService(pset) { }

//**********************************************************************

Index GeoRopChannelGroupService::size() const {
  return m_size;
}

//**********************************************************************

Name GeoRopChannelGroupService::name(Index igrp) const {
  if ( igrp < m_names.size() ) return m_names[igrp];
  return "NoSuchRop";
}

//**********************************************************************

const ChannelVector& GeoRopChannelGroupService::channels(Index igrp) const {
  if ( igrp < m_chanvecs.size() ) return m_chanvecs[igrp];
  static ChannelVector empty;
  return empty;
}

//**********************************************************************

ostream& GeoRopChannelGroupService::print(ostream& out, string prefix) const {
  out << prefix << "GeoRopChannelGroupService: " << endl;
  out << prefix << "  # ROP: " << size() << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(GeoRopChannelGroupService, ChannelGroupService)

//**********************************************************************
