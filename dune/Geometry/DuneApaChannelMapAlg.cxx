////////////////////////////////////////////////////////////////////////
/// \file  DuneApaChannelMapAlg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////

#include "dune/Geometry/DuneApaChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

using std::string;
using std::vector;
using std::find;
using raw::ChannelID_t;
using geo::DuneApaChannelMapAlg;
using geo::PlaneID;
using geo::WireID;
using geo::View_t;
using geo::CryostatGeo;
using geo::TPCID;
using geo::PlaneGeo;
using readout::TPCsetID;
using readout::ROPID;
using geo::SigType_t;

typedef unsigned int Index;
Index badIndex = 999999;

#include <iostream>
using std::cout;
using std::endl;
//----------------------------------------------------------------------------

DuneApaChannelMapAlg::
DuneApaChannelMapAlg(const fhicl::ParameterSet& p)
: fSorter(nullptr) {
  fChannelsPerOpDet = p.get<unsigned int>("ChannelsPerOpDet");
  fOpDetFlag = 0;
  // If DetectorVersion is present, then check if this is 35t.
  string sdet;
  p.get_if_present<string>("DetectorVersion", sdet);
  if ( sdet.substr(0,7) == "dune35t" ) fOpDetFlag = 1;
}

//----------------------------------------------------------------------------

void DuneApaChannelMapAlg::setSorter(const geo::GeoObjectSorter& sort) {
  fSorter = &sort;
}

//----------------------------------------------------------------------------

void DuneApaChannelMapAlg::Initialize(GeometryData_t const& geodata) {
  Uninitialize();
    
  vector<CryostatGeo> const& crygeos = geodata.cryostats;
  Index ncry = crygeos.size();
  fNcryostat = ncry;
  if ( ncry == 0 ) {
    mf::LogError("DuneApaChannelMapAlg") << "No cryostats found.";
    return;
  }
  CryostatGeo const& crygeo = crygeos[0];
    
  mf::LogInfo("DuneApaChannelMapAlg") << "Initializing channel map...";
  fNTpc.resize(ncry);
  fNApa.resize(ncry);
  // The first channel array for each APA plane allow the cryostats to differ.
  fPlanesPerTpc.resize(ncry);
  fWiresPerPlane.resize(ncry);
  fApaTpcs.resize(ncry);
  fFirstChannelInThisPlane.resize(ncry);
  fFirstChannelInNextPlane.resize(ncry);
  fAnchoredPlaneRop.resize(ncry);
  fPlaneApa.resize(ncry);
  fPlaneRop.resize(ncry);
  fPlaneRopIndex.resize(ncry);
  fAnchoredWires.resize(ncry);
  fFirstChannelInThisRop.resize(ncry);
  fFirstChannelInNextRop.resize(ncry);
  fRopsPerApa.resize(ncry);
  fChannelsPerApa.resize(ncry);
  fChannelsPerRop.resize(ncry);
  fPlanesPerRop.resize(ncry);
  fRopTpc.resize(ncry);
  fRopPlane.resize(ncry);
  fNTpcMax = 0;
  fNApaMax = 0;
  fNRopMax = 0;
  for ( Index icry=0; icry<ncry; ++icry ) {
    Index ntpc = crygeos[icry].NTPC();
    Index napa = ntpc/2;  // Assume 1 APA for every two TPCs
    fNTpc[icry] = ntpc;
    fNApa[icry] = napa;
    if ( ntpc > fNTpcMax ) fNTpcMax = ntpc;
    if ( napa > fNApaMax ) fNApaMax = napa;
    fPlanesPerTpc[icry].resize(ntpc);
    fWiresPerPlane[icry].resize(ntpc);
    fFirstChannelInThisPlane[icry].resize(ntpc);
    fFirstChannelInNextPlane[icry].resize(ntpc);
    fAnchoredPlaneRop[icry].resize(ntpc);
    fPlaneApa[icry].resize(ntpc);
    fPlaneRop[icry].resize(ntpc);
    fPlaneRopIndex[icry].resize(ntpc);
    fAnchoredWires[icry].resize(ntpc);
    fRopsPerApa[icry].resize(napa, 4);
    fChannelsPerApa[icry].resize(napa, 0);
    fApaTpcs[icry].resize(napa);
    fChannelsPerRop[icry].resize(napa);
    fPlanesPerRop[icry].resize(napa);
    fFirstChannelInThisRop[icry].resize(napa);
    fFirstChannelInNextRop[icry].resize(napa);
    fRopTpc[icry].resize(napa);
    fRopPlane[icry].resize(napa);
    for ( Index itpc=0; itpc<ntpc; ++itpc ) {
      Index npla = crygeo.TPC(itpc).Nplanes();
      if ( npla != 3 ) {
        mf::LogError("DuneApaChannelMapAlg")
          << "# planes/TPC is " << npla << " rather than 3 " << " for TPC " << itpc;
        return;
      }
      fPlaneApa[icry][itpc].resize(npla, badIndex);
      fPlaneRop[icry][itpc].resize(npla, badIndex);
      fPlaneRopIndex[icry][itpc].resize(npla, badIndex);
      fFirstChannelInThisPlane[icry][itpc].resize(npla, 0);
      fFirstChannelInNextPlane[icry][itpc].resize(npla, 0);
      fAnchoredPlaneRop[icry][itpc].resize(npla, 0);
      fAnchoredWires[icry][itpc].resize(npla, 0);
      fWiresPerPlane[icry][itpc].resize(npla, 0);
      fPlanesPerTpc[icry][itpc] = npla;
      for ( Index ipla=0; ipla<npla; ++ipla ) {
        Index nwir = crygeo.TPC(itpc).Plane(ipla).Nwires();
        fWiresPerPlane[icry][itpc][ipla] = nwir;
      }
    }
    Index itpc = 0;
    for ( Index iapa=0; iapa<napa; ++iapa ) {
      Index nrop = fRopsPerApa[icry][iapa];
      if ( nrop > fNRopMax ) fNRopMax = nrop;
      fApaTpcs[icry][iapa].push_back(itpc);
      fApaTpcs[icry][iapa].push_back(itpc+1);
      fChannelsPerRop[icry][iapa].resize(nrop);
      fFirstChannelInThisRop[icry][iapa].resize(nrop, 0);
      fFirstChannelInNextRop[icry][iapa].resize(nrop, 0);
      fPlanesPerRop[icry][iapa].resize(nrop, 0);
      fRopTpc[icry][iapa].resize(nrop);
      fRopPlane[icry][iapa].resize(nrop);
      Index ipla = 0;
      // Induction planes
      for ( Index irop=0; irop<2; ++irop ) {
        fPlanesPerRop[icry][iapa][irop] = 2;
        fRopTpc[icry][iapa][irop].push_back(itpc);
        fRopTpc[icry][iapa][irop].push_back(itpc+1);
        fRopPlane[icry][iapa][irop].push_back(ipla);
        fRopPlane[icry][iapa][irop].push_back(ipla);
        ++ipla;
      }
      // Collection planes.
      for ( Index irop=2; irop<4; ++irop ) {
        fPlanesPerRop[icry][iapa][irop] = 1;
        fRopTpc[icry][iapa][irop].push_back(itpc);
        fRopPlane[icry][iapa][irop].push_back(ipla);
        ++itpc;
      }
    }
  }
  fPlaneIDs.clear();
  fTopChannel = 0;

  // Map planes to ROPs and wires to channels.
  // We use the geometry to deduce the numbers of channels/plane.
  // For induction planes, we look for the first adjacent pair of wires
  // with the same center z-position and assume the index of the first
  // of those gives the count of wires read out from one side.
  // These are called anchored wires.
  ChannelID_t icha = 0;
  for ( Index icry=0; icry!=ncry; ++icry ) {
    Index napa = fNApa[icry];
    for ( Index iapa=0; iapa!=napa; ++iapa ) {
      Index nrop = fRopsPerApa[icry][iapa];
      for ( Index irop=0; irop!=nrop; ++irop ) {
        Index nrpl = fPlanesPerRop[icry][iapa][irop];
        fFirstChannelInThisRop[icry][iapa][irop] = icha;
        for ( Index irpl=0; irpl!=nrpl; ++irpl ) {
          Index itpc = fRopTpc[icry][iapa][irop][irpl];
          Index ipla = fRopPlane[icry][iapa][irop][irpl];
          fPlaneApa[icry][itpc][ipla] = iapa;
          fPlaneRop[icry][itpc][ipla] = irop;
          fPlaneRopIndex[icry][itpc][ipla] = irpl;
          const PlaneGeo& plageo = crygeos[icry].TPC(itpc).Plane(ipla);
#if 0
const PlaneGeo plageo2 = plageo;
#endif
          View_t view = plageo.View();
          const Vector<View_t> eview = {geo::kU, geo::kV, geo::kZ};
          if ( view != eview[ipla] ) {
            throw cet::exception("DuneApaChannelMapAlg") << __func__
                  << ": View " << view << " is not the expected " << eview[ipla];
          }
          double xyz[3] = {0.};
          double xyz_next[3] = {0.};
          Index nAnchoredWires = 0;  // # wires from this TPC plane contributing to the ROP
          // Collection plane.
          Index nwir = fWiresPerPlane[icry][itpc][ipla];
          if ( view == geo::kZ ) {
            nAnchoredWires = nwir;
          // Induction planes.
          } else {
            for ( unsigned int iwir=0; iwir+1<nwir; ++iwir ) {
              plageo.Wire(iwir).GetCenter(xyz);
              plageo.Wire(iwir+1).GetCenter(xyz_next);
              if ( xyz[2] == xyz_next[2] ) {
                nAnchoredWires = iwir;
                break;
              }
            }
          }
          fAnchoredWires[icry][itpc][ipla] = nAnchoredWires;
          fFirstChannelInThisPlane[icry][itpc][ipla] = icha;
          icha += nAnchoredWires;
          fFirstChannelInNextPlane[icry][itpc][ipla] = icha;
          fAnchoredPlaneRop[icry][itpc][ipla] = irop;
        }
        fFirstChannelInNextRop[icry][iapa][irop] = icha;
      }
    }
  }
  fNchannels = icha;
  fChannelsPerAPA = fNchannels/fNApa[0];

  // Assign first channels for the TPCs.
  fPlaneData.resize(ncry);
  for ( Index icry=0; icry<ncry; ++icry ) {
    Index ntpc = fNTpc[icry];
    fPlaneData[icry].resize(ntpc);
    for ( Index itpc=0; itpc<ntpc; ++itpc ) {
      fPlaneData[icry][itpc].resize(crygeos[icry].TPC(itpc).Nplanes());
      for ( Index ipla=0; ipla<crygeos[icry].TPC(itpc).Nplanes(); ++ipla ) {
        PlaneData_t& PlaneData = fPlaneData[icry][itpc][ipla];
        const PlaneGeo& thePlane = crygeos[icry].TPC(itpc).Plane(ipla);
        double xyz[3];
        fPlaneIDs.emplace(icry, itpc, ipla);
        thePlane.Wire(0).GetCenter(xyz);
        PlaneData.fFirstWireCenterY = xyz[1];
        PlaneData.fFirstWireCenterZ = xyz[2];
        // we are interested in the ordering of wire numbers: we find that a
        // point is N wires left of a wire W: is that wire W + N or W - N?
        // In fact, for TPC #0 it is W + N for V and Z planes, W - N for U
        // plane; for TPC #0 it is W + N for V and Z planes, W - N for U
        PlaneData.fWireSortingInZ = thePlane.WireIDincreasesWithZ()? +1.: -1.;

	  // find boundaries of the outside APAs for this plane by looking at endpoints of wires

	  double endpoint[3];
	  thePlane.Wire(0).GetStart(endpoint);
	  PlaneData.fYmax = endpoint[1];
	  PlaneData.fYmin = endpoint[1];
	  PlaneData.fZmax = endpoint[2];
	  PlaneData.fZmin = endpoint[2];
	  unsigned int nwires = thePlane.Nwires(); 
	  for (unsigned int iwire=0;iwire<nwires;iwire++){
  	    thePlane.Wire(iwire).GetStart(endpoint);
	    PlaneData.fYmax = std::max(PlaneData.fYmax,endpoint[1]);
	    PlaneData.fYmin = std::min(PlaneData.fYmin,endpoint[1]);
	    PlaneData.fZmax = std::max(PlaneData.fZmax,endpoint[2]);
	    PlaneData.fZmin = std::min(PlaneData.fZmin,endpoint[2]);
  	    thePlane.Wire(iwire).GetEnd(endpoint);
	    PlaneData.fYmax = std::max(PlaneData.fYmax,endpoint[1]);
	    PlaneData.fYmin = std::min(PlaneData.fYmin,endpoint[1]);
	    PlaneData.fZmax = std::max(PlaneData.fZmax,endpoint[2]);
	    PlaneData.fZmin = std::min(PlaneData.fZmin,endpoint[2]);	    
	  } // loop on wire 

      } // for plane
    } // for TPC
  } // for cryostat

  Index npla = crygeo.TPC(0).Nplanes();
  fWirePitch.resize(npla);
  fOrientation.resize(npla);
  fSinOrientation.resize(npla);
  fCosOrientation.resize(npla);

  for ( Index ipla=0; ipla<npla; ++ipla ) {
    fWirePitch[ipla] = crygeo.TPC(0).WirePitch(ipla);
    fOrientation[ipla] = crygeo.TPC(0).Plane(ipla).Wire(0).ThetaZ();
    fSinOrientation[ipla] = sin(fOrientation[ipla]);
    fCosOrientation[ipla] = cos(fOrientation[ipla]);
  }

  for ( Index icry=0; icry<ncry; ++icry ) {
    mf::LogVerbatim("DuneApaChannelMapAlg") << "Cryostat " << icry << ":"; 
    mf::LogVerbatim("DuneApaChannelMapAlg") << "  " << fNchannels << " total channels"; 
    mf::LogVerbatim("DuneApaChannelMapAlg") << "  " << fNTpc[icry]/2 << " APAs";       
    mf::LogVerbatim("DuneApaChannelMapAlg") << "  For all identical APA:" ; 
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Number of channels per APA = " << fChannelsPerAPA ; 
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    U channels per APA = " << 2*fAnchoredWires[0][0][0] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    V channels per APA = " << 2*fAnchoredWires[0][0][1] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Z channels per APA = " << 2*fAnchoredWires[0][0][2] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in U Plane = " << fWirePitch[0] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in V Plane = " << fWirePitch[1] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in Z Plane = " << fWirePitch[2] ;
  }

}

//----------------------------------------------------------------------------

void DuneApaChannelMapAlg::Uninitialize() {
  PlaneInfoMap_t<ChannelID_t>().swap(fFirstChannelInThisPlane);
  PlaneInfoMap_t<ChannelID_t>().swap(fFirstChannelInNextPlane);
}

//----------------------------------------------------------------------------

vector<WireID> DuneApaChannelMapAlg::ChannelToWire(ChannelID_t icha) const {
  vector< WireID > wirids;
  if ( icha >= fNchannels ) return wirids;
  // Loop over ROPs to find the one holding this channel.
  Index ncry = fNcryostat;
  Index icry = badIndex;
  Index iapa = badIndex;
  Index irop = badIndex;
  Index ichaRop = badIndex;    // Channel number in the ROP
  bool found = false;
  for ( icry=0; icry<ncry; ++icry ) {
    Index napa = fNApa[icry];
    for ( iapa=0; iapa<napa; ++iapa ) {
      Index nrop = fRopsPerApa[icry][iapa];
      for ( irop=0; irop<nrop; ++irop ) {
        Index icha1 = fFirstChannelInThisRop[icry][iapa][irop];
        Index icha2 = fFirstChannelInNextRop[icry][iapa][irop];
        found = icha >= icha1 && icha < icha2;
        if ( found ) {
          ichaRop = icha - icha1;
          break;
        }
      }
      if ( found ) break;
    }
    if ( found ) break;
  }
  if ( icry >= ncry ) {
    mf::LogError("DuneApaChannelMapAlg") << "Unable to find APA plane for channel " << icha;
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": Unable to find APA plane for channel " << icha;
    return wirids;
  }
  // Extract TPC(s) from ROP
  Index nrpl = fPlanesPerRop[icry][iapa][irop];
  if ( nrpl == 0 ) throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": No TPC planes.";
  if ( nrpl > 2 ) throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": Too many TPC planes.";
  Index itpc1 = fRopTpc[icry][iapa][irop][0];
  Index ipla = fRopPlane[icry][iapa][irop][0];
  Index itpc2 = (nrpl > 1 ) ? fRopTpc[icry][iapa][irop][1] : itpc1;
  Index nAnchored = fAnchoredWires[icry][itpc1][ipla];
  bool wrapped = ipla < 2;
  if ( wrapped && itpc2 == itpc1 )
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": 2nd plane not found for wrapped ROP";
  if ( wrapped && ipla != fRopPlane[icry][iapa][irop][1] )
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": Wrapped planes have inconsistent indices.";
  // For now, assume the second TPC plane has the same # anchored wires.
  // Code will need some work if we want to relax this assumption.
  if (  wrapped && fAnchoredWires[icry][itpc2][ipla] != nAnchored )
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": Planes have inconsistent anchor counts.";
  // If this is a wrapped ROP and the wire number is larger than nAnchored, then
  // the first wire for this channel is in the other TPC plane.
  Index itpc = itpc1;
  Index iwir = ichaRop;
  if ( wrapped && iwir >= nAnchored ) {
    itpc = itpc2;
    iwir -= nAnchored;
  }
  if ( iwir >= nAnchored ) {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": Invalid channel: iwir =" << iwir;
  }
  // Loop over wires and create IDs.
  while ( iwir < fWiresPerPlane[icry][itpc][ipla] ) {
    WireID wirid(icry, itpc, ipla, iwir);
    wirids.push_back(wirid);
    iwir += fAnchoredWires[icry][itpc][ipla];
    itpc = (itpc == itpc1) ? itpc2 : itpc1;
  }
  return wirids;
}

//----------------------------------------------------------------------------

unsigned int DuneApaChannelMapAlg::Nchannels() const {
   return fNchannels;
}
  
//----------------------------------------------------------------------------

unsigned int DuneApaChannelMapAlg::Nchannels(ROPID const& ropid) const {
  Index icry = ropid.Cryostat;
  Index iapa = ropid.TPCset;
  Index irop = ropid.ROP;
  Index icha1 = fFirstChannelInThisRop[icry][iapa][irop];
  Index icha2 = fFirstChannelInNextRop[icry][iapa][irop];
  return icha2 - icha1;
}

//----------------------------------------------------------------------------
  
double DuneApaChannelMapAlg::
WireCoordinate(double YPos, double ZPos, unsigned int PlaneNo, unsigned int TPCNo,
               unsigned int cstat) const {
  return WireCoordinate(YPos, ZPos, PlaneID(cstat, TPCNo, PlaneNo));
}
  
//----------------------------------------------------------------------------

double DuneApaChannelMapAlg::
WireCoordinate(double YPos, double ZPos, PlaneID const& plaid) const {
  const PlaneData_t& PlaneData = AccessElement(fPlaneData, plaid);
  Index icry = plaid.Cryostat;
  Index itpc = plaid.TPC;
  Index ipla = plaid.Plane;
  // The formula used here is geometric:
  // distance = delta_y cos(theta_z) + delta_z sin(theta_z)
  // with a correction for the orientation of the TPC:
  // odd TPCs have supplementary wire angle (pi-theta_z), changing cosine sign
  double backsign = (fPlaneRopIndex[icry][itpc][ipla] == 1) ? -1.0 : 1.0;
  float distance =
    -(YPos - PlaneData.fFirstWireCenterY) * backsign * fCosOrientation[ipla]
    + (ZPos - PlaneData.fFirstWireCenterZ) * fSinOrientation[ipla]
    ;
  // The sign of this formula is correct if the wire with larger ID is on the
  // "right" (intuitively, larger z; rigorously, smaller intercept)
  // than this one.
  // Of course, we are not always that lucky. fWireSortingInZ fixes our luck.
  return PlaneData.fWireSortingInZ * distance/fWirePitch[ipla];
}
  
//----------------------------------------------------------------------------

WireID DuneApaChannelMapAlg::
NearestWireID(const TVector3& xyz, PlaneID const& plaid) const {


    // cap the position to be within the boundaries of the wire endpoints.
    // This simulates charge drifting in from outside of the wire frames inwards towards
    // the first and last collection wires on the side, and towards the right endpoints

  const PlaneData_t& PlaneData = AccessElement(fPlaneData, plaid);
  double ycap = std::max(PlaneData.fYmin,std::min(PlaneData.fYmax,xyz.Y()));
  double zcap = std::max(PlaneData.fZmin,std::min(PlaneData.fZmax,xyz.Z()));

  int iwirSigned = 0.5 + WireCoordinate(ycap, zcap, plaid);
  Index iwir = (iwirSigned < 0) ? 0 : iwirSigned;
  Index icry = plaid.Cryostat;
  Index itpc = plaid.TPC;
  Index ipla = plaid.Plane;
  Index iwirMax = fWiresPerPlane[icry][itpc][ipla] - 1;
  if ( iwir > iwirMax ) iwir = iwirMax;
  return WireID(plaid, iwir);
}

//----------------------------------------------------------------------------

ChannelID_t DuneApaChannelMapAlg::PlaneWireToChannel(WireID const& wirid) const {
  Index icry = wirid.Cryostat;
  Index itpc = wirid.TPC;
  Index ipla = wirid.Plane;
  Index ichaRop = wirid.Wire;
  Index iapa = fPlaneApa[icry][itpc][ipla];
  Index irop = fPlaneRop[icry][itpc][ipla];
  Index irpl = fPlaneRopIndex[icry][itpc][ipla];
  Index nrpl = fPlanesPerRop[icry][iapa][irop];
  Index ncha = fAnchoredWires[icry][itpc][ipla];
  if ( nrpl > 1 ) {      // Wrapped ROP
    Index ipla1 = fRopPlane[icry][iapa][irop][0];
    Index ipla2 = fRopPlane[icry][iapa][irop][1];
    // Wire is in the back TPC.
    if ( irpl == 1 ) {
      ichaRop += fAnchoredWires[icry][itpc][ipla1];
      ncha += fAnchoredWires[icry][itpc][ipla1];
    // Wire is in the front TPC.
    } else {
      ncha += fAnchoredWires[icry][itpc][ipla2];
    }
  }
  // Channel # in ROP is modulus the # channels in the ROP.
  Index icha1 = fFirstChannelInThisRop[icry][iapa][irop];
  Index icha = icha1 + ichaRop%ncha;
  return icha;
}

//----------------------------------------------------------------------------

SigType_t DuneApaChannelMapAlg::SignalTypeForChannelImpl(ChannelID_t const icha) const {
  Index ncry = fNcryostat;
  for ( Index icry=0; icry<ncry; ++icry ) {
    Index napa=fNApa[icry];
    for ( Index iapa=0; iapa<napa; ++iapa ) {
      Index nrop = fRopsPerApa[icry][iapa];
      for ( Index irop=0; irop<nrop; ++irop ) {
        Index icha1 = fFirstChannelInThisRop[icry][iapa][irop];
        Index icha2 = fFirstChannelInNextRop[icry][iapa][irop];
        if ( icha >= icha1 && icha < icha2 ) {
          if      ( irop < 2 ) return geo::kInduction;
          else if ( irop < 4 ) return geo::kCollection;
          else return geo::kMysteryType;
        }
      }
    }
  }
  return geo::kMysteryType;
}

//----------------------------------------------------------------------------
/* // this code should be equivalent to the logic implemented in geo::PlaneGeo::UpdateView()
View_t DuneApaChannelMapAlg::View(ChannelID_t const icha) const {
  Index ncry = fNcryostat;
  for ( Index icry=0; icry<ncry; ++icry ) {
    Index napa=fNApa[icry];
    for ( Index iapa=0; iapa<napa; ++iapa ) {
      Index nrop = fRopsPerApa[icry][iapa];
      for ( Index irop=0; irop<nrop; ++irop ) {
        Index icha1 = fFirstChannelInThisRop[icry][iapa][irop];
        Index icha2 = fFirstChannelInNextRop[icry][iapa][irop];
        if ( icha >= icha1 && icha < icha2 ) {
          if ( irop == 0 ) return geo::kU;
          if ( irop == 1 ) return geo::kV;
          if ( irop == 2 ) return geo::kZ;
          if ( irop == 3 ) return geo::kZ;
          return geo::kUnknown;
        }
      }
    }
  }
  return geo::kUnknown;
}
*/

//----------------------------------------------------------------------------

std::set<PlaneID> const& DuneApaChannelMapAlg::PlaneIDs() const {
  return fPlaneIDs;
}

//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::NOpChannels(unsigned int NOpDets) const {
  if ( fOpDetFlag == 1 ) return 12*NOpDets;
  return fChannelsPerOpDet*NOpDets;
}

//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::NOpHardwareChannels(unsigned int opDet) const {
  if ( fOpDetFlag == 1 ) {
    // CSU 3-sipm design
    if (opDet == 0 || opDet == 4 || opDet == 6) return 8;
    // LSU 2-sipm design
    if (opDet == 2) return 2;
    // IU 12-sipm design
    return 12;
  }
  return fChannelsPerOpDet;
}

//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::OpChannel(unsigned int detNum, unsigned int channel) const {
  if ( fOpDetFlag == 1 ) return (detNum * 12) + channel;
  unsigned int uniqueChannel = (detNum * fChannelsPerOpDet) + channel;
  return uniqueChannel;
}

//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::OpDetFromOpChannel(unsigned int opChannel) const {
  if ( fOpDetFlag == 1 ) return opChannel/12;
  unsigned int detectorNum = opChannel/fChannelsPerOpDet;
  return detectorNum;
}

//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::HardwareChannelFromOpChannel(unsigned int opChannel) const {
  if ( fOpDetFlag == 1 ) return opChannel%12;
  unsigned int channel = opChannel % fChannelsPerOpDet;
  return channel;
}

//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::NTPCsets(CryostatID const& cryid) const {
  if (!HasCryostat(cryid)) return 0;
  Index icry = cryid.Cryostat;
  return fNApa[icry];
}

//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::MaxTPCsets() const {
  return fNApaMax;
}

//----------------------------------------------------------------------------

bool DuneApaChannelMapAlg::HasTPCset(TPCsetID const& tpcsetid) const {
  return HasCryostat(tpcsetid) && (tpcsetid.TPCset < NTPCsets(tpcsetid));
}
  
//----------------------------------------------------------------------------

TPCsetID DuneApaChannelMapAlg::TPCtoTPCset(TPCID const& tpcid) const {
  Index icry = tpcid.Cryostat;
  Index itpc = tpcid.TPC;
  Index ipla = 0;
  Index iapa = fPlaneApa[icry][itpc][ipla];
  return TPCsetID(icry, iapa);
}

//----------------------------------------------------------------------------

vector<TPCID> DuneApaChannelMapAlg::TPCsetToTPCs(TPCsetID const& apaid) const {
  Index icry = apaid.Cryostat;
  Index iapa = apaid.TPCset;
  vector<TPCID> tpcids;
  for ( Index itpc : fApaTpcs[icry][iapa] ) tpcids.push_back(TPCID(icry, itpc));
  return tpcids;
}

//----------------------------------------------------------------------------

TPCID DuneApaChannelMapAlg::FirstTPCinTPCset (TPCsetID const& apaid) const {
  Index icry = apaid.Cryostat;
  Index iapa = apaid.TPCset;
  Index itpc = TPCID::InvalidID;
  if ( fApaTpcs[icry][iapa].size() ) itpc = fApaTpcs[icry][iapa][0];
  return TPCID(icry, itpc);
}

//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::NROPs(TPCsetID const& apaid) const {
  if (!HasTPCset(apaid)) return 0;
  Index icry = apaid.Cryostat;
  Index iapa = apaid.TPCset;
  return fRopsPerApa[icry][iapa];
}
  
//----------------------------------------------------------------------------

Index DuneApaChannelMapAlg::MaxROPs() const {
  return fNRopMax;
}
  
//----------------------------------------------------------------------------

bool DuneApaChannelMapAlg::HasROP(ROPID const& ropid) const {
  return HasTPCset(ropid) && (ropid.ROP < NROPs(ropid));
}

//----------------------------------------------------------------------------

ROPID DuneApaChannelMapAlg::WirePlaneToROP(PlaneID const& plaid) const {
  Index icry = plaid.Cryostat;
  Index itpc = plaid.TPC;
  Index ipla = plaid.Plane;
  Index iapa = fPlaneApa[icry][itpc][ipla];
  Index irop = fPlaneRop[icry][itpc][ipla];
  return ROPID(icry, iapa, irop);
}

//----------------------------------------------------------------------------

vector<PlaneID> DuneApaChannelMapAlg::ROPtoWirePlanes(ROPID const& ropid) const {
  vector<PlaneID> plaids;
  Index icry = ropid.Cryostat;
  Index iapa = ropid.TPCset;
  Index irop = ropid.ROP;
  Index nrpl = fPlanesPerRop[icry][iapa][irop];
  for ( Index irpl=0; irpl<nrpl; ++irpl ) {
    Index itpc = fRopTpc[icry][iapa][irop][irpl];
    Index ipla = fRopPlane[icry][iapa][irop][irpl];
    plaids.push_back(PlaneID(icry, itpc, ipla));
  }
  return plaids;
}

//----------------------------------------------------------------------------

vector<TPCID> DuneApaChannelMapAlg::ROPtoTPCs(ROPID const& ropid) const {
  vector<TPCID> tpcids;
  vector<PlaneID> plaids = ROPtoWirePlanes(ropid);
  for ( const PlaneID& plaid : plaids ) {
    TPCID tpcid = plaid;
    if ( find(tpcids.begin(), tpcids.end(), tpcid) == tpcids.end() ) tpcids.push_back(tpcid);
  }
  return tpcids;
}

//----------------------------------------------------------------------------

ROPID DuneApaChannelMapAlg::ChannelToROP(ChannelID_t icha) const {
  Index ncry = fNcryostat;
  for ( Index icry=0; icry<ncry; ++icry ) {
    Index napa = fNApa[icry];
    for ( Index iapa=0; iapa<napa; ++iapa ) {
      Index nrop = fRopsPerApa[icry][iapa];
      for ( Index irop=0; irop<nrop; ++irop ) {
        Index icha1 = fFirstChannelInThisRop[icry][iapa][irop];
        Index icha2 = fFirstChannelInNextRop[icry][iapa][irop];
        if ( icha >= icha1 && icha < icha2 ) return ROPID(icry, iapa, irop);
      }
    }
  }
  return(ROPID(CryostatID::InvalidID, TPCsetID::InvalidID, ROPID::InvalidID));
}

//----------------------------------------------------------------------------
geo::GeoObjectSorter const& DuneApaChannelMapAlg::Sorter() const {
  if ( fSorter == nullptr )
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": Sorter is missing.";
  return *fSorter;
} // DuneApaChannelMapAlg::Sorter()
  
//----------------------------------------------------------------------------

ChannelID_t DuneApaChannelMapAlg::FirstChannelInROP(ROPID const& ropid) const {
  Index icry = ropid.Cryostat;
  Index iapa = ropid.TPCset;
  Index irop = ropid.ROP;
  Index icha1 = fFirstChannelInThisRop[icry][iapa][irop];
  return icha1;
}

//----------------------------------------------------------------------------

PlaneID DuneApaChannelMapAlg::FirstWirePlaneInROP(ROPID const& ropid) const {
  vector<PlaneID> plaids = ROPtoWirePlanes(ropid);
  if ( plaids.size() ) return plaids[0];
  return PlaneID(CryostatID::InvalidID, TPCID::InvalidID, PlaneID::InvalidID);
}

//----------------------------------------------------------------------------
