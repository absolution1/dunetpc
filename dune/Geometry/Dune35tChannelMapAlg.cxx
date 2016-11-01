// Dune35tChannelMapAlg.cxx

#include "dune/Geometry/Dune35tChannelMapAlg.h"
#include <iomanip>
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

using std::setw;
using geo::CryostatID;
using geo::TPCID;
using geo::PlaneID;
using readout::TPCsetID;
using readout::ROPID;

#include <iostream>
using std::cout;
using std::endl;

namespace geo {

//----------------------------------------------------------------------------

Dune35tChannelMapAlg::Dune35tChannelMapAlg(fhicl::ParameterSet const& p)
: fSorter(geo::GeoObjectSorter35(p)) { }

//----------------------------------------------------------------------------

void Dune35tChannelMapAlg::Initialize( GeometryData_t& geodata ) {

  Uninitialize();
  
  std::vector<geo::CryostatGeo*>& cgeo = geodata.cryostats;
  std::vector<geo::AuxDetGeo*>  & adgeo = geodata.auxDets;
  
  fNcryostat = cgeo.size();
  
  mf::LogInfo("Dune35tChannelMapAlg") << "Sorting...";

  fSorter.SortAuxDets(adgeo);
  fSorter.SortCryostats(cgeo);
  for(size_t c = 0; c < cgeo.size(); ++c) 
    cgeo[c]->SortSubVolumes(fSorter);


  mf::LogInfo("Dune35tChannelMapAlg") << "Initializing...";
    
  fNTPC.resize(fNcryostat);
  fNAPA.resize(fNcryostat);
  fWiresPerRop.resize(fNcryostat);
  fFirstChannelInNextRop.resize(fNcryostat);
  fFirstChannelInThisRop.resize(fNcryostat);
  nAnchoredWires.resize(fNcryostat);
  fViews.clear();
  fPlaneIDs.clear();
  fPlanesPerTPC = cgeo[0]->TPC(0).Nplanes();
  fPlanesPerAPA = fPlanesPerTPC + 1;  // Assume induction are shared but collection is not.
  unsigned int tpcPlaneForApaPlane[4] = {0, 1, 2, 2};
  unsigned int tpcOffsetForApaPlane[4] = {0, 0, 0, 1};

  fTopChannel = 0;

  // Size some vectors and initialize the FirstChannel vectors.
  for ( unsigned int cs = 0; cs != fNcryostat; ++cs ) {
    
    fNTPC[cs] = cgeo[cs]->NTPC();
    fNAPA[cs] = fNTPC[cs]/2;

    nAnchoredWires[cs].resize(fNTPC[cs]);
    fWiresPerRop[cs].resize(fNTPC[cs]);
    fFirstChannelInThisRop[cs].resize(fNAPA[cs]);
    fFirstChannelInNextRop[cs].resize(fNAPA[cs]);

    for ( unsigned int apa = 0; apa != fNAPA[cs]; ++apa ) {
      
      nAnchoredWires[cs][apa].resize(fPlanesPerAPA);
      fWiresPerRop[cs][apa].resize(fPlanesPerAPA);
      fFirstChannelInThisRop[cs][apa].resize(fPlanesPerAPA);
      fFirstChannelInNextRop[cs][apa].resize(fPlanesPerAPA);

    }// end loop over apas
  }// end cryostats

  // Find the number of wires anchored to the frame
  for ( unsigned int icry=0; icry!=fNcryostat; ++icry ) {
    geo::CryostatGeo* pgeocry = cgeo[icry];
    for ( unsigned int iapa=0; iapa!=fNAPA[icry]; ++iapa ) {
      for ( unsigned int irop=0; irop!=fPlanesPerAPA; ++irop ) {
        unsigned int ipla = tpcPlaneForApaPlane[irop];
        unsigned int t = 2*iapa + tpcOffsetForApaPlane[irop];
        fWiresPerRop[icry][iapa][irop] = pgeocry->TPC(t).Plane(ipla).Nwires();
        fViews.emplace(pgeocry->TPC(t).Plane(ipla).View());
        if ( pgeocry->TPC(t).Plane(ipla).View() == geo::kZ ) { 
          nAnchoredWires[icry][iapa][irop] = fWiresPerRop[icry][iapa][irop];      
        } else {
          double xyz[3] = {0.};
          double xyz_next[3] = {0.};
          for ( unsigned int iwir=0; iwir!=fWiresPerRop[icry][iapa][irop]; ++iwir ) {
            pgeocry->TPC(t).Plane(ipla).Wire(iwir).GetCenter(xyz);
            pgeocry->TPC(t).Plane(ipla).Wire(iwir+1).GetCenter(xyz_next);
            if ( xyz[2] == xyz_next[2] ) {
              nAnchoredWires[icry][iapa][irop] = iwir; // w-1(for last)+1(for index) = w    
              break;
            }
          }  // End loop over wires
        }
      }  // end loop over planes
    }
  }
  
  raw::ChannelID_t currentChannel = 0;
 
  for ( unsigned int icry = 0; icry!= fNcryostat; ++icry ) {
    for ( unsigned int iapa = 0; iapa != fNAPA[icry]; ++iapa ) {  
      for ( unsigned int irop = 0; irop != fPlanesPerAPA; ++irop ) {
        fFirstChannelInThisRop[icry][iapa][irop] = currentChannel;
        unsigned int nchan = nAnchoredWires[icry][iapa][irop];
        if ( irop < 2 ) nchan += nchan;
        currentChannel += nchan;
        fFirstChannelInNextRop[icry][iapa][irop] = currentChannel;
        mf::LogVerbatim("Dune35tChannelMapAlg") << "Plane: Nwire NAnchored FIrstThis FirstNext: "
                      << icry << "-" << iapa << "-" << irop << ": "
                      << setw(4) << fWiresPerRop[icry][iapa][irop] << " "
                      << setw(4) << nAnchoredWires[icry][iapa][irop] << " "
                      << setw(4) << fFirstChannelInThisRop[icry][iapa][irop] << " "
                      << setw(4) << fFirstChannelInNextRop[icry][iapa][irop];
      }
    }
  }


  // Save the number of channels
  fNchannels = currentChannel;

  // Save the number of channels
  fChannelsPerAPA = fFirstChannelInNextRop[0][0][fPlanesPerAPA-1];

  fWirePitch.resize(fPlanesPerTPC);
  fOrientation.resize(fPlanesPerTPC);
  fSinOrientation.resize(fPlanesPerTPC);
  fCosOrientation.resize(fPlanesPerTPC);

  //save data into fFirstWireCenterY, fFirstWireCenterZ and fWireSortingInZ
  fPlaneData.resize(fNcryostat);
  for ( unsigned int icry=0; icry<fNcryostat; ++icry ){
    fPlaneData[icry].resize(fNTPC[icry]);
    for ( unsigned int tpc=0; tpc<fNTPC[icry]; tpc++){
      fPlaneData[icry][tpc].resize(fPlanesPerTPC);
      for ( unsigned int plane=0; plane<fPlanesPerTPC; plane++){
        PlaneData_t& PlaneData = fPlaneData[icry][tpc][plane];
        fPlaneIDs.emplace(icry, tpc, plane);
        double xyz[3]={0.0, 0.0, 0.0};
        const geo::PlaneGeo& thePlane = cgeo[icry]->TPC(tpc).Plane(plane);
        thePlane.Wire(0).GetCenter(xyz);
        PlaneData.fFirstWireCenterY = xyz[1];
        PlaneData.fFirstWireCenterZ = xyz[2];
        // we are interested in the ordering of wire numbers: we find that a
        // point is N wires left of a wire W: is that wire W + N or W - N?
        // In fact, for TPC #0 it is W + N for V and Z planes, W - N for U
        // plane; for TPC #0 it is W + N for V and Z planes, W - N for U
        PlaneData.fWireSortingInZ = thePlane.WireIDincreasesWithZ()? +1.: -1.;
        // find boundaries of the APA frame for this plane by looking at endpoints of wires
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
        }  // end loop over plane wires
      }  // end loop over TPC planes
    }  // end loop over TPCs
  }  // end loop over cryostats

  //initialize fWirePitch and fOrientation
  for ( unsigned int plane=0; plane<fPlanesPerTPC; ++plane ) {
    fWirePitch[plane]=cgeo[0]->TPC(0).WirePitch(0,1,plane);
    fOrientation[plane]=cgeo[0]->TPC(0).Plane(plane).Wire(0).ThetaZ();
    fSinOrientation[plane] = sin(fOrientation[plane]);
    fCosOrientation[plane] = cos(fOrientation[plane]);

  } // for
  mf::LogVerbatim("Dune35tChannelMapAlg") << "fNchannels = " << fNchannels ; 
  mf::LogVerbatim("Dune35tChannelMapAlg") << "U channels per APA = " << 2*nAnchoredWires[0][0][0] ;
  mf::LogVerbatim("Dune35tChannelMapAlg") << "V channels per APA = " << 2*nAnchoredWires[0][0][1] ;
  mf::LogVerbatim("Dune35tChannelMapAlg") << "Z channels per APA side = " << nAnchoredWires[0][0][2] ;
}
 
//----------------------------------------------------------------------------

void Dune35tChannelMapAlg::Uninitialize() {
  PlaneInfoMap_t<raw::ChannelID_t>().swap(fFirstChannelInThisRop);
  PlaneInfoMap_t<raw::ChannelID_t>().swap(fFirstChannelInNextRop);
}

//----------------------------------------------------------------------------

std::vector<geo::WireID> Dune35tChannelMapAlg::ChannelToWire(raw::ChannelID_t icha) const {
  // first check if this channel ID is legal
  std::vector<WireID> wirids;
  if ( icha >= fNchannels ) return wirids;
  unsigned int icry = 0;
  unsigned int iapa = 0;
  unsigned int irop = 0;
  unsigned int iwir = 0;
  bool foundPlane = false;
  for ( icry=0; icry!=fNcryostat; ++icry ) {
    for ( iapa=0; iapa!=fNAPA[icry]; ++iapa ) {
      for ( irop = 0; irop!=fPlanesPerAPA; ++irop ) {
        unsigned int icha1 = fFirstChannelInThisRop[icry][iapa][irop];
        unsigned int icha2 = fFirstChannelInNextRop[icry][iapa][irop];
        if ( icha >= icha1 && icha < icha2 ) {
          iwir = icha - icha1;
          foundPlane = true;
          break;
        }// end if break          
        if ( foundPlane ) break;
      }// end plane loop        
      if ( foundPlane ) break;
    }// end apa loop      
    if ( foundPlane ) break;
  }// end cryostat loop
  if ( ! foundPlane ) return wirids;
  unsigned int nAnchored = nAnchoredWires[icry][iapa][irop];
  bool wrapped = irop < 2;
  // Find the TPCs
  unsigned int itpc1 = 2*iapa;
  unsigned int itpc2 = itpc1 + 1;
  unsigned int itpc = itpc1;
  // If this is a wrapped ROP and the wire number is larger than nAnchored, then
  // the first wire for this channel is in the next TPC plane.
  if ( wrapped && iwir >= nAnchored ) {
    itpc = itpc2;
    iwir -= nAnchored;
  }
  if ( irop == 3 ) itpc = itpc2;
  if ( iwir >= nAnchored ) {
    cout << "iwir=" << iwir << ", nAnchored=" << nAnchored << endl;
    abort();
  }
  // Find the TPC plane
  unsigned int ipla = (irop == 3) ? 2 : irop;
  // Loop over wires and create IDs.
  unsigned int nwir = fWiresPerRop[icry][iapa][irop];
  while ( iwir < nwir ) {
    geo::WireID wirid(icry, itpc, ipla, iwir);
    wirids.push_back(wirid);
    iwir += nAnchored;
    itpc = (itpc == itpc1) ? itpc2 : itpc1;
  }
  return wirids;
}

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::Nchannels() const {
  return fNchannels;
}

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::Nchannels(ROPID const& ropid) const {
  unsigned int icry = ropid.Cryostat;
  unsigned int iapa = ropid.TPCset;
  unsigned int irop = ropid.ROP;
  return fFirstChannelInNextRop[icry][iapa][irop] - fFirstChannelInThisRop[icry][iapa][irop];
}
  
//----------------------------------------------------------------------------

double Dune35tChannelMapAlg::
WireCoordinate(double YPos, double ZPos, PlaneID const& planeid) const {
  // Returns the wire number corresponding to a (Y,Z) position in the plane
  // with float precision.
  // Core code ripped from original NearestWireID() implementation
  const PlaneData_t& PlaneData = AccessElement(fPlaneData, planeid);
  // the formula used here is geometric:
  // distance = delta_y cos(theta_z) + delta_z sin(theta_z)
  // with a correction for the orientation of the TPC:
  // odd TPCs have supplementary wire angle (pi-theta_z), changing cosine sign
  const bool bSuppl = (planeid.TPC % 2) == 1;
  float distance =
    -(YPos - PlaneData.fFirstWireCenterY) * (bSuppl? -1.: +1.) * fCosOrientation[planeid.Plane]
    +(ZPos - PlaneData.fFirstWireCenterZ) * fSinOrientation[planeid.Plane];
  // The sign of this formula is correct if the wire with larger ID is on the
  // "right" (intuitively, larger z; rigorously, smaller intercept) than this one.
  // Of course, we are not always that lucky. fWireSortingInZ fixes our luck.
  return PlaneData.fWireSortingInZ * distance/fWirePitch[planeid.Plane];
}
  
  
//----------------------------------------------------------------------------

WireID Dune35tChannelMapAlg::
NearestWireID(const TVector3& xyz, PlaneID const& planeid) const {
  // cap the position to be within the boundaries of the wire endpoints.
  // This simulates charge drifting in from outside of the wire frames inwards towards
  // the first and last collection wires on the side, and towards the right endpoints
  const PlaneData_t& PlaneData = AccessElement(fPlaneData, planeid);
  double ycap = std::max( PlaneData.fYmin, std::min(PlaneData.fYmax, xyz.Y()) );
  double zcap = std::max( PlaneData.fZmin, std::min(PlaneData.fZmax, xyz.Z()) );
  double xwir = WireCoordinate(ycap, zcap, planeid);
  if ( xwir < 0.0 ) xwir = 0.0;
  unsigned int iwir = 0.5 + xwir;
  unsigned int nwir = WiresPerPlane(planeid);
  if ( iwir >= nwir ) iwir = nwir - 1;
  return WireID(planeid, iwir);
}
  
//----------------------------------------------------------------------------

raw::ChannelID_t Dune35tChannelMapAlg::PlaneWireToChannel(geo::WireID const& wireid) const {
  unsigned int icry = wireid.Cryostat;
  unsigned int itpc = wireid.TPC;
  unsigned int ipla = wireid.Plane;
  unsigned int iwir = wireid.Wire;
  unsigned int iapa = itpc/2;
  unsigned int irop = ipla;
  bool induction = irop < 2;
  bool collection = !induction;
  bool otherSide = 2*iapa != itpc;
  if ( collection && otherSide ) irop = 3;
  unsigned int nwirAnchored = nAnchoredWires[icry][iapa][irop];

  // Evaluate the wire number in the ROP. For induction wires in the the second TPC in
  // the APA, we must offset by the # anchored wires in the first TPC.
  unsigned int iwirRop = iwir;
  if ( induction && otherSide ) iwirRop += nwirAnchored;

  // Evaluate the channel number in the ROP accounting for the wires wrapping with
  // a period of 2 time # anchored wires.
  unsigned int ichaRop = iwirRop % (2*nwirAnchored);

  // Add the ROP channel offset to get the detector channel number.
  unsigned int icha = ichaRop + fFirstChannelInThisRop[icry][iapa][irop];

  return icha;
}


//----------------------------------------------------------------------------

SigType_t Dune35tChannelMapAlg::SignalType(raw::ChannelID_t const icha) const {
  ROPID ropid = ChannelToROP(icha);
  unsigned int irop = ropid.ROP;
  if ( irop < 2 ) return geo::kInduction;
  if ( irop < 4 ) return geo::kCollection;
  return geo::kMysteryType;
}

//----------------------------------------------------------------------------

View_t Dune35tChannelMapAlg::View(raw::ChannelID_t const icha) const {
  ROPID ropid = ChannelToROP(icha);
  unsigned int irop = ropid.ROP;
  if ( irop < 1 ) return geo::kU;
  if ( irop < 2 ) return geo::kV;
  if ( irop < 4 ) return geo::kZ;
  return geo::kUnknown;
}  
 
//----------------------------------------------------------------------------

std::set<View_t> const& Dune35tChannelMapAlg::Views() const {
  return fViews;
}

//----------------------------------------------------------------------------

std::set<PlaneID> const& Dune35tChannelMapAlg::PlaneIDs() const {
  return fPlaneIDs;
}

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::NOpChannels(unsigned int NOpDets) const {
  return 12*NOpDets;
}

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::NOpHardwareChannels(unsigned int opDet) const {
  // CSU 3-sipm design
  if (opDet == 0 || opDet == 4 || opDet == 6) return 8;
  // LSU 2-sipm design
  if (opDet == 2) return 2;
  // IU 12-sipm design
  return 12;    
}

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::OpChannel(unsigned int detNum, unsigned int channel) const {
  unsigned int uniqueChannel = (detNum * 12) + channel;
  return uniqueChannel;
}

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::OpDetFromOpChannel(unsigned int opChannel) const {
  unsigned int detectorNum = (unsigned int) opChannel / 12;
  return detectorNum;
}

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::HardwareChannelFromOpChannel(unsigned int opChannel) const {
  unsigned int channel = opChannel % 12;
  return channel;
}

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::
NTPCsets(CryostatID const& cryoid) const { return 4; }

//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::MaxTPCsets() const { return 4; }
  
//----------------------------------------------------------------------------

bool Dune35tChannelMapAlg::HasTPCset(TPCsetID const& tpcsetid) const {
  return tpcsetid.TPCset < NTPCsets(tpcsetid);
}

//----------------------------------------------------------------------------

TPCsetID Dune35tChannelMapAlg::TPCtoTPCset(TPCID const& tpcid) const {
  unsigned int itpc = tpcid.TPC;
  unsigned int iapa = itpc/2;
  CryostatID cryid = tpcid;
  return TPCsetID(cryid, iapa);
  
}
  
//----------------------------------------------------------------------------

std::vector<TPCID>
Dune35tChannelMapAlg::TPCsetToTPCs(TPCsetID const& apaid) const {
  std::vector<TPCID> tpcids;
  CryostatID cryid = apaid;
  unsigned int iapa = apaid.TPCset;
  tpcids.push_back(TPCID(cryid, 2*iapa));
  tpcids.push_back(TPCID(cryid, 2*iapa+1));
  return tpcids;
}
  
//----------------------------------------------------------------------------

TPCID Dune35tChannelMapAlg::
FirstTPCinTPCset(TPCsetID const& apaid) const {
  std::vector<TPCID> tpcids;
  CryostatID cryid = apaid;
  unsigned int iapa = apaid.TPCset;
  return TPCID(cryid, 2*iapa);
}
  
//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::NROPs(TPCsetID const&) const {
  return 4;
}
  
//----------------------------------------------------------------------------

unsigned int Dune35tChannelMapAlg::MaxROPs() const {
  return 4;
}
  
//----------------------------------------------------------------------------

bool Dune35tChannelMapAlg::HasROP(ROPID const& ropid) const {
  return ropid.ROP < NROPs(ropid);
}
  
  
//----------------------------------------------------------------------------

ROPID Dune35tChannelMapAlg::WirePlaneToROP (PlaneID const& planeid) const {
  unsigned int icry = planeid.Cryostat;
  unsigned int ipla = planeid.Plane;
  unsigned int itpc = planeid.TPC;
  unsigned int iapa = itpc/2;
  unsigned int irop = ROPID::InvalidID;
  if ( itpc < 2 ) {
    irop = ipla;
  } else if ( itpc < 4 ) {
    bool tpcIsEven = 2*(itpc/2) == itpc;
    irop = tpcIsEven ? 2 : 3;
  } else {
    iapa = TPCsetID::InvalidID;
  }
  return ROPID(icry, iapa, irop);
}

//----------------------------------------------------------------------------

std::vector<PlaneID> Dune35tChannelMapAlg::ROPtoWirePlanes(ROPID const& ropid) const {
  std::vector<PlaneID> plaids;
  std::vector<TPCID> tpcids = ROPtoTPCs(ropid);
  unsigned int irop = ropid.ROP;
  for ( TPCID tpcid : tpcids ) {
    if ( irop == 0 || irop == 1 ) plaids.push_back(PlaneID(tpcid, irop));
    else if ( irop == 2 || irop == 3 ) plaids.push_back(PlaneID(tpcid, 2));
  }
  return plaids;
}
  
//----------------------------------------------------------------------------
  
std::vector<TPCID> Dune35tChannelMapAlg::ROPtoTPCs(ROPID const& ropid) const {
  std::vector<unsigned int> itpcs;
  unsigned int irop = ropid.ROP;
  unsigned int iapa = ropid.TPCset;
  if ( irop == 0 || irop == 1 || irop == 2 ) itpcs.push_back(2*iapa);
  if ( irop == 0 || irop == 1 || irop == 3 ) itpcs.push_back(2*iapa + 1);
  std::vector<TPCID> tpcids;
  CryostatID cryid(ropid);
  for ( unsigned int itpc : itpcs ) tpcids.push_back(TPCID(cryid, itpc));
  return tpcids;
}
  
//----------------------------------------------------------------------------

ROPID Dune35tChannelMapAlg::ChannelToROP(raw::ChannelID_t icha) const {
  for ( unsigned int icry=0; fNcryostat; ++icry ) {
    for ( unsigned int iapa=0; iapa<fNAPA[icry]; ++iapa ) {
      for ( unsigned int irop=0; irop<fPlanesPerAPA; ++irop ) {
        if ( icha >= fFirstChannelInThisRop[icry][iapa][irop] &&
             icha < fFirstChannelInNextRop[icry][iapa][irop] ) {
          return ROPID(icry, iapa, irop);
        }
      }
    }
  }
  return  ROPID(CryostatID::InvalidID, TPCsetID::InvalidID, ROPID::InvalidID);
}
  
//----------------------------------------------------------------------------

raw::ChannelID_t Dune35tChannelMapAlg::FirstChannelInROP(ROPID const& ropid) const {
  unsigned int icry = ropid.Cryostat;
  unsigned int iapa = ropid.TPCset;
  unsigned int irop = ropid.ROP;
  return fFirstChannelInThisRop[icry][iapa][irop];
}
  
//----------------------------------------------------------------------------

PlaneID Dune35tChannelMapAlg::FirstWirePlaneInROP(ROPID const& ropid) const {
  PlaneID pid = ROPtoWirePlanes(ropid).at(0);
  return pid;
}
  
//----------------------------------------------------------------------------

}  // end geo namespace
