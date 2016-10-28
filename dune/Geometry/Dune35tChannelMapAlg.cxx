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
  fWiresPerPlane.resize(fNcryostat);
  fFirstChannelInNextPlane.resize(fNcryostat);
  fFirstChannelInThisPlane.resize(fNcryostat);
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
    fWiresPerPlane[cs].resize(fNTPC[cs]);
    fFirstChannelInThisPlane[cs].resize(fNAPA[cs]);
    fFirstChannelInNextPlane[cs].resize(fNAPA[cs]);

    for(unsigned int apa = 0; apa != fNAPA[cs]; ++apa){
      
      nAnchoredWires[cs][apa].resize(fPlanesPerAPA);
      fWiresPerPlane[cs][apa].resize(fPlanesPerAPA);
      fFirstChannelInThisPlane[cs][apa].resize(fPlanesPerAPA);
      fFirstChannelInNextPlane[cs][apa].resize(fPlanesPerAPA);

    }// end loop over apas
  }// end cryostats

  // Find the number of wires anchored to the frame
  for(unsigned int c = 0; c != fNcryostat; ++c){
    for(unsigned int a = 0; a != fNAPA[c]; ++a){
      for(unsigned int ap = 0; ap != fPlanesPerAPA; ++ap){
        unsigned int tp = tpcPlaneForApaPlane[ap];
        unsigned int t = a + tpcOffsetForApaPlane[ap];
        fWiresPerPlane[c][a][ap] = cgeo[c]->TPC(t).Plane(tp).Nwires();
        double xyz[3] = {0.};
        double xyz_next[3] = {0.};

        fViews.emplace(cgeo[c]->TPC(t).Plane(tp).View());

        for(unsigned int w = 0; w != fWiresPerPlane[c][a][ap]; ++w){

          // for vertical planes
          if(cgeo[c]->TPC(t).Plane(tp).View() == geo::kZ)   { 
            nAnchoredWires[c][a][ap] = fWiresPerPlane[c][a][ap];      
            break;
          }

          cgeo[c]->TPC(t).Plane(tp).Wire(w).GetCenter(xyz);
          cgeo[c]->TPC(t).Plane(tp).Wire(w+1).GetCenter(xyz_next);

          if(xyz[2]==xyz_next[2]){
            nAnchoredWires[c][a][ap] = w; // w-1(for last)+1(for index) = w    
            break;
          }

        }
      }  // end loop over planes
      // If the wires are missing for the second collection plane, assume
      // it is the same as the first.
      //unsigned int p = fPlanesPerTPC;
      //if ( p < fPlanesPerAPA && fWiresPerPlane[c][a][p] == 0 ) {
      //  fWiresPerPlane[c][a][p] = fWiresPerPlane[c][a][p-1];
      //}
    }
  }
  

  raw::ChannelID_t currentChannel = 0;
 
  for(unsigned int cs = 0; cs != fNcryostat; ++cs){
    for(unsigned int apa = 0; apa != fNAPA[cs]; ++apa){  
      for(unsigned int p = 0; p != fPlanesPerAPA; ++p){

        fFirstChannelInThisPlane[cs][apa][p] = currentChannel;
        unsigned int nchan = nAnchoredWires[cs][apa][p];
        if ( p < 2 ) nchan += nchan;
        currentChannel += nchan;
        fFirstChannelInNextPlane[cs][apa][p] = currentChannel;
        mf::LogVerbatim("Dune35tChannelMapAlg") << "Plane: Nwire NAnchored FIrstThis FirstNext: "
                      << cs << "-" << apa << "-" << p << ": "
                      << setw(4) << fWiresPerPlane[cs][apa][p] << " "
                      << setw(4) << nAnchoredWires[cs][apa][p] << " "
                      << setw(4) << fFirstChannelInThisPlane[cs][apa][p] << " "
                      << setw(4) << fFirstChannelInNextPlane[cs][apa][p];
      }
    }
  }


  // Save the number of channels
  fNchannels = currentChannel;

  // Save the number of channels
  fChannelsPerAPA = fFirstChannelInNextPlane[0][0][fPlanesPerAPA-1];

  
  fWirePitch.resize(fPlanesPerTPC);
  fOrientation.resize(fPlanesPerTPC);
  fSinOrientation.resize(fPlanesPerTPC);
  fCosOrientation.resize(fPlanesPerTPC);

  //save data into fFirstWireCenterY, fFirstWireCenterZ and fWireSortingInZ
  fPlaneData.resize(fNcryostat);
  for (unsigned int cs=0; cs<fNcryostat; cs++){
    fPlaneData[cs].resize(fNTPC[cs]);
    for (unsigned int tpc=0; tpc<fNTPC[cs]; tpc++){
      fPlaneData[cs][tpc].resize(fPlanesPerTPC);
      for (unsigned int plane=0; plane<fPlanesPerTPC; plane++){
        PlaneData_t& PlaneData = fPlaneData[cs][tpc][plane];
        fPlaneIDs.emplace(cs, tpc, plane);
        double xyz[3]={0.0, 0.0, 0.0};
        const geo::PlaneGeo& thePlane = cgeo[cs]->TPC(tpc).Plane(plane);
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

  return;

}
 
//----------------------------------------------------------------------------

void Dune35tChannelMapAlg::Uninitialize() {
  PlaneInfoMap_t<raw::ChannelID_t>().swap(fFirstChannelInThisPlane);
  PlaneInfoMap_t<raw::ChannelID_t>().swap(fFirstChannelInNextPlane);
}

//----------------------------------------------------------------------------

std::vector<geo::WireID> Dune35tChannelMapAlg::ChannelToWire(raw::ChannelID_t channel) const {

  // first check if this channel ID is legal
  if(channel >= fNchannels )
    throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel << "\n";

  std::vector< WireID > AllSegments;
    
  static unsigned int cstat;
  static unsigned int tpc;
  static unsigned int plane;
  static unsigned int wireThisPlane;
  static unsigned int NextPlane;
  static unsigned int ThisPlane;
    
  for(unsigned int csloop = 0; csloop != fNcryostat; ++csloop){
      
    bool breakVariable = false;
      
    for(unsigned int apaloop = 0; apaloop != fNAPA[csloop]; ++apaloop){
      for(unsigned int planeloop = 0; planeloop != fPlanesPerTPC; ++planeloop){
        NextPlane = fFirstChannelInNextPlane[csloop][apaloop][planeloop];
        ThisPlane = fFirstChannelInThisPlane[csloop][apaloop][planeloop];
        if(channel < NextPlane){
          cstat = csloop;
          tpc   = 2*apaloop;
          plane = planeloop;
          wireThisPlane  = channel - ThisPlane;
          breakVariable = true;
          break;
        }// end if break          
        if(breakVariable) break;
      }// end plane loop        
      if(breakVariable) break;
    }// end apa loop      
    if(breakVariable) break;
      
  }// end cryostat loop
    

  int WrapDirection = 1; // go from tpc to (tpc+1) or tpc to (tpc-1)

  // find the lowest wire
  raw::ChannelID_t ChannelGroup = std::floor( wireThisPlane/nAnchoredWires[cstat][tpc/2][plane] );
  unsigned int bottomwire = wireThisPlane-ChannelGroup*nAnchoredWires[cstat][tpc/2][plane];
    
  if(ChannelGroup%2==1){
    // start in the other TPC
    tpc += 1;
    WrapDirection  = -1;         
  }
    
  for(unsigned int WireSegmentCount = 0; WireSegmentCount != 50; ++WireSegmentCount){
      
    tpc += WrapDirection*(WireSegmentCount%2);
      
    geo::WireID CodeWire(cstat, tpc, plane, bottomwire + WireSegmentCount*nAnchoredWires[cstat][std::floor(tpc/2)][plane]);
      
    AllSegments.push_back(CodeWire);
      
    // reset the tcp variable so it doesnt "accumulate value"
    tpc -= WrapDirection*(WireSegmentCount%2);
      
    if( bottomwire + (WireSegmentCount+1)*nAnchoredWires[cstat][std::floor(tpc/2)][plane] > 
      fWiresPerPlane[cstat][std::floor(tpc/2)][plane]-1) break;
      
  } //end WireSegmentCount loop
    
    return AllSegments;
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
  return fFirstChannelInNextPlane[icry][iapa][irop] - fFirstChannelInThisPlane[icry][iapa][irop];
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
    +(ZPos - PlaneData.fFirstWireCenterZ) * fSinOrientation[planeid.Plane]
    ;
    
  // The sign of this formula is correct if the wire with larger ID is on the
  // "right" (intuitively, larger z; rigorously, smaller intercept)
  // than this one.
  // Of course, we are not always that lucky. fWireSortingInZ fixes our luck.

  return PlaneData.fWireSortingInZ * distance/fWirePitch[planeid.Plane];
} // Dune35tChannelMapAlg::WireCoordinate()
  
  
//----------------------------------------------------------------------------

WireID Dune35tChannelMapAlg::
NearestWireID(const TVector3& xyz, PlaneID const& planeid) const {

  // cap the position to be within the boundaries of the wire endpoints.
  // This simulates charge drifting in from outside of the wire frames inwards towards
  // the first and last collection wires on the side, and towards the right endpoints

  const PlaneData_t& PlaneData = AccessElement(fPlaneData, planeid);
  double ycap = std::max(PlaneData.fYmin,std::min(PlaneData.fYmax,xyz.Y()));
  double zcap = std::max(PlaneData.fZmin,std::min(PlaneData.fZmax,xyz.Z()));

  // add 0.5 to have the correct rounding
  int NearestWireNumber
    = int (0.5 + WireCoordinate(ycap, zcap, planeid));
    
  // If we are outside of the wireplane range, throw an exception
  // (this response maintains consistency with the previous
  // implementation based on geometry lookup)
  if(NearestWireNumber < 0 ||
     NearestWireNumber >= (int) WiresPerPlane(planeid)) {
    const int wireNumber = NearestWireNumber; // save for the output
      
    if(wireNumber < 0 ) NearestWireNumber = 0;
    else                NearestWireNumber = WiresPerPlane(planeid) - 1;
    
  /*
    // comment in the following statement to throw an exception instead
    throw InvalidWireIDError("Geometry", wireNumber, NearestWireNumber)
      << "Dune35tChannelMapAlg::NearestWireID(): can't Find Nearest Wire for position (" 
      << xyz.X() << "," << xyz.Y() << "," << xyz.Z() << ")"
      << " approx wire number # " << wireNumber
      << " (capped from " << NearestWireNumber << ")\n";
  */
  } // if invalid wire
    
  return { planeid, (geo::WireID::WireID_t) NearestWireNumber };
} // Dune35tChannelMapAlg::NearestWireID()
  
//----------------------------------------------------------------------------

raw::ChannelID_t Dune35tChannelMapAlg::PlaneWireToChannel (geo::WireID const& wireid) const {

  unsigned int OtherSideWires = 0;

  raw::ChannelID_t Channel = AccessAPAelement(fFirstChannelInThisPlane, wireid);

  // get number of wires starting on the first side of the APA if starting
  // on the other side TPC.
  if (wireid.TPC % 2 == 1) OtherSideWires += AnchoredWires(wireid);
  
  // Lastly, account for the fact that channel number while moving up wire number in one
  // plane resets after 2 times the number of wires anchored -- one for each APA side.
  // At the same time, OtherSideWires accounts for the fact that if a channel starts on
  // the other side, it is offset by the number of wires on the first side.
  Channel += (OtherSideWires + wireid.Wire)%(2*AnchoredWires(wireid));

  return Channel;
}


//----------------------------------------------------------------------------

SigType_t Dune35tChannelMapAlg::SignalType( raw::ChannelID_t const channel ) const {
  raw::ChannelID_t chan = channel % fChannelsPerAPA;
  SigType_t sigt = kInduction;
  if (       chan <  fFirstChannelInThisPlane[0][0][2]     ) {
    sigt = kInduction;
  } else {
    if ( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
         (chan <  fFirstChannelInNextPlane[0][0][2])    ) {
      sigt = kCollection;
    } else {
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel << " (" << chan
                                             << ") not given signal type." << std::endl;
    }
  }
  return sigt;
}

//----------------------------------------------------------------------------

View_t Dune35tChannelMapAlg::View( raw::ChannelID_t const channel ) const {
  raw::ChannelID_t chan = channel % fChannelsPerAPA;
  View_t view = geo::kU;
  if(       chan <  fFirstChannelInNextPlane[0][0][0]     ){ view = geo::kU; }
  else if( (chan >= fFirstChannelInThisPlane[0][0][1]) &&
           (chan <  fFirstChannelInNextPlane[0][0][1])    ){ view = geo::kV; }
  else if( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
           (chan <  fFirstChannelInNextPlane[0][0][2])    ){ view = geo::kZ; }
  else{    mf::LogWarning("BadChannelViewType") << "Channel " << channel 
                                                << " (" << chan << ") not given view type.";}
  return view;
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
  if (opDet == 0 || opDet == 4 || opDet == 6)
    return 8;

  // LSU 2-sipm design
  if (opDet == 2)
    return 2;

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

ROPID Dune35tChannelMapAlg::ChannelToROP (raw::ChannelID_t icha) const {
  for ( unsigned int icry=0; fNcryostat; ++icry ) {
    for ( unsigned int iapa=0; iapa<fNAPA[icry]; ++iapa ) {
      for ( unsigned int irop=0; irop<fPlanesPerAPA; ++irop ) {
        if ( icha >= fFirstChannelInThisPlane[icry][iapa][irop] &&
             icha < fFirstChannelInNextPlane[icry][iapa][irop] ) {
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
  return fFirstChannelInThisPlane[icry][iapa][irop];
}
  
//----------------------------------------------------------------------------

PlaneID Dune35tChannelMapAlg::FirstWirePlaneInROP(ROPID const& ropid) const {
  PlaneID pid = ROPtoWirePlanes(ropid).at(0);
  return pid;
}
  
//----------------------------------------------------------------------------

}  // end geo namespace
