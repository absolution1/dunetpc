////////////////////////////////////////////////////////////////////////
/// \file  DuneApaChannelMapAlg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////

#include "dune/Geometry/DuneApaChannelMapAlg.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

typedef unsigned int Index;

#include <iostream>
using std::cout;
using std::endl;

namespace geo{

//----------------------------------------------------------------------------

DuneApaChannelMapAlg::DuneApaChannelMapAlg(fhicl::ParameterSet const& p)
: fSorter(geo::GeoObjectSorterAPA(p)) {
  fChannelsPerOpDet = p.get<unsigned int>("ChannelsPerOpDet");
}

//----------------------------------------------------------------------------

void DuneApaChannelMapAlg::Initialize( GeometryData_t& geodata ) {
  Uninitialize();
    
  std::vector<geo::CryostatGeo*>& crygeos = geodata.cryostats;
  std::vector<geo::AuxDetGeo*>  & adgeo = geodata.auxDets;
  fNcryostat = crygeos.size();
  if ( fNcryostat == 0 ) {
    mf::LogError("DuneApaChannelMapAlg") << "No cryostats found.";
    return;
  }
  geo::CryostatGeo& crygeo = *crygeos[0];
    
  mf::LogInfo("DuneApaChannelMapAlg") << "Sorting volumes...";

  fSorter.SortAuxDets(adgeo);
  fSorter.SortCryostats(crygeos);
  for ( Index icry=0; icry<crygeos.size(); ++icry ) {
    crygeos[icry]->SortSubVolumes(fSorter);
  }

  mf::LogInfo("DuneApaChannelMapAlg") << "Initializing channel map...";
  fNTPC.resize(fNcryostat);
  fFirstChannelInNextPlane.resize(1);  // Change 1 to Ncryostat if you want to treat each APA
  fFirstChannelInThisPlane.resize(1);  // uniquely,and do // the same with the other resizes.
  fPlanesPerTPC = crygeo.TPC(0).Nplanes();
  if ( fPlanesPerTPC != 3 ) {
    mf::LogError("DuneApaChannelMapAlg") << "# planes/TPC is not 3: " << fPlanesPerTPC;
    return;
  }
  fPlanesPerAPA = fPlanesPerTPC + 1;
  nAnchoredWires.resize(fPlanesPerTPC);
  fWiresInPlane.resize(fPlanesPerTPC);
  fFirstChannelInThisPlane[0].resize(1);  // remember FirstChannel vectors
  fFirstChannelInNextPlane[0].resize(1);  // for first APA only.
  fFirstChannelInThisPlane[0][0].resize(fPlanesPerTPC);  // Make room for info
  fFirstChannelInNextPlane[0][0].resize(fPlanesPerTPC);  // on each plane.
  fViews.clear();
  fPlaneIDs.clear();
  fTopChannel = 0;

  // Size some vectors and initialize the FirstChannel vectors.
  // If making FirstChannel's for every APA uniquely, they also
  // need to be sized here. Not necessary for now
  for ( unsigned int icry = 0; icry != fNcryostat; ++icry ) {
    fNTPC[icry] = crygeos[icry]->NTPC();
  }

  // Find the number of wires anchored to the frame
  for ( unsigned int ipla=0; ipla!=fPlanesPerTPC; ++ipla ) {
    fWiresInPlane[ipla] = crygeo.TPC(0).Plane(ipla).Nwires();
    double xyz[3] = {0.};
    double xyz_next[3] = {0.};
    const PlaneGeo& plageo = crygeo.TPC(0).Plane(ipla);
#if 0
const PlaneGeo plageo2 = plageo;
#endif
    fViews.emplace(plageo.View());
    for ( unsigned int iwir = 0; iwir+1<fWiresInPlane[ipla]; ++iwir ) {
      // for vertical planes
      if( plageo.View()==geo::kZ)   { 
        nAnchoredWires[ipla] = fWiresInPlane[ipla];      
        break;
      }
      plageo.Wire(iwir).GetCenter(xyz);
      plageo.Wire(iwir+1).GetCenter(xyz_next);
      if ( xyz[2] == xyz_next[2] ) {
        nAnchoredWires[ipla] = iwir;      
        break;
      }
    }
  }

    raw::ChannelID_t CurrentChannel = 0;
   
    for(unsigned int PCount = 0; PCount != fPlanesPerTPC; ++PCount){

      fFirstChannelInThisPlane[0][0][PCount] = CurrentChannel;
      CurrentChannel = CurrentChannel + 2*nAnchoredWires[PCount];
      fFirstChannelInNextPlane[0][0][PCount] = CurrentChannel;

    }// end build loop over planes

    // Save the number of channels
    fChannelsPerAPA = fFirstChannelInNextPlane[0][0][fPlanesPerTPC-1];

    fNchannels = 0;
    for(size_t cs = 0; cs < fNcryostat; ++cs){
      fNchannels = fNchannels + fChannelsPerAPA*fNTPC[cs]/2;
    }

    //save data into fFirstWireCenterY and fFirstWireCenterZ
    fPlaneData.resize(fNcryostat);
    for (unsigned int cs=0; cs<fNcryostat; ++cs){
      fPlaneData[cs].resize(crygeos[cs]->NTPC());
      for (unsigned int tpc=0; tpc<crygeos[cs]->NTPC(); ++tpc){
        fPlaneData[cs][tpc].resize(crygeos[cs]->TPC(tpc).Nplanes());
        for (unsigned int plane=0; plane<crygeos[cs]->TPC(tpc).Nplanes(); ++plane){
          PlaneData_t& PlaneData = fPlaneData[cs][tpc][plane];
          const geo::PlaneGeo& thePlane = crygeos[cs]->TPC(tpc).Plane(plane);
          double xyz[3];
          fPlaneIDs.emplace(cs, tpc, plane);
          thePlane.Wire(0).GetCenter(xyz);
          PlaneData.fFirstWireCenterY = xyz[1];
          PlaneData.fFirstWireCenterZ = xyz[2];
          // we are interested in the ordering of wire numbers: we find that a
          // point is N wires left of a wire W: is that wire W + N or W - N?
          // In fact, for TPC #0 it is W + N for V and Z planes, W - N for U
          // plane; for TPC #0 it is W + N for V and Z planes, W - N for U
          PlaneData.fWireSortingInZ = thePlane.WireIDincreasesWithZ()? +1.: -1.;
        } // for plane
      } // for TPC
    } // for cryostat

    //initialize fWirePitch and fOrientation
    Index npla = crygeo.TPC(0).Nplanes();
    fWirePitch.resize(npla);
    fOrientation.resize(npla);
    fSinOrientation.resize(npla);
    fCosOrientation.resize(npla);

    for ( unsigned int ipla=0; ipla<npla; ++ipla ) {
      fWirePitch[ipla] = crygeo.TPC(0).WirePitch(0,1,ipla);
      fOrientation[ipla] = crygeo.TPC(0).Plane(ipla).Wire(0).ThetaZ();
      fSinOrientation[ipla] = sin(fOrientation[ipla]);
      fCosOrientation[ipla] = cos(fOrientation[ipla]);
    }


    for(unsigned int c=0; c<fNcryostat; c++){

      mf::LogVerbatim("DuneApaChannelMapAlg") << "Cryostat " << c << ":"; 

      mf::LogVerbatim("DuneApaChannelMapAlg") << "  " << fNchannels << " total channels"; 
      mf::LogVerbatim("DuneApaChannelMapAlg") << "  " << fNTPC[c]/2 << " APAs";       
      mf::LogVerbatim("DuneApaChannelMapAlg") << "  For all identical APA:" ; 
      mf::LogVerbatim("DuneApaChannelMapAlg") << "    Number of channels per APA = " << fChannelsPerAPA ; 
      
      mf::LogVerbatim("DuneApaChannelMapAlg") << "    U channels per APA = " << 2*nAnchoredWires[0] ;
      mf::LogVerbatim("DuneApaChannelMapAlg") << "    V channels per APA = " << 2*nAnchoredWires[1] ;
      mf::LogVerbatim("DuneApaChannelMapAlg") << "    Z channels per APA = " << 2*nAnchoredWires[2] ;
      
      mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in U Plane = " << fWirePitch[0] ;
      mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in V Plane = " << fWirePitch[1] ;
      mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in Z Plane = " << fWirePitch[2] ;
      
    }

    return;

  }
   
  //----------------------------------------------------------------------------
  void DuneApaChannelMapAlg::Uninitialize()
  {

    PlaneInfoMap_t<raw::ChannelID_t>().swap(fFirstChannelInThisPlane);
    PlaneInfoMap_t<raw::ChannelID_t>().swap(fFirstChannelInNextPlane);

  }

  //----------------------------------------------------------------------------
  std::vector<geo::WireID> DuneApaChannelMapAlg::ChannelToWire(raw::ChannelID_t channel)  const
  {

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
    
    raw::ChannelID_t chan       = channel%fChannelsPerAPA;
    raw::ChannelID_t pureAPAnum = std::floor( channel/fChannelsPerAPA );

    bool breakVariable = false;
    for(unsigned int planeloop = 0; planeloop != fPlanesPerTPC; ++planeloop){
          
      NextPlane = fFirstChannelInNextPlane[0][0][planeloop];
      ThisPlane = fFirstChannelInThisPlane[0][0][planeloop];
          
      if(chan < NextPlane){

        // fNTPC[0] works for now since there are the same number of TPCs per crostat
        cstat = std::floor( channel / ((fNTPC[0]/2)*fChannelsPerAPA) );
        tpc   = (2*pureAPAnum) % fNTPC[0];
        plane = planeloop;
        wireThisPlane  = chan - ThisPlane;
            
        breakVariable = true;
        break;
      }
      if(breakVariable) break;          
    }// end plane loop

    

    int WrapDirection = 1; // go from tpc to (tpc+1) or tpc to (tpc-1)

    // find the lowest wire
    raw::ChannelID_t ChannelGroup = std::floor( wireThisPlane/nAnchoredWires[plane] );
    unsigned int bottomwire = wireThisPlane-ChannelGroup*nAnchoredWires[plane];
    
    if(ChannelGroup%2==1){
      tpc += 1;
      WrapDirection  = -1;         
    }
    
    for(unsigned int SegCount = 0; SegCount != 50; ++SegCount){
      
      tpc += WrapDirection*(SegCount%2);
      geo::WireID CodeWire(cstat, tpc, plane, bottomwire + SegCount*nAnchoredWires[plane]);
      AllSegments.push_back(CodeWire);
      
      // reset the tcp variable so it doesnt "accumulate value"
      tpc -= WrapDirection*(SegCount%2);
      
      if( bottomwire + (SegCount+1)*nAnchoredWires[plane] > fWiresInPlane[plane]-1) break;
    }
    
    return AllSegments;
  }


  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::Nchannels() const
  {
    return fNchannels;
  }
  
  
  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::Nchannels
    (readout::ROPID const& ropid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }

  //----------------------------------------------------------------------------
  
double DuneApaChannelMapAlg::
WireCoordinate(double YPos, double ZPos, unsigned int PlaneNo, unsigned int TPCNo,
               unsigned int cstat) const {
  return WireCoordinate(YPos, ZPos, geo::PlaneID(cstat, TPCNo, PlaneNo));
}
  
  //----------------------------------------------------------------------------
  double DuneApaChannelMapAlg::WireCoordinate
    (double YPos, double ZPos, geo::PlaneID const& planeid) const
  {
    // Returns the wire number corresponding to a (Y,Z) position in the plane
    // with float precision.
    // Core code ripped from original NearestWireID() implementation
    
    const PlaneData_t& PlaneData = AccessElement(fPlaneData, planeid);
    
    //get the orientation angle of a given plane and calculate the distance between first wire
    //and a point projected in the plane
//    const double rotate = (planeid.TPC % 2 == 1)? -1.: +1.;
    
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
  }
  
  //----------------------------------------------------------------------------
  WireID DuneApaChannelMapAlg::NearestWireID
    (const TVector3& xyz, geo::PlaneID const& planeid) const
  {
    // add 0.5 to have the correct rounding
    int NearestWireNumber
      = int (0.5 + WireCoordinate(xyz.Y(), xyz.Z(), planeid));
    
    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if(NearestWireNumber < 0 ||
       NearestWireNumber >= (int) fWiresInPlane[planeid.Plane])
    {
      const int wireNumber = NearestWireNumber; // save for the output
      
      if(wireNumber < 0 ) NearestWireNumber = 0;
      else                NearestWireNumber = fWiresInPlane[planeid.Plane] - 1;
    
    /*
      // comment in the following statement to throw an exception instead
      throw InvalidWireIDError("Geometry", wireNumber, NearestWireNumber)
        << "DuneApaChannelMapAlg::NearestWireID(): can't Find Nearest Wire for position (" 
        << xyz.X() << "," << xyz.Y() << "," << xyz.Z() << ")"
        << " approx wire number # " << wireNumber
        << " (capped from " << NearestWireNumber << ")\n";
    */
    } // if invalid wire
    
    return { planeid, (geo::WireID::WireID_t) NearestWireNumber };
  }
  
  
  //----------------------------------------------------------------------------
  raw::ChannelID_t DuneApaChannelMapAlg::PlaneWireToChannel
    (geo::WireID const& wireid) const
  {
    unsigned int OtherSideWires = 0;

    // start in very first APA:
    raw::ChannelID_t Channel = fFirstChannelInThisPlane[0][0][wireid.Plane];
    // move channel to proper cryostat:
    Channel += wireid.Cryostat*(fNTPC[wireid.Cryostat]/2)*fChannelsPerAPA;
    // move channel to proper APA:
    Channel += std::floor( wireid.TPC/2 )*fChannelsPerAPA;
    // get number of wires on the first side of the APA if starting on the other side TPC:
    OtherSideWires += (wireid.TPC%2)*nAnchoredWires[wireid.Plane];


    // Lastly, account for the fact that channel number while moving up wire number in one
    // plane resets after 2 times the number of wires anchored -- one for each APA side.
    // At the same time, OtherSideWires accounts for the fact that if a channel starts on 
    // the other side, it is offset by the number of wires on the first side.
    Channel += (OtherSideWires + wireid.Wire)%(2*nAnchoredWires[wireid.Plane]);
    
    return Channel;

  }

  //----------------------------------------------------------------------------
  SigType_t DuneApaChannelMapAlg::SignalType( raw::ChannelID_t const channel )  const
  {
    raw::ChannelID_t chan = channel % fChannelsPerAPA;
    SigType_t sigt = kInduction;

    // instead of calling channel to wire, we can make use of the way we 
    // up channels among APAs;
    // the first two planes are induction, and the last one is collection
    if(       chan <  fFirstChannelInThisPlane[0][0][2]     ){ sigt = kInduction;  }
    else if( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
             (chan <  fFirstChannelInNextPlane[0][0][2])    ){ sigt = kCollection; }
    else{    mf::LogWarning("BadChannelSignalType") << "Channel " << channel 
                                                    << " (" << chan << ") not given signal type." << std::endl;         }
  
    return sigt;
  }

  //----------------------------------------------------------------------------
  View_t DuneApaChannelMapAlg::View( raw::ChannelID_t const channel )  const
  {
    raw::ChannelID_t chan = channel % fChannelsPerAPA;
    View_t view = geo::kU;

    // instead of calling channel to wire, we can make use of the way we 
    // up channels among APAs;
    // Plane 0: U, Plane 1: V, Plane 2: W
    if(       chan <  fFirstChannelInNextPlane[0][0][0]     ){ view = geo::kU; }
    else if( (chan >= fFirstChannelInThisPlane[0][0][1]) &&
             (chan <  fFirstChannelInNextPlane[0][0][1])    ){ view = geo::kV; }
    else if( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
             (chan <  fFirstChannelInNextPlane[0][0][2])    ){ view = geo::kZ; }
    else{    mf::LogWarning("BadChannelViewType") << "Channel " << channel 
                                                  << " (" << chan << ") not given view type.";  }
    
    return view;
  }  

  //----------------------------------------------------------------------------
  std::set<View_t> const& DuneApaChannelMapAlg::Views() const
  {
    return fViews;
  }

  //----------------------------------------------------------------------------
  std::set<PlaneID> const& DuneApaChannelMapAlg::PlaneIDs() const
  {
    return fPlaneIDs;
  }

  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::NOpChannels(unsigned int NOpDets) const
  {
    return fChannelsPerOpDet*NOpDets;
  }

  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::NOpHardwareChannels(unsigned int opDet) const
  {
    return fChannelsPerOpDet;
  }

  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::OpChannel(unsigned int detNum, unsigned int channel) const
  {
    unsigned int uniqueChannel = (detNum * fChannelsPerOpDet) + channel;
    return uniqueChannel;
  }

  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::OpDetFromOpChannel(unsigned int opChannel) const
  {
    unsigned int detectorNum = (unsigned int) opChannel / fChannelsPerOpDet;
    return detectorNum;
  }

  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::HardwareChannelFromOpChannel(unsigned int opChannel) const
  {
    unsigned int channel = opChannel % fChannelsPerOpDet;
    return channel;
  }

  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::NTPCsets
    (readout::CryostatID const& cryoid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::MaxTPCsets() const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  bool DuneApaChannelMapAlg::HasTPCset(readout::TPCsetID const& tpcsetid) const
  {
    return tpcsetid.TPCset < NTPCsets(tpcsetid);
  }
  
  
  //----------------------------------------------------------------------------
  readout::TPCsetID DuneApaChannelMapAlg::TPCtoTPCset
    (geo::TPCID const& tpcid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  std::vector<geo::TPCID> DuneApaChannelMapAlg::TPCsetToTPCs
    (readout::TPCsetID const& tpcsetid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  geo::TPCID DuneApaChannelMapAlg::FirstTPCinTPCset
    (readout::TPCsetID const& tpcsetid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::NROPs
      (readout::TPCsetID const& tpcsetid) const
  {
    throw cet::exception("ChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  unsigned int DuneApaChannelMapAlg::MaxROPs() const {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  bool DuneApaChannelMapAlg::HasROP(readout::ROPID const& ropid) const {
    return ropid.ROP < NROPs(ropid);
  }
  
  
  //----------------------------------------------------------------------------
  readout::ROPID DuneApaChannelMapAlg::WirePlaneToROP
    (geo::PlaneID const& planeid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  std::vector<geo::PlaneID> DuneApaChannelMapAlg::ROPtoWirePlanes
    (readout::ROPID const& ropid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  std::vector<geo::TPCID> DuneApaChannelMapAlg::ROPtoTPCs
    (readout::ROPID const& ropid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  readout::ROPID DuneApaChannelMapAlg::ChannelToROP
    (raw::ChannelID_t channel) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  raw::ChannelID_t DuneApaChannelMapAlg::FirstChannelInROP
    (readout::ROPID const& ropid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  geo::PlaneID DuneApaChannelMapAlg::FirstWirePlaneInROP
    (readout::ROPID const& ropid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
} // namespace
