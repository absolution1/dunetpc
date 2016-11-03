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

using std::vector;
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
    
  vector<geo::CryostatGeo*>& crygeos = geodata.cryostats;
  vector<geo::AuxDetGeo*>  & adgeo = geodata.auxDets;
  Index ncry = crygeos.size();
  fNcryostat = ncry;
  if ( ncry == 0 ) {
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
  fNTPC.resize(ncry);
  fNAPA.resize(ncry);
  // The first channel array for each TPC plane assume the cryostats are identical.
  fPlanesPerTPC = crygeo.TPC(0).Nplanes();
  fFirstChannelInThisPlane.resize(1);  // uniquely,and do // the same with the other resizes.
  fFirstChannelInNextPlane.resize(1);  // Change 1 to Ncryostat if you want to treat each APA
  fFirstChannelInThisPlane[0].resize(1);  // remember FirstChannel vectors
  fFirstChannelInNextPlane[0].resize(1);  // for first APA only.
  fFirstChannelInThisPlane[0][0].resize(fPlanesPerTPC);  // Make room for info
  fFirstChannelInNextPlane[0][0].resize(fPlanesPerTPC);  // on each plane.
  if ( fPlanesPerTPC != 3 ) {
    mf::LogError("DuneApaChannelMapAlg") << "# planes/TPC is not 3: " << fPlanesPerTPC;
    return;
  }
  // The first channel array for each APA plane allow the cryostats to differ.
  fFirstChannelInThisRop.resize(ncry);
  fFirstChannelInNextRop.resize(ncry);
  fRopsPerApa.resize(ncry);
  fChannelsPerApa.resize(ncry);
  fPlanesPerRop.resize(ncry);
  fRopTpc.resize(ncry);
  fRopPlane.resize(ncry);
  Index itpc = 0;
  for ( Index icry=0; icry<ncry; ++icry ) {
    Index ntpc = crygeos[icry]->NTPC();
    Index napa = ntpc/2;  // Assume 1 APA for every two TPCs
    fRopsPerApa[icry].resize(napa, 4);
    fChannelsPerApa[icry].resize(napa, 0);
    fPlanesPerRop[icry].resize(napa);
    fFirstChannelInThisRop[icry].resize(napa);
    fFirstChannelInNextRop[icry].resize(napa);
    fRopTpc[icry].resize(napa);
    fRopPlane[icry].resize(napa);
    for ( Index iapa=0; iapa<napa; ++iapa ) {
      Index nrop = fRopsPerApa[icry][iapa];
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
        fRopTpc[icry][iapa][irop].push_back(++itpc);
        fRopPlane[icry][iapa][irop].push_back(ipla);
      }
    }
  }
  nAnchoredWires.resize(fPlanesPerTPC);
  fWiresInPlane.resize(fPlanesPerTPC);
  fViews.clear();
  fPlaneIDs.clear();
  fTopChannel = 0;

  // Size some vectors and initialize the FirstChannel vectors.
  // If making FirstChannel's for every APA uniquely, they also
  // need to be sized here. Not necessary for now
  for ( unsigned int icry = 0; icry != ncry; ++icry ) {
    fNTPC[icry] = crygeos[icry]->NTPC();
    fNAPA[icry] = fNTPC[icry]/2;
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

  // Assign first channels for the TPCs.
  raw::ChannelID_t icha = 0;
  for ( Index ipla=0; ipla!=fPlanesPerTPC; ++ipla ) {
    fFirstChannelInThisPlane[0][0][ipla] = icha;
    icha += 2*nAnchoredWires[ipla];
    fFirstChannelInNextPlane[0][0][ipla] = icha;
  }
  Index nchaTpc = fFirstChannelInNextPlane[0][0][fPlanesPerTPC-1];
  fNchannels = 0;
  for ( size_t icry=0; icry<ncry; ++icry ) {
    fNchannels += nchaTpc*fNTPC[icry]/2;
  }
  fChannelsPerAPA = nchaTpc;

  // Assign first channels for the APAs.
  icha = 0;
  for ( Index icry=0; icry!=ncry; ++icry ) {
    Index napa = fNAPA[icry];
    for ( Index iapa=0; iapa!=napa; ++iapa ) {
      Index nrop = fRopsPerApa[icry][iapa];
      for ( Index irop=0; irop!=nrop; ++irop ) {
        Index nrpl = fPlanesPerRop[icry][iapa][irop];
        fFirstChannelInThisRop[icry][iapa][irop] = icha;
        for ( Index irpl=0; irpl!=nrpl; ++irpl ) {
          //Index itpc = fRopTpc[icry][iapa][irop][irpl];
          Index ipla = fRopPlane[icry][iapa][irop][irpl];
          Index nchaPlane = nAnchoredWires[ipla];
          icha += nchaPlane;
        }
        fFirstChannelInNextRop[icry][iapa][irop] = icha;
      }
    }
  }
  if ( icha != fNchannels ) {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << "TPC and APA channel counts disagree: "
                                                 << fNchannels << " != " << icha << endl;
  }

  fPlaneData.resize(ncry);
  for ( Index icry=0; icry<ncry; ++icry ) {
    fPlaneData[icry].resize(crygeos[icry]->NTPC());
    for ( Index itpc=0; itpc<crygeos[icry]->NTPC(); ++itpc ) {
      fPlaneData[icry][itpc].resize(crygeos[icry]->TPC(itpc).Nplanes());
      for ( Index ipla=0; ipla<crygeos[icry]->TPC(itpc).Nplanes(); ++ipla ) {
        PlaneData_t& PlaneData = fPlaneData[icry][itpc][ipla];
        const geo::PlaneGeo& thePlane = crygeos[icry]->TPC(itpc).Plane(ipla);
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
      }
    }
  }

  Index npla = crygeo.TPC(0).Nplanes();
  fWirePitch.resize(npla);
  fOrientation.resize(npla);
  fSinOrientation.resize(npla);
  fCosOrientation.resize(npla);

  for ( Index ipla=0; ipla<npla; ++ipla ) {
    fWirePitch[ipla] = crygeo.TPC(0).WirePitch(0,1,ipla);
    fOrientation[ipla] = crygeo.TPC(0).Plane(ipla).Wire(0).ThetaZ();
    fSinOrientation[ipla] = sin(fOrientation[ipla]);
    fCosOrientation[ipla] = cos(fOrientation[ipla]);
  }

  for ( Index icry=0; icry<ncry; ++icry ) {
    mf::LogVerbatim("DuneApaChannelMapAlg") << "Cryostat " << icry << ":"; 
    mf::LogVerbatim("DuneApaChannelMapAlg") << "  " << fNchannels << " total channels"; 
    mf::LogVerbatim("DuneApaChannelMapAlg") << "  " << fNTPC[icry]/2 << " APAs";       
    mf::LogVerbatim("DuneApaChannelMapAlg") << "  For all identical APA:" ; 
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Number of channels per APA = " << fChannelsPerAPA ; 
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    U channels per APA = " << 2*nAnchoredWires[0] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    V channels per APA = " << 2*nAnchoredWires[1] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Z channels per APA = " << 2*nAnchoredWires[2] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in U Plane = " << fWirePitch[0] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in V Plane = " << fWirePitch[1] ;
    mf::LogVerbatim("DuneApaChannelMapAlg") << "    Pitch in Z Plane = " << fWirePitch[2] ;
  }

}

//----------------------------------------------------------------------------

void DuneApaChannelMapAlg::Uninitialize() {
  PlaneInfoMap_t<raw::ChannelID_t>().swap(fFirstChannelInThisPlane);
  PlaneInfoMap_t<raw::ChannelID_t>().swap(fFirstChannelInNextPlane);
}

//----------------------------------------------------------------------------

vector<geo::WireID> DuneApaChannelMapAlg::ChannelToWire(raw::ChannelID_t icha) const {
  vector< WireID > wirids;
  if ( icha >= fNchannels ) return wirids;

  // Loop over planes to find the one holding this channel.
  // We assume all cryostats have the same (offset) wire-channel mapping.
  // If not, this must be replaced with a loop over cryostats.
  raw::ChannelID_t ichaApa = icha%fChannelsPerAPA;
  raw::ChannelID_t iapa = icha/fChannelsPerAPA;
  Index ntpc = fNTPC[0];
  Index napa = fNAPA[0];
  Index icry = icha/(napa*fChannelsPerAPA);
  Index itpc = (2*iapa) % ntpc;

  Index ipla = 0;
  Index iwir = 0;
  for ( ipla=0; ipla!=fPlanesPerTPC; ++ipla ) {
    Index icha1 = fFirstChannelInThisPlane[0][0][ipla];
    Index icha2 = fFirstChannelInNextPlane[0][0][ipla];
    iwir = ichaApa - icha1;
    if ( ichaApa >= icha1 && ichaApa < icha2 ) break;
  }
  if ( ipla >= fPlanesPerTPC ) {
    mf::LogError("DuneApaChannelMapAlg") << "Unable to find APA plane for channel " << icha;
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << ": Unable to find APA plane for channel " << icha;
    return wirids;
  }

  int WrapDirection = 1; // go from tpc to (tpc+1) or tpc to (tpc-1)

  // find the lowest wire
  raw::ChannelID_t ChannelGroup = iwir/nAnchoredWires[ipla];
  unsigned int bottomwire = iwir - ChannelGroup*nAnchoredWires[ipla];
    
  if(ChannelGroup%2==1){
    itpc += 1;
    WrapDirection  = -1;         
  }
    
  for(unsigned int SegCount = 0; SegCount != 50; ++SegCount){
      
    itpc += WrapDirection*(SegCount%2);
    geo::WireID CodeWire(icry, itpc, ipla, bottomwire + SegCount*nAnchoredWires[ipla]);
    wirids.push_back(CodeWire);
      
    // reset the tcp variable so it doesnt "accumulate value"
    itpc -= WrapDirection*(SegCount%2);
      
    if( bottomwire + (SegCount+1)*nAnchoredWires[ipla] > fWiresInPlane[ipla]-1) break;
  }
    
  return wirids;
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
  vector<geo::TPCID> DuneApaChannelMapAlg::TPCsetToTPCs
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
  vector<geo::PlaneID> DuneApaChannelMapAlg::ROPtoWirePlanes
    (readout::ROPID const& ropid) const
  {
    throw cet::exception("DuneApaChannelMapAlg") << __func__ << " not implemented!\n";
  }
  
  
  //----------------------------------------------------------------------------
  vector<geo::TPCID> DuneApaChannelMapAlg::ROPtoTPCs
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
