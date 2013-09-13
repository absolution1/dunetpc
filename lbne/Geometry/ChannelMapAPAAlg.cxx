////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapAPAAlg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////

#include "Geometry/ChannelMapAPAAlg.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace geo{

  //----------------------------------------------------------------------------
  ChannelMapAPAAlg::ChannelMapAPAAlg(fhicl::ParameterSet const& p)
    : fSorter(geo::GeoObjectSorterAPA(p))
  {
  }

  //----------------------------------------------------------------------------
  ChannelMapAPAAlg::~ChannelMapAPAAlg()
  {
  }

  //----------------------------------------------------------------------------
  void ChannelMapAPAAlg::Initialize(std::vector<geo::CryostatGeo*> & cgeo)
  {

    if(!fFirstChannelInThisPlane.empty() || !fFirstChannelInNextPlane.empty())
      {
	this->Uninitialize();
      }

    fNcryostat = cgeo.size();
    
    mf::LogInfo("ChannelMapAPAAlg") << "Sorting volumes...";

    fSorter.SortCryostats(cgeo);
    for(size_t c = 0; c < cgeo.size(); ++c) 
      cgeo[c]->SortSubVolumes(fSorter);

    mf::LogInfo("ChannelMapAPAAlg") << "Initializing channel map...";
      
    fNTPC.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(1);  // Change 1 to Ncryostat if you want
    fFirstChannelInThisPlane.resize(1);  // to treat each APA uniquely,and do
					 // the same with the other resizes.
    fPlanesPerAPA = cgeo[0]->TPC(0).Nplanes();
    nAnchoredWires.resize(fPlanesPerAPA);
    fWiresInPlane.resize(fPlanesPerAPA);
    fFirstChannelInThisPlane[0].resize(1);  // remember FirstChannel vectors
    fFirstChannelInNextPlane[0].resize(1);  // for first APA only.
    fFirstChannelInThisPlane[0][0].resize(fPlanesPerAPA);  // Make room for info
    fFirstChannelInNextPlane[0][0].resize(fPlanesPerAPA);  // on each plane.

    fTopChannel = 0;

    // Size some vectors and initialize the FirstChannel vectors.
    // If making FirstChannel's for every APA uniquely, they also
    // need to be sized here. Not necessary for now
    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      
      fNTPC[cs] = cgeo[cs]->NTPC();

    }// end sizing loop over cryostats

    // Find the number of wires anchored to the frame
    for(unsigned int p=0; p!=fPlanesPerAPA; ++p){

      fWiresInPlane[p] = cgeo[0]->TPC(0).Plane(p).Nwires();
      double xyz[3] = {0.};
      double xyz_next[3] = {0.};

      for(unsigned int w=0; w!=fWiresInPlane[p]; ++w){

	// for vertical planes
	if(cgeo[0]->TPC(0).Plane(p).View()==geo::kZ)   { 
	  nAnchoredWires[p] = fWiresInPlane[p];      
	  break;
	}

	cgeo[0]->TPC(0).Plane(p).Wire(w).GetCenter(xyz);
	cgeo[0]->TPC(0).Plane(p).Wire(w+1).GetCenter(xyz_next);

	if(xyz[2]==xyz_next[2]){
	  nAnchoredWires[p] = w-1;      
	  break;
	}
      }// end wire loop

    }// end plane loop

    static uint32_t CurrentChannel = 0;
   
    for(unsigned int PCount = 0; PCount != fPlanesPerAPA; ++PCount){

      fFirstChannelInThisPlane[0][0][PCount] = CurrentChannel;
      CurrentChannel = CurrentChannel + 2*nAnchoredWires[PCount];
      fFirstChannelInNextPlane[0][0][PCount] = CurrentChannel;

    }// end build loop over planes

    // Save the number of channels
    fChannelsPerAPA = fFirstChannelInNextPlane[0][0][fPlanesPerAPA-1];

    fNchannels = 0;
    for(size_t cs = 0; cs < fNcryostat; ++cs){
      fNchannels = fNchannels + fChannelsPerAPA*fNTPC[cs]/2;
    }

    //resize vectors
    fFirstWireCenterY.resize(fNcryostat);
    fFirstWireCenterZ.resize(fNcryostat);
    for (unsigned int cs=0; cs<fNcryostat; cs++){
      fFirstWireCenterY[cs].resize(cgeo[cs]->NTPC());
      fFirstWireCenterZ[cs].resize(cgeo[cs]->NTPC());
      for (unsigned int tpc=0; tpc<cgeo[cs]->NTPC(); tpc++){
        fFirstWireCenterY[cs][tpc].resize(cgeo[cs]->TPC(tpc).Nplanes());
        fFirstWireCenterZ[cs][tpc].resize(cgeo[cs]->TPC(tpc).Nplanes());
      }                                                                   
    }

    fWirePitch.resize(cgeo[0]->TPC(0).Nplanes());
    fOrientation.resize(cgeo[0]->TPC(0).Nplanes());
    fTanOrientation.resize(cgeo[0]->TPC(0).Nplanes());
    fCosOrientation.resize(cgeo[0]->TPC(0).Nplanes());

    //save data into fFirstWireCenterY and fFirstWireCenterZ
    for (unsigned int cs=0; cs<fNcryostat; cs++){
      for (unsigned int tpc=0; tpc<cgeo[cs]->NTPC(); tpc++){
        for (unsigned int plane=0; plane<cgeo[cs]->TPC(tpc).Nplanes(); plane++){
          double xyz[3]={0.0, 0.0, 0.0};
          cgeo[cs]->TPC(tpc).Plane(plane).Wire(0).GetCenter(xyz);
          fFirstWireCenterY[cs][tpc][plane]=xyz[1];
          fFirstWireCenterZ[cs][tpc][plane]=xyz[2];
        }
      }
    }

    //initialize fWirePitch and fOrientation
    for (unsigned int plane=0; plane<cgeo[0]->TPC(0).Nplanes(); plane++){
      fWirePitch[plane]=cgeo[0]->TPC(0).WirePitch(0,1,plane);
      fOrientation[plane]=cgeo[0]->TPC(0).Plane(plane).Wire(0).ThetaZ();
      fTanOrientation[plane] = tan(fOrientation[plane]);
      fCosOrientation[plane] = cos(fOrientation[plane]);
    }


    mf::LogVerbatim("GeometryTest") << "fNchannels = " << fNchannels ; 

    mf::LogVerbatim("GeometryTest") << "For all identical APA:" ; 
    mf::LogVerbatim("GeometryTest") << "Number of channels per APA = " << fChannelsPerAPA ; 

    mf::LogVerbatim("GeometryTest") << "Number of WireIDs in a U plane = " << fWiresInPlane[0] ;
    mf::LogVerbatim("GeometryTest") << "Number of WireIDs in a V plane = " << fWiresInPlane[1] ;
    mf::LogVerbatim("GeometryTest") << "Number of WireIDs in a Z Plane = " << fWiresInPlane[2] ;

    mf::LogVerbatim("GeometryTest") << "U channels per APA = " << 2*nAnchoredWires[0] ;
    mf::LogVerbatim("GeometryTest") << "V channels per APA = " << 2*nAnchoredWires[1] ;
    mf::LogVerbatim("GeometryTest") << "Z channels per APA side = " << nAnchoredWires[2] ;

    mf::LogVerbatim("GeometryTest") << "Pitch in U Plane = " << fWirePitch[0] ;
    mf::LogVerbatim("GeometryTest") << "Pitch in V Plane = " << fWirePitch[1] ;
    mf::LogVerbatim("GeometryTest") << "Pitch in Z Plane = " << fWirePitch[2] ;

    return;

  }
   
  //----------------------------------------------------------------------------
  void ChannelMapAPAAlg::Uninitialize()
  {

    std::vector< std::vector<std::vector<uint32_t> > >().swap(fFirstChannelInThisPlane);
    std::vector< std::vector<std::vector<uint32_t> > >().swap(fFirstChannelInNextPlane);

  }

  //----------------------------------------------------------------------------
  std::vector<geo::WireID> ChannelMapAPAAlg::ChannelToWire(uint32_t channel)  const
  {

    // first check if this channel ID is legal
    if(channel >= fNchannels )
      throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel;

    std::vector< WireID > AllSegments;
    
    static unsigned int cstat;
    static unsigned int tpc;
    static unsigned int plane;
    static unsigned int wireThisPlane;
    static unsigned int NextPlane;
    static unsigned int ThisPlane;
    
    uint32_t chan       = channel%fChannelsPerAPA;
    uint32_t pureAPAnum = std::floor( channel/fChannelsPerAPA );

    bool breakVariable = false;
    for(unsigned int planeloop = 0; planeloop != fPlanesPerAPA; ++planeloop){
	  
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
    uint32_t ChannelGroup = std::floor( wireThisPlane/nAnchoredWires[plane] );
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
  uint32_t ChannelMapAPAAlg::Nchannels() const
  {
    return fNchannels;
  }
  

  //----------------------------------------------------------------------------
  WireID  ChannelMapAPAAlg::NearestWireID(const TVector3& xyz,
					  unsigned int    plane,
					  unsigned int    tpc,
					  unsigned int    cryostat)     const
  {

    //get the position of first wire in a given cryostat, tpc and plane
    double firstxyz[3]={0.0, 0.0, 0.0};
    firstxyz[1]=fFirstWireCenterY[cryostat][tpc][plane];
    firstxyz[2]=fFirstWireCenterZ[cryostat][tpc][plane];

    //get the orientation angle of a given plane and calculate the distance between first wire
    //and a point projected in the plane
    int rotate = 1;
    if (tpc%2 == 1) rotate = -1;

    // old distance formula
    //double distance = std::abs(xyz[1]-firstxyz[1]-rotate*tan(fOrientation[plane])*xyz[2]
    //		   +   rotate*fTanOrientation[plane]*firstxyz[2])/
    //                         std::sqrt(fTanOrientation[plane]*fTanOrientation[plane]+1);

    // simplify and make faster

    double distance = std::abs( (xyz[1]-firstxyz[1] -rotate*fTanOrientation[plane]*(xyz[2]-firstxyz[2]))
                                * fCosOrientation[plane]);
    
    //by dividing distance by wirepitch and given that wires are sorted in increasing order,
    //then the wire that is closest to a given point can be calculated
    //uint32_t iwire=int(distance/fWirePitch[plane]);
    //if the distance between the wire and a given point is greater than the half of wirepitch,
    //then the point is closer to a i+1 wire thus add one
    //double res = distance/fWirePitch[plane] - int( distance/fWirePitch[plane] );
    //if (res > fWirePitch[plane]/2)	iwire+=1;

    double dwire=distance/fWirePitch[plane];
    uint32_t iwire=int(dwire);
    if (dwire-iwire>fWirePitch[plane]*0.5) ++iwire;
    uint32_t maxwireminus1=fWiresInPlane[plane]-1;
    if(iwire>maxwireminus1) iwire=maxwireminus1;

    WireID wid(cryostat, tpc, plane, iwire);
    return wid;

  }
  
  //----------------------------------------------------------------------------
  uint32_t ChannelMapAPAAlg::PlaneWireToChannel(unsigned int plane,
						unsigned int wire,
						unsigned int tpc,
						unsigned int cstat) const
  {
    unsigned int OtherSideWires = 0;

    uint32_t Channel = fFirstChannelInThisPlane[0][0][plane]; // start in very first APA.
    Channel += cstat*(fNTPC[cstat]/2)*fChannelsPerAPA;       // move channel to proper cstat.
    Channel += std::floor( tpc/2 )*fChannelsPerAPA;		  // move channel to proper APA.
    OtherSideWires += (tpc%2)*nAnchoredWires[plane];	          // get number of wires on the first
								  // side of the APA if starting
								  // on the other side TPC.


    // Lastly, account for the fact that channel number while moving up wire number in one
    // plane resets after 2 times the number of wires anchored -- one for each APA side.
    // At the same time, OtherSideWires accounts for the fact that if a channel starts on 
    // the other side, it is offset by the number of wires on the first side.
    Channel += (OtherSideWires + wire)%(2*nAnchoredWires[plane]);
    
    return Channel;

  }

  //----------------------------------------------------------------------------
  const SigType_t ChannelMapAPAAlg::SignalType( uint32_t const channel )  const
  {
    uint32_t chan = channel % fChannelsPerAPA;
    SigType_t sigt;

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
  const View_t ChannelMapAPAAlg::View( uint32_t const channel )  const
  {
    uint32_t chan = channel % fChannelsPerAPA;
    View_t view;

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


} // namespace
