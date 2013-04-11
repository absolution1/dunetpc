////////////////////////////////////////////////////////////////////////
/// \file  ChannelMap35Alg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////

#include "Geometry/ChannelMap35Alg.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace geo{



  //----------------------------------------------------------------------------
  // Define sort order for cryostats in APA configuration
  //   same as standard
  static bool sortCryo35(const CryostatGeo* c1, const CryostatGeo* c2)
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.}; 
    c1->LocalToWorld(local, xyz1);
    c2->LocalToWorld(local, xyz2);

    return xyz1[0] < xyz2[0];   
  }


  //----------------------------------------------------------------------------
  // Define sort order for tpcs in APA configuration.
  static bool sortTPC35(const TPCGeo* t1, const TPCGeo* t2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    // The goal is to number TPCs first in the x direction so that,
    // in the case of APA configuration, TPCs 2c and 2c+1 make up APA c.
    // then numbering will go in y then in z direction.

    // First sort all TPCs into same-z groups
    if(xyz1[2] < xyz2[2]) return true;
 
    // Within a same-z group, sort TPCs into same-y groups
    if(xyz1[2] == xyz2[2] && xyz1[1] < xyz2[1]) return true;
 
    // Within a same-z, same-y group, sort TPCs according to x
    if(xyz1[2] == xyz2[2] && xyz1[1] == xyz2[1] && xyz1[0] < xyz2[0]) return true;      

    // none of those are true, so return false
    return false;
  }


  //----------------------------------------------------------------------------
  // Define sort order for planes in APA configuration
  //   same as standard, but implemented differently
  static bool sortPlane35(const PlaneGeo* p1, const PlaneGeo* p2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    p1->LocalToWorld(local, xyz1);
    p2->LocalToWorld(local, xyz2);

    return xyz1[0] > xyz2[0];
  }

  //----------------------------------------------------------------------------
    // we want the wires to be sorted such that the smallest corner wire
    // on the readout end of a plane is wire zero, with wire number
    // increasing away from that wire.

    // Since 35t has an APA which is both above and below the world origin,
    // we cannot use the APA trick. If we could ask where wire 0 was, we could
    // still do this in a single implimentation, but we aren't sure what wire
    // center we will be getting, so this reversed sorting must be handled
    // at the plane level where there is one vertical center.
    // If the plane center is above, count from top down (the top stacked and
    // largest APAs) If the plane is below (bottom stacked APA) count bottom up

  bool sortWire35(WireGeo* w1, WireGeo* w2){
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};

    fflush(stdout);
    w1->GetCenter(xyz1); w2->GetCenter(xyz2);

    // immedieately take care of vertical wires regardless of which TPC
    // vertical wires should always have same y, and always increase in z direction
    if( xyz1[1]==xyz2[1] && xyz1[2]<xyz2[2] ) return true;


    // we want the wires to be sorted such that the smallest corner wire
    // on the readout end of a plane is wire zero, with wire number
    // increasing away from that wire. 

    // The number 76.35 is hard-coded as the z position of the vertical boundary 
    // between adjacent APAs. This is necessary because the largest TPC in 
    // 35t geometry would be considered "top" AND "bottom" by the APA 
    // configuration sorting method. Reconciling this came down to complicating 
    // the detector-agnostic framework of sorting volumes, or hard coding.
    // Let's keep the 35t messiness inside this class as much as possible.

    // the hard coded nuber should be equivalent to the zposition of 
    // volDetEnclosure in generate_35.pl

    if( xyz1[2]>76.35 ){

      // if a bottom TPC, count from bottom up
      if( xyz1[1]<0 && xyz1[1] < xyz2[1] ) return true; 
     
      // if a top TPC, count from top down
      if( xyz1[1]>0 && xyz1[1] > xyz2[1] ) return true;
      
    } else if( xyz1[2]<76.35 ){

      if( xyz1[1] > xyz2[1] ) return true;

    }

    return false;
  }


  //----------------------------------------------------------------------------
  ChannelMap35Alg::ChannelMap35Alg()
  {
  }

  //----------------------------------------------------------------------------
  ChannelMap35Alg::~ChannelMap35Alg()
  {
  }

  //----------------------------------------------------------------------------
  void ChannelMap35Alg::Initialize(std::vector<geo::CryostatGeo*> & cgeo)
  {

    if(!fFirstChannelInThisPlane.empty() || !fFirstChannelInNextPlane.empty())
      {
	this->Uninitialize();
      }


    fNcryostat = cgeo.size();
    
    mf::LogInfo("ChannelMap35Alg") << "Sorting...";

    std::sort(cgeo.begin(), cgeo.end(), sortCryo35);
    for(size_t c = 0; c < cgeo.size(); ++c) 
      cgeo[c]->SortSubVolumes(sortTPC35, sortPlane35, sortWire35);


    mf::LogInfo("ChannelMap35Alg") << "Initializing...";
      
    fNTPC.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    nAnchoredWires.resize(fNcryostat);

    fPlanesPerAPA = cgeo[0]->TPC(0).Nplanes();

    fTopChannel = 0;

    // Size some vectors and initialize the FirstChannel vectors.
    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      
      fNTPC[cs] = cgeo[cs]->NTPC();

      nAnchoredWires[cs].resize(fNTPC[cs]);
      fWiresPerPlane[cs].resize(fNTPC[cs]);
      fFirstChannelInThisPlane[cs].resize(fNTPC[cs]/2);
      fFirstChannelInNextPlane[cs].resize(fNTPC[cs]/2);

      for(unsigned int apa = 0; apa != fNTPC[cs]/2; ++apa){
	
        nAnchoredWires[cs][apa].resize(fPlanesPerAPA);
	fWiresPerPlane[cs][apa].resize(fPlanesPerAPA);
        fFirstChannelInThisPlane[cs][apa].resize(fPlanesPerAPA);
        fFirstChannelInNextPlane[cs][apa].resize(fPlanesPerAPA);

      }// end loop over apas
    }// end cryostats

    // Find the number of wires anchored to the frame
    for(unsigned int c = 0; c != fNcryostat; ++c){
      for(unsigned int a = 0; a != fNTPC[c]/2; ++a){
        for(unsigned int p = 0; p != fPlanesPerAPA; ++p){

          unsigned int t = 2*a;
          fWiresPerPlane[c][a][p] = cgeo[c]->TPC(t).Plane(p).Nwires();
          double xyz[3] = {0.};
          double xyz_next[3] = {0.};

          for(unsigned int w=0; w!=fWiresPerPlane[c][a][p]; ++w){

	    // for vertical planes
	    if(cgeo[c]->TPC(t).Plane(p).View() == geo::kZ)   { 
	      nAnchoredWires[c][a][p] = fWiresPerPlane[c][a][p];      
	      break;
	    }

	    cgeo[c]->TPC(t).Plane(p).Wire(w).GetCenter(xyz);
	    cgeo[c]->TPC(t).Plane(p).Wire(w+1).GetCenter(xyz_next);

    	    if(xyz[2]==xyz_next[2]){
	      nAnchoredWires[c][a][p] = w-1;      
	      break;
	    }

          }
        }
      }
    }

    static uint32_t CurrentChannel = 0;
 
    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      for(unsigned int apa = 0; apa != fNTPC[cs]/2; ++apa){  
        for(unsigned int p = 0; p != fPlanesPerAPA; ++p){

          fFirstChannelInThisPlane[cs][apa][p] = CurrentChannel;
          CurrentChannel = CurrentChannel + 2*nAnchoredWires[cs][apa][p];
          fFirstChannelInNextPlane[cs][apa][p] = CurrentChannel;

        }// end plane loop
      }// end apa loop
    }// end cs


    // Save the number of channels
    fNchannels = CurrentChannel;

    // Save the number of channels
    fChannelsPerAPA = fFirstChannelInNextPlane[0][0][fPlanesPerAPA-1];

    //resize vectors
    fFirstWireCenterY.resize(fNcryostat);
    fFirstWireCenterZ.resize(fNcryostat);
    for (unsigned int cs=0; cs<fNcryostat; cs++){
      fFirstWireCenterY[cs].resize(fNTPC[cs]);
      fFirstWireCenterZ[cs].resize(fNTPC[cs]);
      for (unsigned int tpc=0; tpc<fNTPC[cs]; tpc++){
        fFirstWireCenterY[cs][tpc].resize(fPlanesPerAPA);
        fFirstWireCenterZ[cs][tpc].resize(fPlanesPerAPA);
      }                                                                   
    }

    fWirePitch.resize(fPlanesPerAPA);
    fOrientation.resize(fPlanesPerAPA);
    fTanOrientation.resize(fPlanesPerAPA);
    fCosOrientation.resize(fPlanesPerAPA);


    //save data into fFirstWireCenterY and fFirstWireCenterZ
    for (unsigned int cs=0; cs<fNcryostat; cs++){
      for (unsigned int tpc=0; tpc<fNTPC[cs]; tpc++){
        for (unsigned int plane=0; plane<fPlanesPerAPA; plane++){
          double xyz[3]={0.0, 0.0, 0.0};
          cgeo[cs]->TPC(tpc).Plane(plane).Wire(0).GetCenter(xyz);
          fFirstWireCenterY[cs][tpc][plane]=xyz[1];
          fFirstWireCenterZ[cs][tpc][plane]=xyz[2];
        }
      }
    }

    //initialize fWirePitch and fOrientation
    for (unsigned int plane=0; plane<fPlanesPerAPA; plane++){
      fWirePitch[plane]=cgeo[0]->TPC(0).WirePitch(0,1,plane);
      fOrientation[plane]=cgeo[0]->TPC(0).Plane(plane).Wire(0).ThetaZ();
      fTanOrientation[plane] = tan(fOrientation[plane]);
      fCosOrientation[plane] = cos(fOrientation[plane]);
    }


    mf::LogVerbatim("GeometryTest") << "fNchannels = " << fNchannels ; 
    mf::LogVerbatim("GeometryTest") << "nAnchoredWires, plane 0 = " << nAnchoredWires[0][0][0] ;
    mf::LogVerbatim("GeometryTest") << "nAnchoredWires, plane 1 = " << nAnchoredWires[0][0][1] ;
    mf::LogVerbatim("GeometryTest") << "nAnchoredWires, plane 2 = " << nAnchoredWires[0][0][2] ;

    return;

  }
   
  //----------------------------------------------------------------------------
  void ChannelMap35Alg::Uninitialize()
  {

    std::vector< std::vector<std::vector<uint32_t> > >().swap(fFirstChannelInThisPlane);
    std::vector< std::vector<std::vector<uint32_t> > >().swap(fFirstChannelInNextPlane);

  }

  //----------------------------------------------------------------------------
  std::vector<geo::WireID> ChannelMap35Alg::ChannelToWire(uint32_t channel)  const
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
    
    for(unsigned int csloop = 0; csloop != fNcryostat; ++csloop){
      
      bool breakVariable = false;
      
      for(unsigned int apaloop = 0; apaloop != fNTPC[csloop]/2; ++apaloop){
	for(unsigned int planeloop = 0; planeloop != fPlanesPerAPA; ++planeloop){
	  
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
    uint32_t ChannelGroup = std::floor( wireThisPlane/nAnchoredWires[cstat][tpc/2][plane] );
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
  uint32_t ChannelMap35Alg::Nchannels() const
  {
    return fNchannels;
  }
  

  //----------------------------------------------------------------------------
  uint32_t    ChannelMap35Alg::NearestWire(const TVector3& xyz,
					   unsigned int    plane,
					   unsigned int    tpc,
					   unsigned int    cryostat)     const
  {

    //get the position of first wire in a given cryostat, tpc and plane
    double firstxyz[3]={0.0, 0.0, 0.0};
    firstxyz[1]=fFirstWireCenterY[cryostat][tpc][plane];
    firstxyz[2]=fFirstWireCenterZ[cryostat][tpc][plane];

    double distance = 0.;

    if (plane==2){  
      distance = xyz[2] - firstxyz[2];
    } else {

      //get the orientation angle of a given plane and calculate the distance between first wire
      //and a point projected in the plane
      int rotate = 1;
      if (tpc%2 == 1) rotate = -1;
 
      distance = std::abs( (xyz[1]-firstxyz[1] -rotate*fTanOrientation[plane]*(xyz[2]-firstxyz[2]))
                                * fCosOrientation[plane]);

    }

    //if the distance between the wire and a given point is greater than the half of wirepitch,
    //then the point is closer to a i+1 wire thus add one
    //double res = distance/fWirePitch[plane] - int( distance/fWirePitch[plane] );
    //if (res > fWirePitch[plane]/2)	iwire+=1;

    // do it, but also check to see if we are on the edge

    double dwire=distance/fWirePitch[plane];
    uint32_t iwire=int(dwire);
    if (dwire-iwire>fWirePitch[plane]*0.5) ++iwire;
    uint32_t maxwireminus1=fWiresPerPlane[0][tpc/2][plane]-1;
    if(iwire>maxwireminus1) iwire=maxwireminus1;

    return iwire;

  }
  
  //----------------------------------------------------------------------------
  uint32_t ChannelMap35Alg::PlaneWireToChannel(unsigned int plane,
					       unsigned int wire,
					       unsigned int tpc,
					       unsigned int cstat) const
  {

    unsigned int OtherSideWires = 0;

    uint32_t Channel = fFirstChannelInThisPlane[cstat][std::floor(tpc/2)][plane];

    // get number of wires starting on the first side of the APA if starting
    // on the other side TPC.
    OtherSideWires += (tpc%2)*nAnchoredWires[cstat][std::floor(tpc/2)][plane];
    
    // Lastly, account for the fact that channel number while moving up wire number in one
    // plane resets after 2 times the number of wires anchored -- one for each APA side.
    // At the same time, OtherSideWires accounts for the fact that if a channel starts on
    // the other side, it is offset by the number of wires on the first side.
    Channel += (OtherSideWires + wire)%(2*nAnchoredWires[cstat][std::floor(tpc/2)][plane]);

    return Channel;

  }


  //----------------------------------------------------------------------------
  const SigType_t ChannelMap35Alg::SignalType( uint32_t const channel )  const
  {
    uint32_t chan = channel % fChannelsPerAPA;
    SigType_t sigt;

    if(       chan <  fFirstChannelInThisPlane[0][0][2]     ){ sigt = kInduction;  }
    else if( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
             (chan <  fFirstChannelInNextPlane[0][0][2])    ){ sigt = kCollection; }
    else{    mf::LogWarning("BadChannelSignalType") << "Channel " << channel 
						    << " (" << chan << ") not given signal type." << std::endl;         }
  
    return sigt;
  }

  //----------------------------------------------------------------------------
  const View_t ChannelMap35Alg::View( uint32_t const channel )  const
  {
    uint32_t chan = channel % fChannelsPerAPA;
    View_t view;

    if(       chan <  fFirstChannelInNextPlane[0][0][0]     ){ view = geo::kU; }
    else if( (chan >= fFirstChannelInThisPlane[0][0][1]) &&
             (chan <  fFirstChannelInNextPlane[0][0][1])    ){ view = geo::kV; }
    else if( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
             (chan <  fFirstChannelInNextPlane[0][0][2])    ){ view = geo::kZ; }
    else{    mf::LogWarning("BadChannelViewType") << "Channel " << channel 
						  << " (" << chan << ") not given view type.";}
    
    return view;
  }  
 


} // namespace
