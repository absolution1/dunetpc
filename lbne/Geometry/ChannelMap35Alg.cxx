////////////////////////////////////////////////////////////////////////
/// \file  ChannelMap35Alg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////

#include "lbne/Geometry/ChannelMap35Alg.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace geo{

  //----------------------------------------------------------------------------
  ChannelMap35Alg::ChannelMap35Alg(fhicl::ParameterSet const& p)
    : fSorter(geo::GeoObjectSorter35(p))
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

    fSorter.SortCryostats(cgeo);
    for(size_t c = 0; c < cgeo.size(); ++c) 
      cgeo[c]->SortSubVolumes(fSorter);


    mf::LogInfo("ChannelMap35Alg") << "Initializing...";
      
    fNTPC.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    nAnchoredWires.resize(fNcryostat);
    fViews.clear();
    fPlaneIDs.clear();
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

	  fViews.emplace(cgeo[c]->TPC(t).Plane(p).View());

          for(unsigned int w = 0; w != fWiresPerPlane[c][a][p]; ++w){

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

    
    fWirePitch.resize(fPlanesPerAPA);
    fOrientation.resize(fPlanesPerAPA);
    fSinOrientation.resize(fPlanesPerAPA);
    fCosOrientation.resize(fPlanesPerAPA);


    //save data into fFirstWireCenterY, fFirstWireCenterZ and fWireSortingInZ
    fPlaneData.resize(fNcryostat);
    for (unsigned int cs=0; cs<fNcryostat; cs++){
      fPlaneData[cs].resize(fNTPC[cs]);
      for (unsigned int tpc=0; tpc<fNTPC[cs]; tpc++){
        fPlaneData[cs][tpc].resize(fPlanesPerAPA);
        for (unsigned int plane=0; plane<fPlanesPerAPA; plane++){
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
        } // for plane
      } // for TPC
    } // for cryostat

    //initialize fWirePitch and fOrientation
    for (unsigned int plane=0; plane<fPlanesPerAPA; plane++){
      fWirePitch[plane]=cgeo[0]->TPC(0).WirePitch(0,1,plane);
      fOrientation[plane]=cgeo[0]->TPC(0).Plane(plane).Wire(0).ThetaZ();
      fSinOrientation[plane] = sin(fOrientation[plane]);
      fCosOrientation[plane] = cos(fOrientation[plane]);

    } // for
    
    
    mf::LogVerbatim("GeometryTest") << "fNchannels = " << fNchannels ; 
    mf::LogVerbatim("GeometryTest") << "U channels per APA = " << 2*nAnchoredWires[0][0][0] ;
    mf::LogVerbatim("GeometryTest") << "V channels per APA = " << 2*nAnchoredWires[0][0][1] ;
    mf::LogVerbatim("GeometryTest") << "Z channels per APA side = " << nAnchoredWires[0][0][2] ;

    mf::LogVerbatim("GeometryTest") << "Pitch in U Plane = " << fWirePitch[0] ;
    mf::LogVerbatim("GeometryTest") << "Pitch in V Plane = " << fWirePitch[1] ;
    mf::LogVerbatim("GeometryTest") << "Pitch in Z Plane = " << fWirePitch[2] ;

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
  double ChannelMap35Alg::WireCoordinate(double YPos, double ZPos,
                                         unsigned int PlaneNo,
                                         unsigned int TPCNo,
                                         unsigned int cstat) const
  {
    // Returns the wire number corresponding to a (Y,Z) position in PlaneNo
    // with float precision.
    // Core code ripped from original NearestWireID() implementation
    
    const PlaneData_t& PlaneData = fPlaneData[cstat][TPCNo][PlaneNo];
    
    //get the orientation angle of a given plane and calculate the distance between first wire
    //and a point projected in the plane
//    const double rotate = (TPCNo % 2 == 1)? -1.: +1.;
    
    // the formula used here is geometric:
    // distance = delta_y cos(theta_z) + delta_z sin(theta_z)
    // with a correction for the orientation of the TPC:
    // odd TPCs have supplementary wire angle (pi-theta_z), changing cosine sign
    
    const bool bSuppl = (TPCNo % 2) == 1;
    float distance =
      -(YPos - PlaneData.fFirstWireCenterY) * (bSuppl? -1.: +1.) * fCosOrientation[PlaneNo]
      +(ZPos - PlaneData.fFirstWireCenterZ) * fSinOrientation[PlaneNo]
      ;
    
    // The sign of this formula is correct if the wire with larger ID is on the
    // "right" (intuitively, larger z; rigorously, smaller intercept)
    // than this one.
    // Of course, we are not always that lucky. fWireSortingInZ fixes our luck.

    return PlaneData.fWireSortingInZ * distance/fWirePitch[PlaneNo];
  } // ChannelMap35Alg::WireCoordinate()
  
  
  //----------------------------------------------------------------------------
  WireID ChannelMap35Alg::NearestWireID(const TVector3& xyz,
                                        unsigned int    plane,
                                        unsigned int    tpc,
                                        unsigned int    cryostat)     const
  {
    // add 0.5 to have the correct rounding
    int NearestWireNumber
      = int (0.5 + WireCoordinate(xyz.Y(), xyz.Z(), plane, tpc, cryostat));
    
    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if(NearestWireNumber < 0 ||
       NearestWireNumber >= (int) fWiresPerPlane[cryostat][tpc/2][plane])
    {
      const int wireNumber = NearestWireNumber; // save for the output
      
      if(wireNumber < 0 ) NearestWireNumber = 0;
      else                NearestWireNumber = fWiresPerPlane[cryostat][tpc/2][plane] - 1;
    
    /*
      // comment in the following statement to throw an exception instead
      throw InvalidWireIDError("Geometry", wireNumber, NearestWireNumber)
        << "ChannelMap35Alg::NearestWireID(): can't Find Nearest Wire for position (" 
        << xyz.X() << "," << xyz.Y() << "," << xyz.Z() << ")"
        << " approx wire number # " << wireNumber
        << " (capped from " << NearestWireNumber << ")\n";
    */
    } // if invalid wire
    
    return { cryostat, tpc, plane, (unsigned int) NearestWireNumber };
  } // ChannelMap35Alg::NearestWireID()
  
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
  SigType_t ChannelMap35Alg::SignalType( uint32_t const channel )  const
  {
    uint32_t chan = channel % fChannelsPerAPA;
    SigType_t sigt = kInduction;

    if(       chan <  fFirstChannelInThisPlane[0][0][2]     ){ sigt = kInduction;  }
    else if( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
             (chan <  fFirstChannelInNextPlane[0][0][2])    ){ sigt = kCollection; }
    else{    mf::LogWarning("BadChannelSignalType") << "Channel " << channel 
						    << " (" << chan << ") not given signal type." << std::endl;         }
  
    return sigt;
  }

  //----------------------------------------------------------------------------
  View_t ChannelMap35Alg::View( uint32_t const channel )  const
  {
    uint32_t chan = channel % fChannelsPerAPA;
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
  std::set<View_t> const& ChannelMap35Alg::Views() const
  {
    return fViews;
  }

  //----------------------------------------------------------------------------
  std::set<PlaneID> const& ChannelMap35Alg::PlaneIDs() const
  {
    return fPlaneIDs;
  }

} // namespace
