////////////////////////////////////////////////////////////////////////
// $Id: GeometryTest_module.cc,v 1.1 2011/02/17 01:45:48 brebel Exp $
//
//
// geometry unit tests
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef GEO_GEOMETRYTEST_H
#define GEO_GEOMETRYTEST_H
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

// ROOT includes
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TNtuple.h"
#include "TGeoManager.h"
#include "TStopwatch.h"
#include "TMath.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/OpDetGeo.h"
#include "Geometry/geo.h"
#include "SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDAnalyzer.h"

namespace geo { class Geometry; }

///tracking algorithms
namespace geo {
  class GeometryTest : public art::EDAnalyzer {
  public:
    explicit GeometryTest(fhicl::ParameterSet const& pset);
    virtual ~GeometryTest();

    virtual void analyze(art::Event const&) {}
    virtual void beginJob();

  private:

    void printChannelSummary();
    void printVolBounds();
    void printDetDim();
    void printWirePos();
    void printWiresInTPC(const TPCGeo& tpc, std::string indent = "") const;
    void printAllGeometry() const;
    void testCryostat();
    void testTPC(unsigned int const& c);
    void testChannelToWire();
    void testFindPlaneCenters();
    void testProject();
    void testWirePitch();
    void testPlanePitch();
    void testStandardWirePos();
    void testAPAWirePos();
    void testNearestWire();
    void testStepping();

    bool fCheckOverlaps;  ///< do the overlap check or not
    bool fPrintWires;  ///< print all the wires in geometry (really: all!)
  };
}

namespace geo{

  //......................................................................
  GeometryTest::GeometryTest(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset)
    , fCheckOverlaps( pset.get<bool>("CheckForOverlaps", false) )
    , fPrintWires( pset.get<bool>("PrintWires", false) )
  {
  }

  //......................................................................
  GeometryTest::~GeometryTest()
  {
  }

  //......................................................................
  void GeometryTest::beginJob()
  {
    art::ServiceHandle<geo::Geometry> geom;

    // change the printed version number when changing the output
    mf::LogVerbatim("GeometryTest") << "GeometryTest version 1.0";
    
    try{
      mf::LogVerbatim("GeometryTest") << "Wire Rmax  "         << geom->Plane(1).Wire(10).RMax()    ;
      mf::LogVerbatim("GeometryTest") << "Wire length "        << 2.*geom->Plane(1).Wire(10).HalfL();
      mf::LogVerbatim("GeometryTest") << "Wire Rmin  "         << geom->Plane(1).Wire(10).RMin()    ;
      mf::LogVerbatim("GeometryTest") << "Total mass "         << geom->TotalMass()                 ;
      mf::LogVerbatim("GeometryTest") << "Number of views "    << geom->Nviews()                    ;
      mf::LogVerbatim("GeometryTest") << "Number of channels " << geom->Nchannels()                 ;

      LOG_DEBUG("GeometryTest") << "print channel information ...";
      printChannelSummary();
      LOG_DEBUG("GeometryTest") << "done printing.";
      //mf::LogVerbatim("GeometryTest") << "print Cryo/TPC boundaries in world coordinates ...";
      //printVolBounds();
      //mf::LogVerbatim("GeometryTest") << "done printing.";
      //mf::LogVerbatim("GeometryTest") << "print Cryo/TPC dimensions ...";
      //printDetDim();
      //mf::LogVerbatim("GeometryTest") << "done printing.";
      //mf::LogVerbatim("GeometryTest") << "print wire center positions in world coordinates ...";
      //printWirePos();
      //mf::LogVerbatim("GeometryTest") << "done printing.";

      if(fCheckOverlaps){
        LOG_DEBUG("GeometryTest") << "test for overlaps ...";
        gGeoManager->CheckOverlaps(1e-5);
        gGeoManager->PrintOverlaps();
        LOG_DEBUG("GeometryTest") << "complete.";
      }

      LOG_DEBUG("GeometryTest") << "test Cryostat methods ...";
      testCryostat();
      LOG_DEBUG("GeometryTest") << "complete.";

      LOG_DEBUG("GeometryTest") << "test channel to plane wire and back ...";
      testChannelToWire();
      LOG_DEBUG("GeometryTest") << "complete.";

      LOG_DEBUG("GeometryTest") << "test find plane centers...";
      testFindPlaneCenters();
      LOG_DEBUG("GeometryTest") << "complete.";

      LOG_DEBUG("GeometryTest") << "testProject...";
      testProject();
      LOG_DEBUG("GeometryTest") << "complete.";

      //LOG_DEBUG("GeometryTest") << "testWirePos...";
      // There is a contradiction here, and these must be tested differently
      // Testing based on detector ID should NOT become common practice
      //if( geom->DetId()==geo::kLBNE10kt 
      //         || geom->DetId()==geo::kLBNE35t 
      //	   || geom->DetId()==geo::kLBNE34kt) testAPAWirePos();
      //else testStandardWirePos();
      //LOG_DEBUG("GeometryTest") << "complete.";

      LOG_DEBUG("GeometryTest") << "testNearestWire...";
      testNearestWire();
      LOG_DEBUG("GeometryTest") << "complete.";
      
      LOG_DEBUG("GeometryTest") << "testWirePitch...";
      testWirePitch();
      LOG_DEBUG("GeometryTest") << "complete.";

      LOG_DEBUG("GeometryTest") << "testPlanePitch...";
      testPlanePitch();
      LOG_DEBUG("GeometryTest") << "complete.";

      LOG_DEBUG("GeometryTest") << "testStepping...";
      testStepping();
      LOG_DEBUG("GeometryTest") << "complete.";
      
      if (fPrintWires) {
        LOG_DEBUG("GeometryTest") << "printAllGeometry...";
        printAllGeometry();
        LOG_DEBUG("GeometryTest") << "complete.";
      }
    }
    catch (cet::exception &e) {
      mf::LogWarning("GeometryTest") << "exception caught: \n" << e;
    }
    
    return;
  }


  //......................................................................
  void GeometryTest::printChannelSummary()
  {
    art::ServiceHandle<geo::Geometry> geom;
    
    static unsigned int OneSeg = 0;
    static unsigned int TwoSegs = 0;
    static unsigned int ThreeSegs = 0;
    static unsigned int FourSegs = 0;
    uint32_t channels = geom->Nchannels();
    if(geom->NTPC() > 1) channels /= geom->NTPC()/2;

    for(uint32_t c = 0; c < channels; c++){

      unsigned int ChanSize = geom->ChannelToWire(c).size();

       if     (ChanSize==1) ++OneSeg;
       else if(ChanSize==2) ++TwoSegs;
       else if(ChanSize==3) ++ThreeSegs;
       else if(ChanSize==4) ++FourSegs;

    }

     mf::LogVerbatim("GeometryTest") << "OneSeg: "       << OneSeg 
				     << ",  TwoSegs: "   << TwoSegs
				     << ",  ThreeSegs: " << ThreeSegs
				     << ",  FourSegs: "  << FourSegs;

  }

  //......................................................................
  void GeometryTest::printVolBounds()
  {
    art::ServiceHandle<geo::Geometry> geom;

      double origin[3] = {0.};
      double world[3] = {0.};
      for(unsigned int c = 0; c < geom->Ncryostats(); ++c){
	geom->Cryostat(c).LocalToWorld(origin, world);

        mf::LogVerbatim("GeometryTest") << "Cryo " << c;
	mf::LogVerbatim("GeometryTest") << "    -x: " << world[0] - geom->Cryostat(c).HalfWidth();
	mf::LogVerbatim("GeometryTest") << "    +x: " << world[0] + geom->Cryostat(c).HalfWidth();
	mf::LogVerbatim("GeometryTest") << "    -y: " << world[1] - geom->Cryostat(c).HalfHeight();
	mf::LogVerbatim("GeometryTest") << "    +y: " << world[1] + geom->Cryostat(c).HalfHeight();
	mf::LogVerbatim("GeometryTest") << "    -z: " << world[2] - geom->Cryostat(c).Length()/2;
	mf::LogVerbatim("GeometryTest") << "    +z: " << world[2] + geom->Cryostat(c).Length()/2;

        for(unsigned int t = 0; t < geom->NTPC(c); ++t){
          geom->Cryostat(c).TPC(t).LocalToWorld(origin, world);

          mf::LogVerbatim("GeometryTest") << "  TPC " << t;
          mf::LogVerbatim("GeometryTest") << "    -x: " << world[0] - geom->Cryostat(c).TPC(t).HalfWidth();
          mf::LogVerbatim("GeometryTest") << "    +x: " << world[0] + geom->Cryostat(c).TPC(t).HalfWidth();
          mf::LogVerbatim("GeometryTest") << "    -y: " << world[1] - geom->Cryostat(c).TPC(t).HalfHeight();
          mf::LogVerbatim("GeometryTest") << "    +y: " << world[1] + geom->Cryostat(c).TPC(t).HalfHeight();
          mf::LogVerbatim("GeometryTest") << "    -z: " << world[2] - geom->Cryostat(c).TPC(t).Length()/2;
          mf::LogVerbatim("GeometryTest") << "    +z: " << world[2] + geom->Cryostat(c).TPC(t).Length()/2;
        }
      }

  }



  //......................................................................
  // great sanity check for geometry, only call in analyze when debugging
  void GeometryTest::printDetDim()
  {
    art::ServiceHandle<geo::Geometry> geom;

    for(unsigned int c = 0; c < geom->Ncryostats(); ++c){

      mf::LogVerbatim("GeometryTest") << "Cryo " << c;
      mf::LogVerbatim("GeometryTest") << "    width: "
				      << geom->CryostatHalfWidth(c);
      mf::LogVerbatim("GeometryTest") << "    height: "
				      << geom->CryostatHalfHeight(c);
      mf::LogVerbatim("GeometryTest") << "    length: "
				      << geom->CryostatLength(c);

      mf::LogVerbatim("GeometryTest") << "  TPC 0";
      mf::LogVerbatim("GeometryTest") << "    width: "
				      << geom->DetHalfWidth(0,c);
      mf::LogVerbatim("GeometryTest") << "    height: "
				      << geom->DetHalfHeight(0,c);
      mf::LogVerbatim("GeometryTest") << "    length: "
				      << geom->DetLength(0,c);      
    }
  }

  //......................................................................
  // great sanity check for volume sorting, only call in analyze when debugging
  void GeometryTest::printWirePos()
  {
    art::ServiceHandle<geo::Geometry> geom;

    unsigned int cs = 0;

    for(unsigned int t=0; t<std::floor(geom->NTPC()/12)+1; ++t){
      for(unsigned int p=0; p<3; ++p){
        for(unsigned int w=0; w<geom->Cryostat(0).TPC(t).Plane(p).Nwires(); w++){
        
          double xyz[3] = {0.};
          geom->Cryostat(0).TPC(t).Plane(p).Wire(w).GetCenter(xyz);

          std::cout << "WireID (" << cs << ", " << t << ", " << p << ", " << w
                    << "):  x = " << xyz[0] 
                    << ", y = " << xyz[1]
                    << ", z = " << xyz[2] << std::endl;
        }
      }
    }
  }

  //......................................................................
  // great insanity: print all wires in a TPC
  void GeometryTest::printWiresInTPC
    (const geo::TPCGeo& tpc, std::string indent /* = "" */) const
  {
    const unsigned int nPlanes = tpc.Nplanes();
    const double Origin[3] = { 0., 0., 0. };
    double TPCpos[3];
    tpc.LocalToWorld(Origin, TPCpos);
    mf::LogVerbatim("GeometryTest") << indent << "TPC at ("
      << TPCpos[0] << ", " << TPCpos[1] << ", " << TPCpos[2]
      << ") cm has " << nPlanes << " wire planes:";
    for(unsigned int p = 0; p < nPlanes; ++p) {
      const geo::PlaneGeo& plane = tpc.Plane(p);
      const unsigned int nWires = plane.Nwires();
      double PlanePos[3];
      plane.LocalToWorld(Origin, PlanePos);
      std::string coord, orientation;
      switch (plane.View()) {
        case geo::kU:       coord = "U direction"; break;
        case geo::kV:       coord = "V direction"; break;
        case geo::kZ:       coord = "Z direction"; break;
        case geo::k3D:      coord = "3D coordinate"; break;
        case geo::kUnknown: coord = "an unknown direction"; break;
        default:            coord = "unexpected direction"; break;
      } // switch
      switch (plane.Orientation()) {
        case geo::kHorizontal: orientation = "horizontal"; break;
        case geo::kVertical:   orientation = "vertical"; break;
        default:               orientation = "unexpected"; break;
      }
      mf::LogVerbatim("GeometryTest") << indent << "  plane #" << p << " at ("
        << PlanePos[0] << ", " << PlanePos[1] << ", " << PlanePos[2] << ") cm"
        " has " << orientation << " orientation and "
        << nWires << " wires measuring " << coord << ":";
      for(unsigned int w = 0;  w < nWires; ++w) {
        const geo::WireGeo& wire = plane.Wire(w);
        double xyz[3] = { 0. };
        double WireS[3],  WireM[3], WireE[3]; // start, middle point and end
        wire.GetCenter(xyz);
        // the wire should be aligned on z axis, half on each side of 0,
        // in its local frame
        double WirePoint[3] = { 0., 0., 0. };
        wire.LocalToWorld(WirePoint, WireM);
        WirePoint[2] = wire.HalfL();
        wire.LocalToWorld(WirePoint, WireE);
        WirePoint[2] = -wire.HalfL();
        wire.LocalToWorld(WirePoint, WireS);
        mf::LogVerbatim("GeometryTest") << indent
          << "    wire #" << w
          << " at (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")"
          << "\n" << indent << "       start at (" << WireS[0] << ", " << WireS[1] << ", " << WireS[2] << ")"
          << "\n" << indent << "      middle at (" << WireM[0] << ", " << WireM[1] << ", " << WireM[2] << ")"
          << "\n" << indent << "         end at (" << WireE[0] << ", " << WireE[1] << ", " << WireE[2] << ")"
          ;
      } // for wire
    } // for plane
  } // GeometryTest::printWiresInTPC()

  
  void GeometryTest::printAllGeometry() const {
    art::ServiceHandle<geo::Geometry> geom;
    const unsigned int nCryostats = geom->Ncryostats();
    const double Origin[3] = { 0., 0., 0. };
    mf::LogVerbatim("GeometryTest") << "Detector " << geom->GetDetectorName()
      << " has " << nCryostats << " cryostats:";
    for(unsigned int c = 0; c < nCryostats; ++c) {
      const geo::CryostatGeo& cryostat = geom->Cryostat(c);
      const unsigned int nTPCs = cryostat.NTPC();
      double CryoPos[3];
      cryostat.LocalToWorld(Origin, CryoPos);
      mf::LogVerbatim("GeometryTest") << "  cryostat #" << c << " at ("
        << CryoPos[0] << ", " << CryoPos[1] << ", " << CryoPos[2] << ") cm has "
        << nTPCs << " TPC(s):";
      for(unsigned int t = 0;  t < nTPCs; ++t) {
        const geo::TPCGeo& tpc = cryostat.TPC(t);
        if (nTPCs > 1) mf::LogVerbatim("GeometryTest") << "    TPC #" << t;
        printWiresInTPC(tpc, "    ");
      } // for TPC
    } // for cryostat
    mf::LogVerbatim("GeometryTest") << "End of detector "
      << geom->GetDetectorName() << " geometry.";
  } // GeometryTest::printAllGeometry()

  //......................................................................
  void GeometryTest::testCryostat()
  {
    art::ServiceHandle<geo::Geometry> geom;

    mf::LogVerbatim("GeometryTest") << "\tThere are " << geom->Ncryostats() << " cryostats in the detector";

    for(unsigned int c = 0; c < geom->Ncryostats(); ++c){

      mf::LogVerbatim("GeometryTest") << "\n\t\tCryostat " << c 
				      << " " << geom->Cryostat(c).Volume()->GetName()
				      << " Dimensions: " << 2.*geom->Cryostat(c).HalfWidth()
				      << " x "           << 2.*geom->Cryostat(c).HalfHeight() 
				      << " x "           << geom->Cryostat(c).Length()
				      << "\n\t\t mass: "   << geom->Cryostat(c).Mass();

      double cryobound[6] = {0.};
      geom->CryostatBoundaries(cryobound, c);
      mf::LogVerbatim("GeometryTest") << "Cryostat boundaries are at:\n"
				      << "\t-x:" << cryobound[0] << " +x:" << cryobound[1]
				      << "\t-y:" << cryobound[2] << " +y:" << cryobound[3]
				      << "\t-z:" << cryobound[4] << " +z:" << cryobound[5];

      // pick a position in the middle of the cryostat in the world coordinates
      double worldLoc[3] = {0.5*(cryobound[1] - cryobound[0]) + cryobound[0],
			    0.5*(cryobound[3] - cryobound[2]) + cryobound[2],
			    0.5*(cryobound[5] - cryobound[4]) + cryobound[4]};
		
      LOG_DEBUG("GeometryTest") << "\t testing Geometry::PoitionToCryostat....";
      try{
	unsigned int cstat = 0;
	geom->PositionToCryostat(worldLoc, cstat);
      }
      catch(cet::exception &e){
	mf::LogWarning("FailedToLocateCryostat") << "\n exception caught:" << e;
      }
      LOG_DEBUG("GeometryTest") << "done";

      LOG_DEBUG("GeometryTest") << "\t Now test the TPCs associated with this cryostat";
      this->testTPC(c);
    }

    return;
  }

  //......................................................................
  void GeometryTest::testTPC(unsigned int const& c)
  {
    art::ServiceHandle<geo::Geometry> geom;

    mf::LogVerbatim("GeometryTest") << "\tThere are " << geom->Cryostat(c).NTPC() 
				    << " TPCs in the detector";
    
    for(size_t t = 0; t < geom->Cryostat(c).NTPC(); ++t){
      mf::LogVerbatim("GeometryTest") << std::endl << "\t\tTPC " << t 
				      << " " << geom->GetLArTPCVolumeName(t, c) 
				      << " has " 
				      << geom->Cryostat(c).TPC(t).Nplanes() << " planes.";
      for(size_t p = 0; p < geom->Cryostat(c).TPC(t).Nplanes(); ++p)
	mf::LogVerbatim("GeometryTest") << std::endl << "\t\tPlane " << p << " has " 
					<< geom->Cryostat(c).TPC(t).Plane(p).Nwires() 
					<< " wires and is at (x,y,z) = (" 
					<< geom->Cryostat(c).TPC(t).PlaneLocation(p)[0] << "," 
					<< geom->Cryostat(c).TPC(t).PlaneLocation(p)[1] << "," 
					<< geom->Cryostat(c).TPC(t).PlaneLocation(p)[2] 
					<< "); \n\t\tpitch from plane 0 is "
					<< geom->Cryostat(c).TPC(t).Plane0Pitch(p) << "; \n\t\tOrientation "
					<< geom->Cryostat(c).TPC(t).Plane(p).Orientation() << ", View "
					<< geom->Cryostat(c).TPC(t).Plane(p).View() << ", Wire angle "
					<< geom->Cryostat(c).TPC(t).Plane(p).Wire(0).ThetaZ()
					<< "\n\t\tTPC Dimensions: " << 2.*geom->TPC(t).HalfWidth()
					<< " x " << 2.*geom->Cryostat(c).TPC(t).HalfHeight() 
					<< " x " << geom->Cryostat(c).TPC(t).Length()
					<< "\n\t\tTPC Active Dimensions: " 
					<< 2.*geom->Cryostat(c).TPC(t).ActiveHalfWidth()
					<< " x " << 2.*geom->Cryostat(c).TPC(t).ActiveHalfHeight() 
					<< " x " << geom->Cryostat(c).TPC(t).ActiveLength()
					<< "\n\t\tTPC mass: " << geom->Cryostat(c).TPC(t).ActiveMass()
					<< "\n\t\tTPC drift distance: " 
					<< geom->Cryostat(c).TPC(t).DriftDistance();

      geo::DriftDirection_t dir = geom->Cryostat(c).TPC(t).DriftDirection();
      if     (dir == geo::kNegX) 
	mf::LogVerbatim("GeometryTest") << "\t\tdrift direction is towards negative x values";
      else if(dir == geo::kPosX) 
	mf::LogVerbatim("GeometryTest") << "\t\tdrift direction is towards positive x values";
      else{
	throw cet::exception("UnknownDriftDirection") << "\t\tdrift direction is unknown\n";
      }

      LOG_DEBUG("GeometryTest") << "\t testing PositionToTPC...";
      // pick a position in the middle of the cryostat in the world coordinates
      double worldLoc[3] = {0.};
      double localLoc[3] = {0.};
      geom->Cryostat(c).TPC(t).LocalToWorld(localLoc, worldLoc);

      unsigned int tpc   = UINT_MAX;
      geom->Cryostat(c).PositionToTPC(worldLoc,tpc,1+1.e-4);

      if(tpc != t)
	throw cet::exception("BadTPCLookupFromPosition") << "TPC look up returned tpc = "
							 << tpc << " should be " << t << "\n";

      LOG_DEBUG("GeometryTest") << "done.";
    }
    
    return;
  }


  //......................................................................
  void GeometryTest::testChannelToWire()
  {
    art::ServiceHandle<geo::Geometry> geom;


    for(unsigned int cs = 0; cs < geom->Ncryostats(); ++cs){
      for(unsigned int tpc = 0; tpc < geom->Cryostat(cs).NTPC(); ++tpc){
	for(unsigned int plane = 0; plane < geom->Cryostat(cs).TPC(tpc).Nplanes(); ++plane){
	  for(unsigned int wire = 0; wire < geom->Cryostat(cs).TPC(tpc).Plane(plane).Nwires(); ++wire){

	    uint32_t channel = geom->PlaneWireToChannel(plane, wire, tpc, cs);
	    //std::cout << "WireID (" << cs << ", " << tpc << ", " << plane << ", " << wire 
	    //	<< ") --> Channel " << channel << std::endl;    
	    std::vector< geo::WireID > wireIDs = geom->ChannelToWire(channel);
	    

	    if ( wireIDs.size() == 0 ) 
	      throw cet::exception("BadChannelLookup") << "requested channel: " << channel 
						       << ";" << cs << "," << tpc
						       << "," << plane << "," << wire << "\n"
						       << "got back an empty vector of WireID " << "\n";

	    bool goodLookup = false;
	    for( auto const& wid : wireIDs){
	      if(wid.Cryostat == cs    && 
		 wid.TPC      == tpc   && 
		 wid.Plane    == plane && 
		 wid.Wire     == wire) goodLookup = true;
	    }
	    
	    if(!goodLookup)
	    {
	      std::cout << "Returned: " << std::endl;
              for(unsigned int id=0; id<wireIDs.size(); ++id)
	      {
		std::cout << "wireIDs[" << id << "] = ("
		          << wireIDs[id].Cryostat << ", "
		          << wireIDs[id].TPC      << ", "
		          << wireIDs[id].Plane    << ", "
		          << wireIDs[id].Wire     << ")" << std::endl;
              }
	      throw cet::exception("BadChannelLookup") << "requested channel " << channel 
						       << "expected to return" << cs << "," << tpc
						       << "," << plane << "," << wire << "\n"
						       << "no returned geo::WireID structs matched\n";
            }

	    if(geom->SignalType(channel) != geom->Plane(plane, tpc, cs).SignalType() )
	      throw cet::exception("BadChannelLookup") << "expected signal type: SignalType(channel) = " 
						       << geom->SignalType(channel)
						       << " for channel " 
						       << channel << ", WireID ("  
						       << cs << ", " << tpc << ", " << plane << ", " << wire
						       << "), got: Plane(" << plane << ", " << tpc 
						                           << ", " << cs << ").SignalType() = "
						       << geom->Plane(plane, tpc, cs).SignalType() << "\n";


	    if(geom->View(channel) != geom->Plane(plane, tpc, cs).View() )
	      throw cet::exception("BadChannelLookup") << "expected view type: View(channel) = " 
						       << geom->View(channel)
						       << " for channel " 
						       << channel << ", WireID ("  
						       << cs << ", " << tpc << ", " << plane << ", " << wire
						       << "), got: Plane(" << plane << ", " << tpc 
						                           << ", " << cs << ").View() = "
						       << geom->Plane(plane, tpc, cs).View() << "\n";

	  }
	}
      }
    }

    return;
  }

  //......................................................................
  void GeometryTest::testFindPlaneCenters()
  {
    art::ServiceHandle<geo::Geometry> geom;

    double xyz[3] = {0.},   xyzW[3] = {0.};
    for(size_t i = 0; i < geom->Nplanes(); ++i){ 
      geom->Plane(i).LocalToWorld(xyz,xyzW);
      mf::LogVerbatim("GeometryTest") << "\n\tplane " 
				      << i << " is centered at (x,y,z) = (" 
				      << xyzW[0] << "," << xyzW[1]
				      << "," << xyzW[2] << ")";
    } 
  } 

  //......................................................................
  void GeometryTest::testStandardWirePos() 
  {
    art::ServiceHandle<geo::Geometry> geom;

    double xyz[3] = {0.};
    double xyzprev[3] = {0.};
    for(size_t cs = 0; cs < geom->Ncryostats(); ++cs){
      for(size_t t = 0; t < geom->Cryostat(cs).NTPC(); ++t){
	const geo::TPCGeo* tpc = &geom->Cryostat(cs).TPC(t); 

	for (size_t i=0; i < tpc->Nplanes(); ++i) {
	  const geo::PlaneGeo* plane = &tpc->Plane(i);

	  for (size_t j = 1; j < plane->Nwires(); ++j) {

	    const geo::WireGeo wire = plane->Wire(j);
	    const geo::WireGeo wireprev = plane->Wire(j-1);

	    wire.GetCenter(xyz);
	    wireprev.GetCenter(xyzprev);

	    // wires increase in +z order
	    if(xyz[2] < xyzprev[2])
	      throw cet::exception("WireOrderProblem") 	<< "\n\twires do not increase in +z order in"
							<< "Cryostat " << cs
							<< ", TPC " << t
							<< ", Plane " << i
							<< ";  at wire " << j << "\n";

	  }// end loop over wires
	}// end loop over planes
      }// end loop over tpcs
    }// end loop over cryostats

}

  //......................................................................
  void GeometryTest::testAPAWirePos() 
  {
    art::ServiceHandle<geo::Geometry> geom;

    double origin[3] = {0.};
    double tpcworld[3] = {0.};
    double xyz[3] = {0.};
    double xyzprev[3] = {0.};
    for(size_t cs = 0; cs < geom->Ncryostats(); ++cs){
      for(size_t t = 0; t < geom->Cryostat(cs).NTPC(); ++t){
	const geo::TPCGeo* tpc = &geom->Cryostat(cs).TPC(t);
	tpc->LocalToWorld(origin, tpcworld);

	for (size_t i=0; i < tpc->Nplanes(); ++i) {
	  const geo::PlaneGeo* plane = &tpc->Plane(i);

	  for (size_t j = 1; j < plane->Nwires(); ++j) {
	    const geo::WireGeo wire = plane->Wire(j);
	    const geo::WireGeo wireprev = plane->Wire(j-1);

	    wire.GetCenter(xyz);
	    wireprev.GetCenter(xyzprev);

            // top TPC wires increase in -y
	    if(tpcworld[1] > 0 && xyz[1] > xyzprev[1])
	      throw cet::exception("WireOrderProblem") 	<< "\n\ttop TPC wires do not increase in -y order in"
							<< "Cryostat " << cs
							<< ", TPC " << t
							<< ", Plane " << i
							<< ";  at wire " << j << "\n";
            // bottom TPC wires increase in +y
	    else if(tpcworld[1] < 0 && xyz[1] < xyzprev[1])
	      throw cet::exception("WireOrderProblem") 	<< "\n\tbottom TPC wires do not increase in +y order in"
                                                        << "Cryostat " << cs
							<< ", TPC " << t
                                                        << ", Plane " << i 
                                                        << ";  at wire " << j << "\n";
	  }// end loop over wires
	}// end loop over planes
      }// end loop over tpcs
    }// end loop over cryostats

}


  //......................................................................
  void GeometryTest::testNearestWire()
  {
    art::ServiceHandle<geo::Geometry> geom;
                                                       
    // Even if you comment it out, please leave the TStopWatch code
    // in this code for additional testing. The NearestChannel routine
    // is the most frequently called in the simulation, so its execution time
    // is an important component of LArSoft's speed.
    TStopwatch stopWatch;
    stopWatch.Start();

    // get a wire and find its center
    for(unsigned int cs = 0; cs < geom->Ncryostats(); ++cs){
      for(unsigned int t = 0; t < geom->Cryostat(cs).NTPC(); ++t){
	for(unsigned int p = 0; p < geom->Cryostat(cs).TPC(t).Nplanes(); ++p){
	  for(unsigned int w = 0; w < geom->Cryostat(cs).TPC(t).Plane(p).Nwires(); ++w){
	
	    const geo::WireGeo& wire = geom->Cryostat(cs).TPC(t).Plane(p).Wire(w);
	    const double pos[3] = {0., 0.0, 0.};
	    double posWorld[3] = {0.};
	    wire.LocalToWorld(pos, posWorld);

	    uint32_t nearest = 0;
	    std::vector< geo::WireID > wireIDs;

	    try{
	      // The double[] version tested here falls back on the
	      // TVector3 version, so this test both.
	      nearest = geom->NearestChannel(posWorld, p, t, cs);

	      // We also want to test the std::vector<duoble> version
	      std::vector<double> posWorldV(3);
	      for (int i=0; i<3; ++i) {
		posWorldV[i] = posWorld[i] + 0.001;
	      }
	      nearest = geom->NearestChannel(posWorldV, p, t, cs);
	    }
	    catch(cet::exception &e){
	      mf::LogWarning("GeoTestCaughtException") << e;
	    }

	    try{
	      wireIDs = geom->ChannelToWire(nearest);

	      if ( wireIDs.size() == 0 ) {
		throw cet::exception("BadPositionToChannel") << "test point is at " 
							     << posWorld[0] << " " 
							     << posWorld[1] << " " 
							     << posWorld[2] << "\n"
							     << "nearest channel is " 
							     << nearest << " for " 
							     << cs << " " << t << " "
							     << p << " " << w << "\n";
	      }
	    }
	    catch(cet::exception &e){
	      mf::LogWarning("GeoTestCaughtException") << e;
	    }

            bool goodLookup = false;
            for( auto const& wid : wireIDs){
              if(wid.Cryostat == cs    &&
                 wid.TPC      == t     &&
                 wid.Plane    == p     &&
                 wid.Wire     == w   ) goodLookup = true;
            }

	    if(!goodLookup){
	      throw cet::exception("BadPositionToChannel") << "Current WireID ("
							   << cs << "," << t << "," << p << "," << w << ") "
							   << "has a world position at "
							   << posWorld[0] << " " 
							   << posWorld[1] << " " 
							   << posWorld[2] << "\n"
							   << "NearestWire for this position is "
							   << geom->NearestWire(posWorld,p,t,cs) << "\n"
							   << "NearestChannel is " 
							   << nearest << " for " 
							   << cs << " " << t << " " << p << " " << w << "\n"
							   << "Should be channel "
							   << geom->PlaneWireToChannel(p,w,t,cs);
	    } // if good lookup fails
	  } // end loop over wires
	} // end loop over planes
      }// end loop over tpcs
    }// end loop over cryostats

    stopWatch.Stop();
    LOG_DEBUG("GeometryTest") << "\tdone testing closest channel";
    stopWatch.Print();
    
    // trigger an exception with NearestChannel
    mf::LogVerbatim("GeometryTest") << "\tattempt to cause an exception to be caught "
				    << "when looking for a nearest channel";

    double posWorld[3] = {geom->CryostatHalfWidth()*2,
			  geom->CryostatHalfHeight()*2,
			  geom->CryostatLength()*2};

    try{
      geom->NearestChannel(posWorld, 0, 0, 0);
    }
    catch(cet::exception &e){
      mf::LogWarning("GeoTestCaughtException") << e;
    }

  }

  //......................................................................
  void GeometryTest::testWirePitch()
  {
    art::ServiceHandle<geo::Geometry> geom;

    // loop over all planes and wires to be sure the pitch is consistent

    // hard code the value we think it should be for each detector
    double shouldbe[3];
    if(geom->DetId() == geo::kArgoNeuT){
      shouldbe[0] = 0.4;
      shouldbe[1] = 0.4;
      shouldbe[2] = 0.4; 
    }
    else if(geom->DetId() == geo::kMicroBooNE
	    || geom->DetId() == geo::kICARUS){	
      shouldbe[0] = 0.3;
      shouldbe[1] = 0.3;
      shouldbe[2] = 0.3; 
    }
    else if(geom->DetId() == geo::kLBNE35t 
	    || geom->DetId() == geo::kLBNE10kt
	    || geom->DetId() == geo::kLBNE34kt){
      shouldbe[0] = 0.49;
      shouldbe[1] = 0.5;
      shouldbe[2] = 0.45;  
    }
    else if(geom->DetId() == geo::kBo){	
      shouldbe[0] = 0.46977;
      shouldbe[1] = 0.46977;
      shouldbe[2] = 0.46977; 
    }

    for(size_t c = 0; c < geom->Ncryostats(); ++c){
      for(size_t t = 0; t < geom->Cryostat(c).NTPC(); ++t){
        for(size_t p = 0; p < geom->Cryostat(c).TPC(t).Nplanes(); ++p){
	  for(size_t w = 0; w < geom->Cryostat(c).TPC(t).Plane(p).Nwires()-1; ++w){
	    // get the wire pitch
	    double pitch = geom->Cryostat(c).TPC(t).WirePitch(w, w+1, p);
	    if(std::abs(pitch - shouldbe[p]) > 0.01*shouldbe[p]){
	      throw cet::exception("UnexpectedWirePitch") << "\n\tpitch is " 
	 						  << pitch << " instead of " << shouldbe[p] 
							  << " for cryostat " << c
							  << ", tpc " << t
							  << ", plane " << p
							  << "; wires: " << w << ", " << w+1 << "\n";
	    }// end if pitch is wrong
	  }// end loop over wires
        }// end loop over planes
      }// end loop over TPCs
    }// end loop over cryostats

  }

  //......................................................................
  void GeometryTest::testPlanePitch()
  {
    art::ServiceHandle<geo::Geometry> geom;

    // loop over all planes to be sure the pitch is consistent

    // hard code the value we think it should be for each detector
    double shouldbe = 0.4; // true for ArgoNeuT
    if(geom->DetId() == geo::kMicroBooNE)         shouldbe = 0.3;
    else if(geom->DetId() == geo::kLBNE35t 
	    || geom->DetId() == geo::kLBNE10kt
	    || geom->DetId() == geo::kLBNE34kt)  shouldbe = 0.5;
    else if(geom->DetId() == geo::kBo)           shouldbe = 0.65;
    else if(geom->DetId() == geo::kICARUS)       shouldbe = 0.476;

    for(size_t t = 0; t < geom->NTPC(); ++t){
      for(size_t p = 0; p < geom->TPC(t).Nplanes()-1; ++p){
	double pitch = std::abs(geom->TPC(t).PlanePitch(p, p+1));
	if(std::abs(pitch - shouldbe) > 0.1*shouldbe){
	  throw cet::exception("UnexpectedPlanePitch") << "\n\tunexpected pitch: " 
						       << pitch << "/" << shouldbe << "\n"; 
	}// end if wrong pitch
      }// end loop over planes
    }// end loop over TPCs

  }

  //......................................................................

  void GeometryTest::testStepping()
  {
    art::ServiceHandle<geo::Geometry> geom;

    //
    // Test stepping. Example is similar to what one would do for photon
    // transport. Rattles photons around inside the scintillator
    // bouncing them off walls.
    //
    double xyz[3]      = {0.};
    double xyzWire[3]  = {0.};
    double dxyz[3]     = {0.};
    double dxyzWire[3] = {0, sin(0.1), cos(0.1)};

    geom->Plane(1).Wire(0).LocalToWorld(xyzWire,xyz);
    geom->Plane(1).Wire(0).LocalToWorldVect(dxyzWire,dxyz);

    mf::LogVerbatim("GeometryTest") << "\n\t" << xyz[0]  << "\t" << xyz[1]  << "\t" << xyz[2] ;
    mf::LogVerbatim("GeometryTest") << "\t"   << dxyz[0] << "\t" << dxyz[1] << "\t" << dxyz[2];

    gGeoManager->InitTrack(xyz, dxyz);
    for (int i=0; i<10; ++i) {
      const double* pos = gGeoManager->GetCurrentPoint();
      const double* dir = gGeoManager->GetCurrentDirection();
      mf::LogVerbatim("GeometryTest") << "\tnode = " 
				      << gGeoManager->GetCurrentNode()->GetName()
				      << "\n\t\tpos=" << "\t"
				      << pos[0] << "\t"
				      << pos[1] << "\t"
				      << pos[2]
				      << "\n\t\tdir=" << "\t"
				      << dir[0] << "\t"
				      << dir[1] << "\t"
				      << dir[2]
				      << "\n\t\tmat = " 
				      << gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->GetName();
      
      gGeoManager->FindNextBoundary();
      gGeoManager->FindNormal();
      gGeoManager->Step(kTRUE,kTRUE);
    }

    xyz[0] = 306.108; xyz[1] = -7.23775; xyz[2] = 856.757;
    gGeoManager->InitTrack(xyz, dxyz);
    mf::LogVerbatim("GeometryTest") << "\tnode = " 
				    << gGeoManager->GetCurrentNode()->GetName()
				    << "\n\tmat = " 
				    << gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->GetName();

    gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->Print();

  }

  //......................................................................

  void GeometryTest::testProject() 
  {
    art::ServiceHandle<geo::Geometry> geom;

    double xlo, xhi;
    double ylo, yhi;
    double zlo, zhi;
    geom->WorldBox(&xlo, &xhi, &ylo, &yhi, &zlo, &zhi);
  
    double xyz[3]   = { 0.0, 0.0, 0.0};
    double dxyz1[3] = { 1.0, 0.0, 0.0};
    double dxyz2[3] = {-1.0, 0.0, 0.0};
    double dxyz3[3] = { 0.0, 1.0, 0.0};
    double dxyz4[3] = { 0.0,-1.0, 0.0};
    double dxyz5[3] = { 0.0, 0.0, 1.0};
    double dxyz6[3] = { 0.0, 0.0,-1.0};

    double xyzo[3];
    geo::ProjectToBoxEdge(xyz, dxyz1, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[0]-xhi)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz2, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[0]-xlo)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz3, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[1]-yhi)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz4, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[1]-ylo)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz5, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[2]-zhi)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz6, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[2]-zlo)>1.E-6) abort();
  }


}//end namespace


namespace geo{

  DEFINE_ART_MODULE(GeometryTest)

}

#endif
