/**
 * @file CRTGen_module.cc
 * @brief Producer generating Monte Carlo truth record in LArSoft format to simulate protodunedp CRT Trigger
 * @author Jose Soto
 */
/**
 * @class evgen::CRTGen
 *  This module assumes muons cross uniformly both CRT pannels. One muon will be generated per event.
 */
#include <string>
#include <fstream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "TLorentzVector.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "CLHEP/Random/RandFlat.h"

namespace evgen {
  class CRTGen;
}

class evgen::CRTGen : public art::EDProducer {
public:
  explicit CRTGen(fhicl::ParameterSet const & p);

  void produce(art::Event & e)                    override;
  void beginJob()               		  override;
  void beginRun(art::Run & run) 		  override;

private:
  short int driftcoordinate=0;
  std::vector<double> CRT_TOP_centerY={-581,305.7,135};//drift in Y geometry in cm
  std::vector<double> CRT_BOT_centerY={+581,-202.7,135}; //drift in Y geometry in cm
  double CRTHeight=8*14.0; //size in Y
  double CRTLengthZ=144; //size in Z
//  std::ifstream* fInputFile;
//  std::string    fInputFileName; ///< Name of text file containing events to simulate
};

//------------------------------------------------------------------------------
evgen::CRTGen::CRTGen(fhicl::ParameterSet const & p)
  : EDProducer{p}
//  , fInputFile(0)
//  , fInputFileName{p.get<std::string>("InputFileName")}
//  , fMoveY{p.get<double>("MoveY", -1e9)}

{
//  if (fMoveY>-1e8){
//    mf::LogWarning("CRTGen")<<"Particles will be moved to a new plane y = "<<fMoveY<<" cm.\n";
//  }

  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();
}

//------------------------------------------------------------------------------
void evgen::CRTGen::beginJob()
{
  //fInputFile = new std::ifstream(fInputFileName.c_str());

  // check that the file is a good one
  //if( !fInputFile->good() )
  //  throw cet::exception("CRTGen") << "input text file "
  //					<< fInputFileName
//					<< " cannot be read.\n";
}

//------------------------------------------------------------------------------
void evgen::CRTGen::beginRun(art::Run& run)
{
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
    driftcoordinate = (geo->TPC(0)).DetectDriftDirection();

    if( driftcoordinate==1 || driftcoordinate==2 )
      {
	std::cout<<" drift coordinate: "<<driftcoordinate<<std::endl;
      }else{ throw cet::exception("CRTGen") << "unknown drift coordinate "                                                                            << driftcoordinate << " \n"; }  

    /*    +1: positive x
      +2: positive y
      +3: positive z
      -1: negative x
      -2: negative y
      -3: negative z
      0: other (or algorithm failed)
    */
 }

//------------------------------------------------------------------------------
void evgen::CRTGen::produce(art::Event & e)
{
// check that the file is still good
//  if( !fInputFile->good() )
//    throw cet::exception("CRTGen") << "input text file "
//					<< fInputFileName
//					<< " cannot be read in produce().\n";

  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
  simb::MCTruth truth;

  //We simulate crossing muons from CRT TOP to CRT BOT

  std::vector<double> CRT_TOP_center=CRT_TOP_centerY; //{-581,305.7,135};//drift in Y geometry in cm
  std::vector<double> CRT_BOT_center=CRT_BOT_centerY;//{+581,-202.7,135}; //drift in Y geometry in cm
  double CRTHeightY=CRTHeight; //8*14.0; //size in Y ..TODO..
  double CRTHeightX=0.0;
  double CRTLengthZ=144; //size in Z
  //drift in x-coordinate//
  if(driftcoordinate==1) // Y->X, X->-Y
    {
      CRT_TOP_center[0] = CRT_TOP_centerY[1];
      CRT_TOP_center[1] = CRT_TOP_centerY[0];
      CRT_BOT_center[0] = CRT_BOT_centerY[1];
      CRT_BOT_center[1] = CRT_BOT_centerY[0];
      CRTHeightX=CRTHeight;
      CRTHeightY=0;
    }

  // declare the variables for reading in the event record
//  int            status         = 0;
  int 	 	 pdg            = 13;
//  int 	 	 firstMother    = 0;
  double 	 energy      	= CLHEP::RandFlat::shoot(10,10); //uniform distribution among 2-3GeV
  double 	 mass        	= 0.1056583745;//GeV
  double 	 xPosition   	= CLHEP::RandFlat::shoot(CRT_TOP_center[0]-0.5*CRTHeightX,CRT_TOP_center[0]+0.5*CRTHeightX);
  double 	 yPosition   	= CLHEP::RandFlat::shoot(CRT_TOP_center[1]-0.5*CRTHeightY,CRT_TOP_center[1]+0.5*CRTHeightY);
  double 	 zPosition   	= CLHEP::RandFlat::shoot(CRT_TOP_center[2]-0.5*CRTLengthZ,CRT_TOP_center[2]+0.5*CRTLengthZ);
  double 	 time        	= 0.;
  double 	 xPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[0]-0.5*CRTHeightX,CRT_BOT_center[0]+0.5*CRTHeightX);
  double 	 yPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[1]-0.5*CRTHeightY,CRT_BOT_center[1]+0.5*CRTHeightY);
  double 	 zPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[2]-0.5*CRTLengthZ,CRT_BOT_center[2]+0.5*CRTLengthZ);

  double totmom = sqrt(pow(energy,2)-pow(mass,2));


  double dx=xPositionEnd-xPosition;
  double dy=yPositionEnd-yPosition;
  double dz=zPositionEnd-zPosition;
  double norm=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
  dx/=norm;
  dy/=norm;
  dz/=norm;

  double 	 xMomentum      = dx*totmom;
  double 	 yMomentum   	= dy*totmom;
  double 	 zMomentum   	= dz*totmom;
  std::cout << "Shooting muon on " << xPosition << " " << yPosition << " " << zPosition <<  "to "<< xPositionEnd << " " << yPositionEnd << " " << zPositionEnd <<" With momentum: " << xMomentum << " " << yMomentum << " " << zMomentum << " E=" << energy << " m=" << mass << std::endl;

  TLorentzVector pos(xPosition, yPosition, zPosition, time);
  TLorentzVector mom(xMomentum, yMomentum, zMomentum, energy);

  simb::MCParticle part(-1, pdg, "primary");//, firstMother, mass, status);
  part.AddTrajectoryPoint(pos, mom);

  truth.Add(part);

  truthcol->push_back(truth);

  e.put(std::move(truthcol));
}

DEFINE_ART_MODULE(evgen::CRTGen)
