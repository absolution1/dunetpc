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
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TLorentzVector.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "CLHEP/Random/RandFlat.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

#include"TH1D.h"
#include"TH2D.h"
#include"TFile.h"
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
  std::vector<double> CRT_TOP_center_driftY={-581,305.7,135};//drift in Y geometry in cm
  std::vector<double> CRT_BOT_center_driftY={+581,-202.7,135}; //drift in Y geometry in cm
  double CRTHeight=8*14.0; //size in the vertical direction
  double CRTLength=144; //size in Z - Along the scintillator bars.

  std::vector<double> CRT_TOP_center; //corrected to the actual geometry
  std::vector<double> CRT_BOT_center; //corrected to the actual geometry
  double CRTSizeX; //size of the CRT in x,y and z corrected to the actual geometry
  double CRTSizeY;
  double CRTSizeZ;

  TH2D CRTTop, CRTBot;
  TH1D EnergyDistribution;

  int fmode;
  int fEnergyDistributionMode;
  std::pair<float,float> fEnergyRange;
  std::string    fInputFileNameCRT; ///< Name of text file containing events to simulate
  std::string    fInputFileNameEnergy; ///< Name of text file containing events to simulate


  TH2F *fTH2CRTTop;
  TH2F *fTH2CRTBot;
  TH1F *fTH1Energy;
  /* We might have some muons that are pointing out of the CRT bottom,
  but they end up reaching CRT bottom due to the scattering, to take those into account, we include a buffer */
  double BufferLengthOnCRTBottom; //buffer in cm to use when using a uniform distribution
};

//------------------------------------------------------------------------------
evgen::CRTGen::CRTGen(fhicl::ParameterSet const & p)
  : EDProducer{p}
  , fmode{p.get<int>("Mode",0)} // 0 for uniform distribution on CRT geoometry, 1 to get the distribution from TH2D
  , fEnergyDistributionMode{p.get<int>("EnergyDistribution",0)} // 0 for uniform distribution on CRT geoometry, 1 to get the distribution from TH2D
  , fEnergyRange{p.get<std::pair<float,float>>("EnergyRange",std::make_pair(2,3))}
  , fInputFileNameCRT{p.get<std::string>("InputFileNameCRT","CRT_RawInputs.root")}
  , fInputFileNameEnergy{p.get<std::string>("InputFileNameEnergy","MuonEnergy.root")}
  , BufferLengthOnCRTBottom{p.get<float>("BufferLengthOnCRTBottom",30.0)}
{
  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();
}

//------------------------------------------------------------------------------
void evgen::CRTGen::beginJob()
{
  if( fmode==0 )
  {
    std::cout<<" fmode=0 - Sending muons uniformly distributed on the CRT panels." << std::endl;
  }else if (fmode==1)
  {
    std::cout<<" fmode=1 - Sending muons following the distribution contained on file: " << fInputFileNameCRT << std::endl;}
  else{ throw cet::exception("CRTGen") << "unknown fmode value " << fmode << " \n"; }  

  if( fEnergyDistributionMode==0 )
  {
    std::cout<<" fEnergyDistributionMode=0 - Muons will follow a uniform distribution between the values: " << fEnergyRange.first<< " " << fEnergyRange.second << std::endl;
  }else if (fEnergyDistributionMode==1)
  {
    std::cout<<" fEnergyDistributionMode=1 - Sending muons following the distribution contained on file: " << fInputFileNameEnergy << std::endl;}
  else{ throw cet::exception("CRTGen") << "unknown fEnergyDistributionMode value " << fEnergyDistributionMode << " \n"; }  


  if(fmode==1)
  {
    TFile *fInputFileCRT = new TFile(fInputFileNameCRT.c_str(),"READ");
    // check that the file is a good one
    if( !fInputFileCRT->IsOpen() )
      throw cet::exception("CRTGen") << "input text file "
  					<< fInputFileNameCRT
					<< " cannot be read.\n";

    TH2D *h =(TH2D*)fInputFileCRT->Get("CRTTop");
    if (!h) throw cet::exception("CRTGen") << "TH2D named CRTTop not found in "
  					<< fInputFileNameCRT
					<< ".\n";
    CRTTop = *h;
    h =(TH2D*)fInputFileCRT->Get("CRTBot");
    if (!h) throw cet::exception("CRTGen") << "TH2D named CRTBot not found in "
  					<< fInputFileNameCRT
					<< ".\n";
    CRTBot = *h;
    fInputFileCRT->Close();

  }

  if(fEnergyDistributionMode==1)
  {
    TFile *fInputFileEnergy = new TFile(fInputFileNameEnergy.c_str(),"READ");
    // check that the file is a good one
    if( !fInputFileEnergy->IsOpen() )
      throw cet::exception("CRTGen") << "input root file "
  					<< fInputFileNameEnergy
					<< " cannot be read.\n";

    TH1D *h =(TH1D*)fInputFileEnergy->Get("EnergyDistribution");
    if (!h) throw cet::exception("CRTGen") << "TH1D named EnergyDistribution not found in "
  					<< fInputFileNameEnergy
					<< ".\n";
    EnergyDistribution = *h;
    fInputFileEnergy->Close();
  }

  art::ServiceHandle<art::TFileService> tfs;

  fTH2CRTTop = tfs->make<TH2F>("CRTTop","Start position at CRT Top",10,0,0,10,0,0);
  fTH2CRTBot = tfs->make<TH2F>("CRTBot","Expected end track position at CRT Bot",10,0,0,10,0,0);
  fTH1Energy = tfs->make<TH1F>("TH1Energy","Muon energy GeV",100,0,0);

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
      }else{ throw cet::exception("CRTGen") << "unknown drift coordinate " << driftcoordinate << " \n"; }  


    /*    +1: positive x
      +2: positive y
      +3: positive z
      -1: negative x
      -2: negative y
      -3: negative z
      0: other (or algorithm failed)
    */


  // We fix the CRT coordinates assuming drift in Y geometry
  CRT_TOP_center=CRT_TOP_center_driftY; //{-581,305.7,135};//drift in Y geometry in cm
  CRT_BOT_center=CRT_BOT_center_driftY;//{+581,-202.7,135}; //drift in Y geometry in cm
  CRTSizeY=CRTHeight; 
  CRTSizeX=0.0;
  CRTSizeZ=CRTLength; //size in Z for drift in Y and X

  //We correct for the drift in X geometry
  if(driftcoordinate==1) // Y->X, X->-Y
  {
     // CRT_TOP_center {305.7,581,135};//drift in X geometry in cm
     // CRT_BOT_center {-202.7,-581,135};//drift in X geometry in cm

    CRT_TOP_center[0] = CRT_BOT_center_driftY[1];
    CRT_TOP_center[1] = CRT_BOT_center_driftY[0];
    CRT_BOT_center[0] = CRT_BOT_center_driftY[1];
    CRT_BOT_center[1] = -CRT_BOT_center_driftY[0];
    CRTSizeX=CRTHeight;
    CRTSizeY=0;
    CRTSizeZ=CRTLength;
  }

 }

//------------------------------------------------------------------------------
void evgen::CRTGen::produce(art::Event & e)
{

  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
  simb::MCTruth truth;

  // declare the variables for reading in the event record
  int 	 	 pdg            = 13;
  double 	 energy;

  if(fEnergyDistributionMode==1) energy = EnergyDistribution.GetRandom();
  else energy = CLHEP::RandFlat::shoot(fEnergyRange.first,fEnergyRange.second); //uniform distribution among values set in the range

  double 	 mass        	= 0.1056583745;//muon mass in GeV
  double xPosition, yPosition, zPosition, xPositionEnd, yPositionEnd, zPositionEnd;
  if(fmode==1)
  {
    if(driftcoordinate==1)
    { //drift in X
      CRTTop.GetRandom2(zPosition,xPosition);
      CRTBot.GetRandom2(zPositionEnd,xPositionEnd);
      yPosition   	= CRT_TOP_center[1];
      yPositionEnd  	= CRT_BOT_center[1];

      fTH2CRTTop->Fill(zPosition,xPosition);
      fTH2CRTBot->Fill(zPositionEnd,xPositionEnd);
    }
    if(driftcoordinate==2)
    { //drift in Y
      CRTTop.GetRandom2(zPosition,yPosition);
      CRTBot.GetRandom2(zPositionEnd,yPositionEnd);
      xPosition   	= CRT_TOP_center[0];
      xPositionEnd  	= CRT_BOT_center[0];
      fTH2CRTTop->Fill(zPosition,yPosition);
      fTH2CRTBot->Fill(zPositionEnd,yPositionEnd);

    }
  }
  else
  {
    xPosition   	= CLHEP::RandFlat::shoot(CRT_TOP_center[0]-0.5*CRTSizeX,CRT_TOP_center[0]+0.5*CRTSizeX);
    yPosition   	= CLHEP::RandFlat::shoot(CRT_TOP_center[1]-0.5*CRTSizeY,CRT_TOP_center[1]+0.5*CRTSizeY);
    zPosition   	= CLHEP::RandFlat::shoot(CRT_TOP_center[2]-0.5*CRTSizeZ,CRT_TOP_center[2]+0.5*CRTSizeZ);

    if(driftcoordinate==1)
    { //drift in X
      xPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[0]-0.5*CRTSizeX-BufferLengthOnCRTBottom,CRT_BOT_center[0]+0.5*CRTSizeX+BufferLengthOnCRTBottom);
      yPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[1]-0.5*CRTSizeY,CRT_BOT_center[1]+0.5*CRTSizeY);
      zPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[2]-0.5*CRTSizeZ-BufferLengthOnCRTBottom,CRT_BOT_center[2]+0.5*CRTSizeZ+BufferLengthOnCRTBottom);
    }
    if(driftcoordinate==2)
    { //drift in Y
      xPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[0]-0.5*CRTSizeX,CRT_BOT_center[0]+0.5*CRTSizeX);
      yPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[1]-0.5*CRTSizeY-BufferLengthOnCRTBottom,CRT_BOT_center[1]+0.5*CRTSizeY+BufferLengthOnCRTBottom);
      zPositionEnd  	= CLHEP::RandFlat::shoot(CRT_BOT_center[2]-0.5*CRTSizeZ-BufferLengthOnCRTBottom,CRT_BOT_center[2]+0.5*CRTSizeZ+BufferLengthOnCRTBottom);
    }
  }

  if(driftcoordinate==1)
  { //drift in X
    fTH2CRTTop->Fill(zPosition,xPosition);
    fTH2CRTBot->Fill(zPositionEnd,xPositionEnd);
  }
  if(driftcoordinate==2)
  { //drift in Y
    fTH2CRTTop->Fill(zPosition,yPosition);
    fTH2CRTBot->Fill(zPositionEnd,yPositionEnd);
  }
  
  double 	 time        	= 0.;
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
   mf::LogPrint("CRTGen") << "Shooting muon on " << xPosition << " " << yPosition << " " << zPosition <<  "to "<< xPositionEnd << " " << yPositionEnd << " " << zPositionEnd <<" With momentum: " << xMomentum << " " << yMomentum << " " << zMomentum << " E=" << energy << " m=" << mass;

  TLorentzVector pos(xPosition, yPosition, zPosition, time);
  TLorentzVector mom(xMomentum, yMomentum, zMomentum, energy);

  fTH1Energy->Fill(energy);

  simb::MCParticle part(-1, pdg, "primary");
  part.AddTrajectoryPoint(pos, mom);

  truth.Add(part);

  truthcol->push_back(truth);

  e.put(std::move(truthcol));
}

DEFINE_ART_MODULE(evgen::CRTGen)
