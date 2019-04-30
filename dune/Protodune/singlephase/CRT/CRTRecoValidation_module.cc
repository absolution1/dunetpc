////////////////////////////////////////////////////////////////////////
// Class:       CRTReco
// Plugin Type: analyzer (art v2_11_02)
// File:        CRTReco_module.cc
//
// Written by: Richard Diurba
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

//LArSoft includes

#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


//Local includes
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TTree.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TString.h"

//c++ includes
#include <numeric> //std::accumulate was moved from <algorithm> to <numeric> in c++14
#include <iostream>
#include <cmath>
using namespace std;   // Namespaces established to make life easier
using namespace ROOT::Math;
namespace CRT {
  class CRTReco;
}


class CRT::CRTReco: public art::EDAnalyzer {
  public: // Setup functions and variables
    explicit CRTReco(fhicl::ParameterSet
      const & p);
  std::string fTrackModuleLabel = "pmtrack";
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  CRTReco(CRTReco
    const & ) = delete;
  CRTReco(CRTReco && ) = delete;
  CRTReco & operator = (CRTReco
    const & ) = delete;
  CRTReco & operator = (CRTReco && ) = delete;
  void analyze(art::Event
    const & e) override;

  bool moduleMatcher(int module1, int module2);

  void beginJob() override;
  void endJob();
  void createPNG(TH1D * histo);
  double setAngle(double angle);

  int nEvents = 0;
  int nHaloMuons = 0;
  private: ofstream logFile;

  //Parameters for reading in CRT::Triggers and associated AuxDetSimChannels.
  art::InputTag fCRTLabel; //Label for the module that produced 

  TTree * fCRTTree;
  TTree * fCRTTree2;
  double hitDisplX_F=-1000, hitDisplY_F=-1000, hitDisplX_B=-1000, hitDisplY_B=-1000, htimeDiff=0, tDiff_F=0, tDiff_B=0;
  double hitY_B, hitX_B, mccX_B, mccY_B;
  int moduleX_B, moduleY_B, moduleX_F, moduleY_F;
  int moduleTimeDiff, mccCheatTime_F, mccRecoTime_F, mccCheatTime_B, mccRecoTime_B;
	double t_F;


  typedef struct // Structures for arrays to move hits from raw to reco
  {

    int channel;
    int module;
    int adc;
    int triggerTime;
  }
  tempHits;

  typedef struct {
    int tempId;

    double hitPositionX;
    double hitPositionY;
    double hitPositionZ;
    double timeAvg;
    double timeModuleDiff;
    int moduleX;
    int channelX;
    int moduleY;
    int channelY;
  }
  recoHits;

  typedef struct {
    double hitPositionX1;
    double hitPositionY1;
    double hitPositionZ1;

    double hitPositionX2;
    double hitPositionY2;
    double hitPositionZ2;
 
	
  }
  combHits;

  typedef struct // These are displacement metrics for track and hit reco
  {
    int tempId;
    int CRTTrackId;
    int recoId;
    double deltaX;

    double averageSignedDistanceXY;
    double averageSignedDistanceYZ;
    double averageSignedDistanceXZ;
    double averageSignedDistance;
    double deltaAngleYZ;
    double deltaAngleXZ;
    double deltaAngleXY;
    double X1;
    double Y1;
    double Z1;
    double X2;
    double Y2;
    double Z2;

  }
  tracksPair;

  struct removePairIndex // Struct to remove tracks after sorting
  {
    const tracksPair tracksPair1;
    removePairIndex(const tracksPair & tracksPair0): tracksPair1(tracksPair0) {}

    bool operator()(const tracksPair & tracksPair2) {
      return (tracksPair1.recoId == tracksPair2.recoId || tracksPair1.tempId == tracksPair2.tempId);
    }
  };

  struct sortPair // Struct to sort to find best CRT track for TPC track
  {
    bool operator()(const tracksPair & pair1,
      const tracksPair & pair2) {
      return (fabs(pair1.averageSignedDistanceYZ) < fabs(pair2.averageSignedDistanceYZ));
    }
  };
  std::vector < combHits > combTrackHits;
  std::vector < recoHits > primaryHits_F;
  std::vector < recoHits > primaryHits_B;

  std::vector < tempHits > tempHits_F;
  std::vector < tempHits > tempHits_B;
  std::vector < tracksPair > allTracksPair;

};

CRT::CRTReco::CRTReco(fhicl::ParameterSet
    const & p):
  EDAnalyzer(p), fCRTLabel(p.get < art::InputTag > ("CRTLabel")) {
    consumes < std::vector < CRT::Trigger >> (fCRTLabel);
    consumes < std::vector < art::Assns < sim::AuxDetSimChannel, CRT::Trigger >>> (fCRTLabel); // CRT art consumables
  }
//Flipped downstream
bool CRT::CRTReco::moduleMatcher(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 0 && (module2 == 5 || module2 == 4)) || (module1 == 12 && (module2 == 5 || module2 == 4)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 1 && (module2 == 6 || module2 == 7)) || (module1 == 13 && (module2 == 6 || module2 == 7)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23)) || (module1 == 2 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 24 || module2 == 25)) || (module1 == 31 && (module2 == 24 || module2 == 25)) || (module1 == 3 && (module2 == 8 || module2 == 9)) || (module1 == 15 && (module2 == 8 || module2 == 9)) || (module1 == 18 && (module2 == 26 || module2 == 27)) || (module1 == 30 && (module2 == 26 || module2 == 27))) return 1;
  else return 0;

}
/*
bool CRT::CRTReco::moduleMatcher(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 0 && (module2 == 5 || module2 == 4)) || (module1 == 12 && (module2 == 5 || module2 == 4)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 1 && (module2 == 6 || module2 == 7)) || (module1 == 13 && (module2 == 6 || module2 == 7)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23)) || (module1 == 2 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 24 || module2 == 25)) || (module1 == 31 && (module2 == 24 || module2 == 25)) || (module1 == 3 && (module2 == 8 || module2 == 9)) || (module1 == 15 && (module2 == 8 || module2 == 9)) || (module1 == 18 && (module2 == 26 || module2 == 27)) || (module1 == 30 && (module2 == 26 || module2 == 27))) return 1;
  else return 0;

}*/

/* Function for data
bool CRT::CRTReco::moduleMatcher(int module1, int module2){
// Function checking if two hits could reasonably be matched into a 2D hit
if ((module1==5 && (module2==6 || module1==7)) || (module1==4 && (module2==6 || module2==7)) || (module1==3 && (module2==1 || module2==0)) || (module1==2 && (module2==1 || module2==0)) || (module1==13 && (module2==15 || module2==14)) || (module1==13 && (module2==14 || module2==15)) || (module1==12 && (module2==14 || module2==15)) || (module1==11 && (module2==8 || module2==9)) || (module1==10 && (module2==8 || module2==9)) || (module1==18 && (module2==17 || module2==16)) || (module1==19 && (module2==17 || module2==16)) || (module1==20 && (module2==22 || module2==23)) || (module1==21 && (module2==22 || module2==23)) || (module1==29 && (module2==31 || module2==30)) || (module1==28 && (module2==31 || module2==30)) || (module1==27 && (module2==24 || module2==25)) || (module1==26 && (module2==24 || module2==25)))     return 1;
else return 0;


}
*/
void CRT::CRTReco::createPNG(TH1D * histo) {
  //Save important histograms as PNG
  TCanvas * c = new TCanvas;
  TH1D * h = histo;
  h -> Draw();
  TImage * img = TImage::Create();
  img -> FromPad(c);
  img -> WriteImage((std::string(histo -> GetName()) + ".png").c_str());
  delete h;
  delete c;
  delete img;

}
double CRT::CRTReco::setAngle(double angle) {
  if (angle < 0) {
    angle += 3.14159265359;
  }
  angle *= 180.00 / 3.14159265359;
  return angle;
}


void CRT::CRTReco::analyze(art::Event
  const & event) // Analysis module
{

  int nHits = 0;
  art::ServiceHandle < cheat::BackTrackerService > backTracker;
  art::ServiceHandle < cheat::ParticleInventoryService > partInventory;

  const sim::ParticleList & plist = partInventory -> ParticleList();

  // Loop over all the particles
  int nPrimaryMuons = 0;

  primaryHits_F.clear();
  primaryHits_B.clear();
  allTracksPair.clear();
  tempHits_F.clear();
  tempHits_B.clear(); // Arrays to compile hits and move them through

  logFile.open("ProtoDUNE.log"); // Logfile I don't use right now

  //Get triggers
  //cout << "Getting triggers" << endl;
  const auto & triggers = event.getValidHandle < std::vector < CRT::Trigger >> (fCRTLabel);

  art::FindManyP < sim::AuxDetSimChannel > trigToSim(triggers, event, fCRTLabel);

  //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
  art::ServiceHandle < geo::Geometry > geom;

for (unsigned int i=0; i<32; i++){
const auto & trigGeo = geom -> AuxDet(i);
const auto & hit1Geo = trigGeo.SensitiveVolume(0);
for (unsigned int j=0; j<32; j++){ 
const auto & trigGeo2 = geom -> AuxDet(j);
const auto & hit2Geo = trigGeo2.SensitiveVolume(0);
if ((i==2 || i==3 || i==14 || i==15 || i==18 || i==19 || i==30 || i==31) && ( j==2 || j==3 || j==14 || j==15 || j==18 || j==19 || j==30 || j==31)){
      const auto hit1Center = hit1Geo.GetCenter();
      const auto hit2Center = hit2Geo.GetCenter();
	cout<<i<<','<<j<<','<<hit1Center.Y()-hit2Center.Y()<<endl;

//cout<<"Module "<<i<<" Channel "<<j<<" Position X: "<<hit1Center.X()<<" Position Y: "<<hit1Center.Y()<<" Position Z: "<<hit1Center.Z()<<endl;
}
//cout<<"Module "<<i<<" Position X: "<<hit2Center.X()<<" Position Y: "<<hit2Center.Y()<<" Position Z: "<<hit2Center.Z()<<endl;
}
}



  //Mapping from channel to trigger
  std::unordered_map < size_t, double > prevTimes;
  int hitID = 0;
  //cout << "Looking for hits in Triggers" << endl;

  for (const auto & trigger: * triggers) {
    const auto & hits = trigger.Hits();
    for (const auto & hit: hits) { // Collect hits on all modules

      if (hit.ADC() > 600) { // Keep if they are above threshold

        tempHits tHits;

        tHits.module = trigger.Channel(); // Values to add to array
        tHits.channel = hit.Channel();
        tHits.adc = hit.ADC();
	tHits.triggerTime=trigger.Timestamp();
	//cout<<trigger.Timestamp()<<endl;
        nHits++;

        const auto & trigGeo = geom -> AuxDet(trigger.Channel()); // Get geo  
        const auto & csens = trigGeo.SensitiveVolume(hit.Channel());
        const auto center = csens.GetCenter();
        if (center.Z() < 100) tempHits_F.push_back(tHits); // Sort F/B from Z
        else tempHits_B.push_back(tHits);
        hitID++;
      }
    }
  }
  nEvents++;
  //cout << "Hits compiled for event: " << nEvents << endl;
  //cout << "Number of Hits above Threshold:  " << hitID << endl;

  for (unsigned int f = 0; f < tempHits_F.size(); f++) {
    for (unsigned int f_test = 0; f_test < tempHits_F.size(); f_test++) {
      const auto & trigGeo = geom -> AuxDet(tempHits_F[f].module);
      const auto & trigGeo2 = geom -> AuxDet(tempHits_F[f_test].module);
      const auto & hit1Geo = trigGeo.SensitiveVolume(tempHits_F[f].channel);
      const auto hit1Center = hit1Geo.GetCenter();
      // Create 2D hits from geo of the Y and X modules
      const auto & hit2Geo = trigGeo2.SensitiveVolume(tempHits_F[f_test].channel);
      const auto hit2Center = hit2Geo.GetCenter();

      if (moduleMatcher(tempHits_F[f_test].module, tempHits_F[f].module) == 1) {
        // Get the center of the hits (CRT_Res=2.5 cm)
        double hitX = hit1Center.X();
        double hitY = hit2Center.Y();
        double hitZ = (hit1Center.Z() + hit2Center.Z()) / 2.f;

        recoHits rHits;

        rHits.hitPositionX = hitX;
        rHits.hitPositionY = hitY;
        rHits.hitPositionZ = hitZ;
	rHits.moduleX=tempHits_F[f].module;
	rHits.moduleY=tempHits_F[f_test].module;
	rHits.timeAvg = (tempHits_F[f_test].triggerTime+tempHits_F[f].triggerTime)/2.f;
	rHits.timeModuleDiff = fabs(tempHits_F[f_test].triggerTime-tempHits_F[f].triggerTime)/1.f;
	rHits.channelX=tempHits_F[f].channel;
	rHits.channelY=tempHits_F[f_test].channel;
	if (fabs(tempHits_F[f_test].triggerTime-tempHits_F[f].triggerTime)<2) primaryHits_F.push_back(rHits); // Add array
      }
    }
  }
  for (unsigned int f = 0; f < tempHits_B.size(); f++) {
    for (unsigned int f_test = 0; f_test < tempHits_B.size(); f_test++) { // Same as above but for back CRT


      const auto & trigGeo = geom -> AuxDet(tempHits_B[f].module);
      const auto & trigGeo2 = geom -> AuxDet(tempHits_B[f_test].module);
      const auto & hit1Geo = trigGeo.SensitiveVolume(tempHits_B[f].channel);
      const auto hit1Center = hit1Geo.GetCenter();

      const auto & hit2Geo = trigGeo2.SensitiveVolume(tempHits_B[f_test].channel);
      const auto hit2Center = hit2Geo.GetCenter();

      if (moduleMatcher(tempHits_B[f_test].module, tempHits_B[f].module) == 1) {

        double hitX = hit1Center.X();
        double hitY = hit2Center.Y();
        double hitZ = (hit1Center.Z() + hit2Center.Z()) / 2.f;

        recoHits rHits;

        rHits.hitPositionX = hitX;
        rHits.hitPositionY = hitY;
        rHits.hitPositionZ = hitZ;
	rHits.moduleX=tempHits_B[f].module;
	rHits.moduleY=tempHits_B[f_test].module;
	rHits.timeAvg = (tempHits_B[f_test].triggerTime+tempHits_B[f].triggerTime)/2.f;
	rHits.timeModuleDiff = fabs(tempHits_B[f_test].triggerTime-tempHits_B[f].triggerTime)/1.f;
	rHits.channelX=tempHits_B[f].channel;
	rHits.channelY=tempHits_B[f_test].channel;
        if (fabs(tempHits_B[f_test].triggerTime-tempHits_B[f].triggerTime)<2) primaryHits_B.push_back(rHits);
      }
    }
  }


  for (const auto & PartPair: plist) {
    const simb::MCParticle & particle = * (PartPair.second);

    // Locate the MCC muon and compare hit to CRT hit if a CRT is reasonably nearby

    if (particle.Process() == "primary" && fabs(particle.PdgCode()) == fabs(13) && (particle.Position(0).Z() < -250  && particle.EndPosition().Z() > 10) && fabs(particle.Position(0).X()) < 400 && fabs(particle.Position(0).Y()) < 600 &&
      fabs(particle.EndPosition().X()) < 400 && fabs(particle.EndPosition().Y()) < 600) {
      nPrimaryMuons++;
const TLorentzVector& primaryMomentumStart = particle.Momentum(0);
	double trackSlopeXZ=primaryMomentumStart.X()/primaryMomentumStart.Z();
	double trackSlopeYZ=primaryMomentumStart.Y()/primaryMomentumStart.Z();
    //unsigned int  mccTrackPoints=particle.NumberTrajectoryPoints();
	int mccTime_F=particle.T(0);
	//int mccTi
//int mccTime_F_close; 
//int mccTime_F_far;
 /*       for (unsigned int i=0; i<mccTrackPoints; i++){
        cout<<particle.T(0)<<','<<particle.T(i)<<endl;
	//if (fabs(particle.Position(i).Z()-267)<30) mccTime_F_close=100*(particle.T(i)-mccTimeVsTPC);
	if (fabs(particle.Position(i).Z()-971)<30) mccTime_F_far=100*(particle.T(i)-mccTimeVsTPC);

	if ((particle.Position(i).Z())>-250) break; }
*/
// Compare CRT hit with MCC
     //int mccTime_F=0;
      hitDisplX_F = -1000, hitDisplY_F = -1000;
      for (unsigned int f = 0; f < primaryHits_F.size(); f++) {
        double hitX = primaryHits_F[f].hitPositionX;
        double hitY = primaryHits_F[f].hitPositionY;
        double hitZ = primaryHits_F[f].hitPositionZ;
	//mccTime_F=mccTime_F_close;
	//if (hitZ<-900) mccTime_F=mccTime_F_far;	
        double mccHitX = trackSlopeXZ*(hitZ-particle.Position(0).Z())+particle.Position(0).X() ;
        double mccHitY = trackSlopeYZ*(hitZ-particle.Position(0).Z())+particle.Position(0).Y();
        if ((fabs(mccHitX - hitX) < 999 && fabs(mccHitY - hitY) < 999) && (pow(pow(fabs(mccHitX - hitX),2)+pow(fabs(mccHitY - hitY),2),0.5) < pow(pow(hitDisplY_F,2)+pow(hitDisplX_F,2),0.5))) {
          hitDisplX_F = mccHitX - hitX;
          hitDisplY_F = mccHitY - hitY;
	  //cout<<"MCC Front Timing : "<<t_F<<','<<mccHitY<<endl;
	  moduleX_F=primaryHits_F[f].moduleX;
	  moduleY_F=primaryHits_F[f].moduleY;
	  t_F=primaryHits_F[f].timeAvg;
	  mccRecoTime_F=t_F;	  
	  tDiff_F=primaryHits_F[f].timeModuleDiff;
	  mccCheatTime_F=mccTime_F;
        }
      }
  if (fabs(hitDisplX_F)<20 && fabs(hitDisplY_F)<20) fCRTTree->Fill();
    } if (particle.Process() == "primary" && fabs(particle.PdgCode()) == fabs(13) &&  (particle.Position(0).Z() < -200  && particle.EndPosition().Z() > 1200) && fabs(particle.Position(0).X()) < 400 && fabs(particle.Position(0).Y()) < 600 &&
      fabs(particle.EndPosition().X()) < 400 && fabs(particle.EndPosition().Y()) < 600) {
      nPrimaryMuons++;
      nHaloMuons++;

      cout << "Total number of Back Cosmic Muons: " << nHaloMuons << endl;
	int nTrajectory=particle.NumberTrajectoryPoints();

	int approxExit;
	// Find Exit
	for (int i=0; i<nTrajectory-2; i++){
	approxExit=i;
	if (particle.Position(i).Z()>700) break;
	}
  	cout<<"Approximate Exiting Position: "<<particle.Position(approxExit).X()<<','<<particle.Position(approxExit).Y()<<','<<particle.Position(approxExit).Z()<<endl;
	const TLorentzVector& primaryMomentumEnd = particle.Momentum(approxExit);
        double trackSlopeXZ=primaryMomentumEnd.X()/primaryMomentumEnd.Z();
	double trackSlopeYZ=primaryMomentumEnd.Y()/primaryMomentumEnd.Z();

// With the track starting before going into the TPC, compare where it would go when it hits the back CRT
      hitDisplX_B = -1000000, hitDisplY_B= -100000;
      htimeDiff=100000;
      if (primaryHits_B.size()<1) break;

       
      for (unsigned int f = 0; f < primaryHits_B.size(); f++) {
        double hitX = primaryHits_B[f].hitPositionX;
        double hitY = primaryHits_B[f].hitPositionY;
        double hitZ = primaryHits_B[f].hitPositionZ;

    double mccHitX = trackSlopeXZ*(hitZ-particle.Position(approxExit).Z())+particle.Position(approxExit).X() ;
        double mccHitY = trackSlopeYZ*(hitZ-particle.Position(approxExit).Z())+particle.Position(approxExit).Y();

        if ((fabs(mccHitX - hitX) < 999 && fabs(mccHitY - hitY) < 999) && (pow(pow(fabs(mccHitX - hitX),2)+pow(fabs(mccHitY - hitY),2),0.5) < pow(pow(hitDisplY_B,2)+pow(hitDisplX_B,2),0.5))) {
          hitDisplX_B = mccHitX - hitX;
          hitDisplY_B = mccHitY - hitY;
	  //cout<<"Channel No."<<primaryHits_B[f].channelX<<','<<primaryHits_B[f].channelY<<endl;
	 // cout<<nEvents<<endl;
	 // cout<<"MCC Back: "<<mccHitX<<','<<mccHitY<<endl;
	  mccCheatTime_B=particle.T(approxExit);
	  mccRecoTime_B=primaryHits_B[f].timeAvg;
	  moduleTimeDiff=primaryHits_B[f].timeAvg-t_F;
	  moduleX_B=primaryHits_B[f].moduleX;
	  moduleY_B=primaryHits_B[f].moduleY;
	  cout<<moduleTimeDiff<<endl;
        }
      }
      fCRTTree2->Fill();
   
}
  }
}
// CRT Metrics
void CRT::CRTReco::beginJob() {
	art::ServiceHandle<art::TFileService> fileServiceHandle;
       fCRTTree = fileServiceHandle->make<TTree>("Front CRT Match", "event by event info");
	fCRTTree->Branch("event", &nEvents, "fnEvents/I");
	fCRTTree->Branch("hitDisplX_F", &hitDisplX_F, "hitDisplX_F/D");
	fCRTTree->Branch("hitDisplY_F", &hitDisplY_F, "hitDisplY_F/D");
	fCRTTree->Branch("htimeModuleDiff_F", &tDiff_F, "tDiff_F/D");
	fCRTTree->Branch("moduleX_F", &moduleX_F, "moduleX_F/I");
	fCRTTree->Branch("moduleY_F", &moduleY_F, "moduleY_F/I");       
	fCRTTree->Branch("mccCheatTime_F", &mccCheatTime_F, "mccCheatTime_F/I");
	fCRTTree->Branch("mccRecoTime_F", &mccRecoTime_F, "mccRecoTime_F/I");
fCRTTree2 = fileServiceHandle->make<TTree>("Two CRTs Match", "event by event info");
       fCRTTree2->Branch("event", &nEvents, "fnEvents/I");
	fCRTTree2->Branch("hitDisplX_B", &hitDisplX_B, "hitDisplX_B/D");
	fCRTTree2->Branch("hitDisplY_B", &hitDisplY_B, "hitDisplY_B/D");
	fCRTTree2->Branch("moduleTimeDiff",&moduleTimeDiff, "moduleTimeDiff/I");
	fCRTTree2->Branch("moduleX_B", &moduleX_B, "moduleX_B/I");
	fCRTTree2->Branch("moduleY_B", &moduleY_B, "moduleY_B/I");
	fCRTTree2->Branch("mccCheatTime_B", &mccCheatTime_B, "mccCheatTime_B/I");
	fCRTTree2->Branch("mccRecoTime_B", &mccRecoTime_B, "mccRecoTime_B/I");

}
void CRT::CRTReco::endJob() // Create PNG plots after completion
{



}














DEFINE_ART_MODULE(CRT::CRTReco)
