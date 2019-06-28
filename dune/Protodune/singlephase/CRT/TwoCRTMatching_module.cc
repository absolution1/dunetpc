////////////////////////////////////////////////////////////////////////
// Class:       TwoCRTMatching
// Plugin Type: Analyzer (art v2_11_02)
// File:        TwoCRTMatching_module.cc
//
// Written by: Richard Diurba
// Adapted from MCC Code by: Arbin Timilsina 
// CRT Trigger Architecture by: Andrew Olivier 
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
//#include "art/Framework/Core/EDAnalyzer.h"
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
#include "art_root_io/TFileService.h"

//LArSoft includes

#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"


#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"


#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"


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
  class TwoCRTMatching;
}


class CRT::TwoCRTMatching: public art::EDAnalyzer {
  public: // Setup functions and variables
    explicit TwoCRTMatching(fhicl::ParameterSet
      const & p);
  std::string fTrackModuleLabel = "pandoraTrack";
  //std::string fTrackModuleLabel = "pmtrack";
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  TwoCRTMatching(TwoCRTMatching
    const & ) = delete;
  TwoCRTMatching(TwoCRTMatching && ) = delete;
  TwoCRTMatching & operator = (TwoCRTMatching
    const & ) = delete;
  TwoCRTMatching & operator = (TwoCRTMatching && ) = delete;
  void analyze(art::Event
    const & e) override;
  //void analyze(art::Event
  //  const & e) override;
// Declare functions and variables for validation
  bool moduleMatcherData(int module1, int module2);
  bool moduleMatcherMCC(int module1, int module2);
  double signedPointToLineDistance(double firstPoint1, double firstPoint2,   double secondPoint1, double secondPoint2, double trackPoint1, double trackPoint2);
  double signed3dDistance(double firstPoint1, double firstPoint2, double firstPoint3, double secondPoint1, double secondPoint2, double secondPoint3, TVector3 trackPos);
  void beginJob() override;
  void endJob() override;
  void createPNG(TH1D * histo);
  double setAngle(double angle);
int moduletoCTB(int module2, int module1);
  int nEvents = 0;
  int nHaloMuons = 0;
  int nHitsPerEvent=0;
  int nMCCMuons=0;
  private: ofstream logFile;

  //Parameters for reading in CRT::Triggers and associated AuxDetSimChannels.
  art::InputTag fCRTLabel; //Label for the module that analyzed 
  art::InputTag fCTBLabel;
  TTree * fCRTTree;
  TTree * fMCCMuon;
    int run, subRun;
    bool fMCCSwitch;
    bool fCTBTriggerOnly;
    bool fSCECorrection;
    bool fModuleSwitch;
    int fADCThreshold;
    int fModuletoModuleTimingCut;
    int fFronttoBackTimingCut;


  //double averageSignedDistanceXY;
    double averageSignedDistanceYZ;
    double averageSignedDistanceXZ;
    double averageSignedDistance;
    double displAngleXZ;
    double displAngleYZ;
    double displAngleXY;
    double CRT_TOF;
    double deltaX_F;
    double deltaX_B;
    double deltaY_F;
    double deltaY_B;
    double dotCos;
    int adcX_F, adcX_B, adcY_F, adcY_B;
    double X_F, X_B, Y_F, Y_B, Z_F, Z_B;
    double trackX1, trackX2, trackY1, trackY2, trackZ1, trackZ2;
    int moduleX_F, moduleX_B, moduleY_F, moduleY_B;
    int stripX_F, stripX_B, stripY_F, stripY_B;
    double measuredT0;
    double measuredXOffset;
    bool recoPandoraT0Check;
    int trackID;
    double truthEnergy;
    int mccTrackId;
  typedef struct // Structures for arrays to move hits from raw to reco to validation
  {

    int channel;
    int module;
    int channelGeo;
    int adc;
    int triggerTime;
  }
  tempHits;

  typedef struct {
    int tempId;	
    int adcX;
    int adcY;

    double hitPositionX;
    double hitPositionY;
    double hitPositionZ;
    double timeAvg;
    int geoX;
    int geoY;
    int stripX;
    int stripY;
  }
  recoHits;

  typedef struct {
    int adcX1;
    int adcX2;
    int adcY1;
    int adcY2;
    double hitPositionX1;
    double hitPositionY1;
    double hitPositionZ1;
    double t0;
    double hitPositionX2;
    double hitPositionY2;
    double hitPositionZ2;
    double timeDiff;
    int moduleX1, moduleX2, moduleY1, moduleY2;
    int stripX1, stripY1, stripY2, stripX2;

  }
  combHits;

  typedef struct // These are displacement metrics for track and hit reco
  {
    int tempId;
    int CRTTrackId;
    int recoId;
    int adcX1;
    int adcY1;
    int adcX2;
    int adcY2;
    double deltaX_F;
    double deltaX_B;
    double deltaY_F;
    double deltaY_B;
    double deltaX;
    double deltaY;
    double averageSignedDistanceXY;
    double averageSignedDistanceYZ;
    double averageSignedDistanceXZ;
    double averageSignedDistance;
    double deltaAngleYZ;
    double deltaAngleXZ;
    double deltaAngleXY;
    double dotProductCos;
    double X1;
    double Y1;
    double Z1;
    double X2;
    double Y2;
    double Z2;
    int trackID;
    double timeDiff;
    double t0;
    double xOffset;
    bool pandoraT0Check;
    TVector3 trackStartPosition;
    TVector3 trackEndPosition;
    int moduleX1, moduleX2, moduleY1, moduleY2;
    int stripX1, stripX2, stripY1, stripY2;
    
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

	return (fabs(pair1.dotProductCos)>fabs(pair2.dotProductCos));
        //return (fabs(pair1.deltaX_F)+fabs(pair1.deltaY_F)+fabs(pair1.deltaX_B)+fabs(pair1.deltaY_B)<fabs(pair2.deltaX_F)+fabs(pair2.deltaY_F)+fabs(pair2.deltaX_B)+fabs(pair2.deltaY_B));
	//return ((fabs(pair1.dotProductCos)>.998 && pair1.deltaY<pair2.deltaY && pair1.deltaX<pair2.deltaX));

  }
  };
  std::vector < combHits > combTrackHits;
  std::vector < recoHits > primaryHits_F;
  std::vector < recoHits > primaryHits_B;

  std::vector < tempHits > tempHits_F;
  std::vector < tempHits > tempHits_B;
  std::vector < tracksPair > allTracksPair;

};

CRT::TwoCRTMatching::TwoCRTMatching(fhicl::ParameterSet
    const & p):
  EDAnalyzer(p), fCRTLabel(p.get < art::InputTag > ("CRTLabel")),  fCTBLabel(p.get<art::InputTag>("CTBLabel")) {
    consumes < std::vector < CRT::Trigger >> (fCRTLabel);
    consumes < std::vector < art::Assns < sim::AuxDetSimChannel, CRT::Trigger >>> (fCRTLabel); 
  fMCCSwitch=(p.get<bool>("MCC"));
  fCTBTriggerOnly=(p.get<bool>("CTBOnly"));
  fSCECorrection=(p.get<bool>("SCECorrection"));
  }


//Displacement functions for validation
double CRT::TwoCRTMatching::signedPointToLineDistance(double firstPoint1, double firstPoint2, double secondPoint1, double secondPoint2, double trackPoint1, double trackPoint2){
        double numerator = (secondPoint2-firstPoint2)*trackPoint1 - (secondPoint1-firstPoint1) * trackPoint2 + secondPoint1*firstPoint2 - firstPoint1*secondPoint2; //removed the absolute value here, so will give signed distance //the sign indicates a right-hand ruled sign: cross product of line-to-point vector and the direction of the vector (in that order) gives the sign of the result
        double denominator = sqrt( (secondPoint2-firstPoint2)*(secondPoint2-firstPoint2) + (secondPoint1-firstPoint1)*(secondPoint1-firstPoint1) );

        return numerator/denominator;

}
double CRT::TwoCRTMatching::signed3dDistance(double firstPoint1, double firstPoint2, double firstPoint3, double secondPoint1, double secondPoint2, double secondPoint3, TVector3 point){

double denominator = sqrt( (secondPoint2-firstPoint2)*(secondPoint2-firstPoint2) + (secondPoint1-firstPoint1)*(secondPoint1-firstPoint1)+ (secondPoint3-firstPoint3)*(secondPoint3-firstPoint3));

double X1=point.X()-firstPoint1;
double Y1=point.Y()-firstPoint2;
double Z1=point.Z()-firstPoint3;

double X2=point.X()-secondPoint1;
double Y2=point.Y()-secondPoint2;
double Z2=point.Z()-secondPoint3;

double numerator=(X1*Y2-Y1*X2)-(X1*Z2-X2*Z1)+(Y1*Z2-Z1*Y2);

return numerator/denominator;




}


// v6 Geo Channel Map
bool CRT::TwoCRTMatching::moduleMatcherMCC(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 6 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 26 || module2 == 27)) || (module1 == 31 && (module2 == 26 || module2 == 27)) || (module1 == 7 && (module2 == 12 || module2 == 13)) || (module1 == 15 && (module2 == 12 || module2 == 13)) || (module1 == 18 && (module2 == 24 || module2 == 25)) || (module1 == 30 && (module2 == 24 || module2 == 25)) || (module1 == 1 && (module2 == 4 || module2 == 5)) || (module1 == 9 && (module2 == 4 || module2 == 5)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 0 && (module2 == 2 || module2 == 3)) || (module1 == 8 && (module2 == 2 || module2 == 3)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23))) return 1;
  else return 0;

}

bool CRT::TwoCRTMatching::moduleMatcherData(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 6 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 26 || module2 == 27)) || (module1 == 31 && (module2 == 26 || module2 == 27)) || (module1 == 7 && (module2 == 12 || module2 == 13)) || (module1 == 15 && (module2 == 12 || module2 == 13)) || (module1 == 18 && (module2 == 24 || module2 == 25)) || (module1 == 30 && (module2 == 24 || module2 == 25)) || (module1 == 1 && (module2 == 4 || module2 == 5)) || (module1 == 9 && (module2 == 4 || module2 == 5)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 0 && (module2 == 2 || module2 == 3)) || (module1 == 8 && (module2 == 2 || module2 == 3)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23))) return 1;
  else return 0;

}

/* v5 Geo
// Function to match CRT modules below is for MCC and the 2nd is for data
bool CRT::TwoCRTMatching::moduleMatcherMCC(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 0 && (module2 == 5 || module2 == 4)) || (module1 == 12 && (module2 == 5 || module2 == 4)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 1 && (module2 == 6 || module2 == 7)) || (module1 == 13 && (module2 == 6 || module2 == 7)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23)) || (module1 == 2 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 24 || module2 == 25)) || (module1 == 31 && (module2 == 24 || module2 == 25)) || (module1 == 3 && (module2 == 8 || module2 == 9)) || (module1 == 15 && (module2 == 8 || module2 == 9)) || (module1 == 18 && (module2 == 26 || module2 == 27)) || (module1 == 30 && (module2 == 26 || module2 == 27))) return 1;
  else return 0;

}

// Function to match CRT modules below is for MCC and the 2nd is for data
bool CRT::TwoCRTMatching::moduleMatcherData(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 0 && (module2 == 5 || module2 == 4)) || (module1 == 12 && (module2 == 5 || module2 == 4)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 1 && (module2 == 6 || module2 == 7)) || (module1 == 13 && (module2 == 6 || module2 == 7)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23)) || (module1 == 2 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 24 || module2 == 25)) || (module1 == 31 && (module2 == 24 || module2 == 25)) || (module1 == 3 && (module2 == 8 || module2 == 9)) || (module1 == 15 && (module2 == 8 || module2 == 9)) || (module1 == 18 && (module2 == 26 || module2 == 27)) || (module1 == 30 && (module2 == 26 || module2 == 27))) return 1;
  else return 0;

}*/

int CRT::TwoCRTMatching::moduletoCTB(int module2, int module1){
  if (module1 == 14 && module2 == 11 ) return 15;
  else if (module1 == 14 &&  module2 == 10) return 10;
  else if (module1 == 6 &&  module2 == 11) return 8;
  else if (module1 == 6 &&  module2 == 10) return 9;
  else if (module1 == 18 &&  module2 == 25) return 4;
  else if (module1 == 18 &&  module2 == 24) return 13;
  else if (module1 == 30 &&  module2 == 25) return 3;
  else if (module1 == 30 &&  module2 == 24) return 2;
  else if (module1 == 31 &&  module2 == 27) return 1;
  else if (module1 == 31 &&  module2 == 26) return 0;
  else if (module1 == 19 &&  module2 == 27) return 12;
  else if (module1 == 19 &&  module2 == 26) return 11;
  else if (module1 == 7  &&  module2 == 12) return 7;
  else if (module1 == 7 &&  module2 == 13) return 6;
  else if (module1 == 15  &&  module2 == 12) return 14;
  else if (module1 == 15 &&  module2 == 13) return 5;
  else if (module1 == 1 &&  module2 == 4) return 25;
  else if (module1 == 1 &&  module2 == 5) return 24;
  else if (module1 == 9 &&  module2 == 4) return 26;
  else if (module1 == 9 &&  module2 == 5) return 31;
  else if (module1 == 16 &&  module2 == 20) return 27;
  else if (module1 == 16 &&  module2 == 21) return 28;
  else if (module1 == 28 &&  module2 == 20) return 16;
  else if (module1 == 28 &&  module2 == 21) return 17;
  else if (module1 == 29 &&  module2 == 22) return 18;
  else if (module1 == 29 &&  module2 == 23) return 19;
  else if (module1 == 17 &&  module2 == 22) return 29;
  else if (module1 == 17 &&  module2 == 23) return 20;
  else if (module1 == 8 &&  module2 == 2) return 30;
  else if (module1 == 8 &&  module2 == 3) return 21;
  else if (module1 == 0 &&  module2 == 2) return 23;
  else if (module1 == 0 &&  module2 == 3) return 22;
  else return -1;
}
void CRT::TwoCRTMatching::createPNG(TH1D * histo) {
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
double CRT::TwoCRTMatching::setAngle(double angle) {
  if (angle < 0) {
    angle += 3.14159265359;
  }
  angle *= 180.00 / 3.14159265359;
  return angle;
}


void CRT::TwoCRTMatching::analyze(art::Event
  const & event) // Analysis module

{

    nEvents++;	
    run=event.run();
    subRun=event.subRun();
  mccTrackId=-1;
  if (fMCCSwitch){
    fModuleSwitch=1;
    fADCThreshold=200;
    fModuletoModuleTimingCut=4;
    fFronttoBackTimingCut=100;
    
}
  else {
    fModuleSwitch=0;
    fADCThreshold=10;
    fModuletoModuleTimingCut=5;
    fFronttoBackTimingCut=8;


} 
   if(!fMCCSwitch){
   //const auto& pdspctbs = event.getValidHandle<std::vector<raw::ctb::pdspctb>>(fCTB_tag);
	art::ValidHandle<std::vector<raw::RDTimeStamp>> timingHandle = event.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");

	const raw::RDTimeStamp& timeStamp = timingHandle->at(0);
	if (fCTBTriggerOnly){
	if(timeStamp.GetFlags()!= 13) return;}
		for(const auto& time: *timingHandle)
      {	cout<<time.GetTimeStamp()<<endl;
	}
  }
  int nHits = 0;




	//Detector properties service
	auto const* detectorPropertiesService = lar::providerFrom<detinfo::DetectorPropertiesService>();

  primaryHits_F.clear();
  primaryHits_B.clear();
  allTracksPair.clear();
  tempHits_F.clear();
  tempHits_B.clear(); // Arrays to compile hits and move them through
  combTrackHits.clear();
    //allTracksPair.clear();
  logFile.open("ProtoDUNE.log"); // Logfile I don't use right now

  //Get triggers
  cout << "Getting triggers" << endl;
  const auto & triggers = event.getValidHandle < std::vector < CRT::Trigger >> (fCRTLabel);

  art::FindManyP < sim::AuxDetSimChannel > trigToSim(triggers, event, fCRTLabel);

  //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
  art::ServiceHandle < geo::Geometry > geom;

  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  //Mapping from channel to trigger
  std::unordered_map < size_t, double > prevTimes;
  int hitID = 0;
  cout << "Looking for hits in Triggers" << endl;

  for (const auto & trigger: * triggers) {
    const auto & hits = trigger.Hits();
    for (const auto & hit: hits) { // Collect hits on all modules
      if (hit.ADC() > fADCThreshold) { // Keep if they are above threshold

        tempHits tHits;
	if (!fMCCSwitch){
	art::ValidHandle<std::vector<raw::RDTimeStamp>> timingHandle = event.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");

        tHits.module = trigger.Channel(); // Values to add to array
        tHits.channelGeo = hit.Channel();
	tHits.channel=hit.Channel();
        tHits.adc = hit.ADC();
	tHits.triggerTime=trigger.Timestamp()-timingHandle->at(0).GetTimeStamp();
	}
	else{
        tHits.module = trigger.Channel(); // Values to add to array
        tHits.channelGeo = hit.Channel();
	tHits.channel=hit.Channel();
        tHits.adc = hit.ADC();
	tHits.triggerTime=trigger.Timestamp();
	}
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

  nHitsPerEvent=nHits;
  cout << "Hits compiled for event: " << nEvents << endl;
  cout << "Number of Hits above Threshold:  " << hitID << endl;

  for (unsigned int f = 0; f < tempHits_F.size(); f++) {
    for (unsigned int f_test = 0; f_test < tempHits_F.size(); f_test++) {
      const auto & trigGeo = geom -> AuxDet(tempHits_F[f].module);
      const auto & trigGeo2 = geom -> AuxDet(tempHits_F[f_test].module);
	int flipChannel=tempHits_F[f].channelGeo;
	int flipX=1;
	// v5 geo fixes
	//if (tempHits_F[f].module==21 && !fMCCSwitch){flipX=-1; flipChannel=flipChannel^63;}
	//if (!fMCCSwitch && (tempHits_F[f_test].module==13 || tempHits_F[f_test].module==1)) {flipY=-1; flipChannel=flipChannel^63;}
	//cout<<"Channel flip: "<<flipChannel<<','<<tempHits_F[f_test].channelGeo;
	//if (fabs(tempHits_F[f].triggerTime-tempHits_F[f_test].triggerTime)<1 && tempHits_F[f].module!=tempHits_F[f_test].module){
	//cout<<tempHits_F[f].module<<','<<tempHits_F[f_test].module<<endl;}
      const auto & hit1Geo = trigGeo.SensitiveVolume(flipChannel);
      const auto hit1Center = hit1Geo.GetCenter();
      // Create 2D hits from geo of the Y and X modules
	flipChannel=tempHits_F[f_test].channelGeo;
	int flipY=1;
       const auto & hit2Geo = trigGeo2.SensitiveVolume(flipChannel);
      const auto hit2Center = hit2Geo.GetCenter();
      bool moduleMatched;
      if(fModuleSwitch) moduleMatched=moduleMatcherMCC(tempHits_F[f_test].module, tempHits_F[f].module);
      else moduleMatched=moduleMatcherData(tempHits_F[f_test].module, tempHits_F[f].module);

      if (moduleMatched) {
        // Get the center of the hits (CRT_Res=2.5 cm)
        double hitX = hit1Center.X();
	for (unsigned int a = 0; a < tempHits_F.size(); a++)
	{
	if(tempHits_F[a].module==tempHits_F[f].module && (tempHits_F[a].channelGeo-flipX)==tempHits_F[f].channelGeo) hitX=hit1Center.X()+1.25;
	}
	double hitYPrelim=hit2Center.Y();
	for (unsigned int a = 0; a < tempHits_F.size(); a++)
	{
	if(tempHits_F[a].module==tempHits_F[f_test].module && (tempHits_F[a].channelGeo-flipY)==tempHits_F[f_test].channelGeo) hitYPrelim=hit2Center.Y()+1.25;
	}
	

	
	double hitY=hitYPrelim;
        double hitZ = (hit1Center.Z() + hit2Center.Z()) / 2.f;

        recoHits rHits;
        rHits.adcX=tempHits_F[f].adc;
	rHits.adcY=tempHits_F[f_test].adc;
        rHits.hitPositionX = hitX;
        rHits.hitPositionY = hitY;
        rHits.hitPositionZ = hitZ;
	rHits.geoX=tempHits_F[f].module;
	rHits.geoY=tempHits_F[f_test].module;

	rHits.stripX=tempHits_F[f].channel;
	rHits.stripY=tempHits_F[f_test].channel;
	rHits.timeAvg = (tempHits_F[f_test].triggerTime+tempHits_F[f].triggerTime)/2.0;
	if (fabs(tempHits_F[f_test].triggerTime-tempHits_F[f].triggerTime)<fModuletoModuleTimingCut) primaryHits_F.push_back(rHits); // Add array
    }
    }
  }
  for (unsigned int f = 0; f < tempHits_B.size(); f++) {
    for (unsigned int f_test = 0; f_test < tempHits_B.size(); f_test++) { // Same as above but for back CRT
	int channelFlipCheck=tempHits_B[f].module;
	/* Code to fix v5 geo issues in data
	if (!fMCCSwitch){   
	if (channelFlipCheck==8) channelFlipCheck=11;
	else if (channelFlipCheck==11) channelFlipCheck=8;
         else if (channelFlipCheck==10) channelFlipCheck=9;
	else if (channelFlipCheck==9) channelFlipCheck=10;

	else if (channelFlipCheck==26) channelFlipCheck=25;
        else if (channelFlipCheck==25) channelFlipCheck=26;
	else if (channelFlipCheck==24) channelFlipCheck=27;
	else if (channelFlipCheck==27) channelFlipCheck=24;
	}

     if (!fMCCSwitch && (tempHits_B[f].module==25 || tempHits_B[f].module==11 || tempHits_B[f].module==24 || tempHits_B[f].module==10)){flipX=-1; flipChannel=flipChannel^63;}

	//if (!fMCCSwitch && (tempHits_B[f_test].module==2 || tempHits_B[f_test].module==3  || tempHits_B[f_test].module==14 || tempHits_B[f_test].module==15)) {flipY=-1; flipChannel=flipChannel^63;}
	*/
     int flipX=1;
     int flipY=1;
     int flipChannel=tempHits_B[f].channelGeo;

 
      const auto & trigGeo = geom -> AuxDet(channelFlipCheck);
      const auto & trigGeo2 = geom -> AuxDet(tempHits_B[f_test].module);
      const auto & hit1Geo = trigGeo.SensitiveVolume(flipChannel);
      const auto hit1Center = hit1Geo.GetCenter();
      flipChannel=tempHits_B[f_test].channelGeo;

      const auto & hit2Geo = trigGeo2.SensitiveVolume(flipChannel);
      const auto hit2Center = hit2Geo.GetCenter();
      bool moduleMatched;
      if(fModuleSwitch) moduleMatched=moduleMatcherMCC(tempHits_B[f_test].module, tempHits_B[f].module);
      else moduleMatched=moduleMatcherData(tempHits_B[f_test].module, tempHits_B[f].module);

      if (moduleMatched) {
        double hitX = hit1Center.X();
	
	
	for (unsigned int a = 0; a < tempHits_B.size(); a++)
	{
	if(tempHits_B[a].module==tempHits_B[f].module && (tempHits_B[a].channelGeo-flipX)==tempHits_B[f].channelGeo) hitX=hit1Center.X()+1.25;
	}
	
        double hitYPrelim = hit2Center.Y();
	
	for (unsigned int a = 0; a < tempHits_B.size(); a++)
	{
	if(tempHits_B[a].module==tempHits_B[f_test].module && (tempHits_B[a].channel-flipY)==tempHits_B[f_test].channel) hitYPrelim=hit2Center.Y()+1.25;
	}
	double hitY=hitYPrelim;

	
        double hitZ = (hit1Center.Z() + hit2Center.Z()) / 2.f;

        recoHits rHits;
        rHits.adcX=tempHits_B[f].adc;
	rHits.adcY=tempHits_B[f_test].adc;
        rHits.hitPositionX = hitX;
        rHits.hitPositionY = hitY;
        rHits.hitPositionZ = hitZ;
	rHits.geoX=tempHits_B[f].module;
	rHits.geoY=tempHits_B[f_test].module;
	rHits.stripX=tempHits_B[f].channel;
	rHits.stripY=flipChannel;
	rHits.timeAvg = (tempHits_B[f_test].triggerTime+tempHits_B[f].triggerTime)/2.0;
       if (fabs(tempHits_B[f_test].triggerTime-tempHits_B[f].triggerTime)<fModuletoModuleTimingCut) primaryHits_B.push_back(rHits); 
     //primaryHits_B.push_back(rHits);
	 }
    }
  }
        int pixel0 = -1;
        int pixel1 = -1;
	if (!fMCCSwitch)
	{
       const auto& pdspctbs = *event.getValidHandle<std::vector<raw::ctb::pdspctb>>(fCTBLabel);
       std::vector<int> uS, dS;
	const size_t npdspctbs = pdspctbs.size();
	for(size_t j=0;j<npdspctbs;++j)
	  {
	    const std::vector<raw::ctb::Trigger> HLTriggers = pdspctbs[j].GetHLTriggers();
	    const std::vector<raw::ctb::ChStatus> chs = pdspctbs[j].GetChStatusAfterHLTs();
for (size_t k=0; k<HLTriggers.size(); ++k)
	      { 
		//cout<<chs[k].timestamp<<endl;
		int num = chs[k].crt;
		//cout<<num<<endl;
	
        const std::string binary = std::bitset<32>(num).to_string();	
	const auto crtmask=chs[k].crt;
         pixel0 = -1;
         pixel1 = -1;
	//cout<<crtmask<<endl;
        for (int i = 0; i<32; ++i){
          if (crtmask & (1<<i)){
            if (i<16){
              pixel0 = i;
            }
            else {
              pixel1 = i;
            }
          }
   	}
        if (pixel0!=-1 && pixel1!=-1) {
	cout<<nEvents<<" TJYang Pixels: "<<pixel0<<","<<pixel1<<endl;
	}
	else if (fCTBTriggerOnly) return;
	      }
	  }
	}
// Make tracks from all front and back CRT hits
for (unsigned int f = 0; f < primaryHits_F.size(); f++) {
    for (unsigned int b = 0; b < primaryHits_B.size(); b++) {
   //if (fabs(primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg)<1000) cout<<primaryHits_F[f].geoX<<','<<primaryHits_B[b].geoX<<','<<primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg<<endl;
   if (!fMCCSwitch && fCTBTriggerOnly){
   if (fabs(primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg)<fFronttoBackTimingCut && moduletoCTB(primaryHits_F[f].geoX, primaryHits_F[f].geoY)==pixel0 && moduletoCTB(primaryHits_B[b].geoX, primaryHits_B[b].geoY)==pixel1 ){
	//cout<<"Reconstructed Hits Converted to CTB: "<<moduletoCTB(primaryHits_F[f].geoX, primaryHits_F[f].geoY)<<','<<moduletoCTB(primaryHits_B[b].geoX, primaryHits_B[b].geoY)<<endl;
      combHits tempHits;
      tempHits.hitPositionX1 = primaryHits_F[f].hitPositionX;
      tempHits.hitPositionY1 = primaryHits_F[f].hitPositionY;
      tempHits.hitPositionZ1 = primaryHits_F[f].hitPositionZ;
      tempHits.hitPositionX2 = primaryHits_B[b].hitPositionX;
      tempHits.hitPositionY2 = primaryHits_B[b].hitPositionY;
      tempHits.hitPositionZ2 = primaryHits_B[b].hitPositionZ;
      tempHits.timeDiff=primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg;
      tempHits.t0=primaryHits_F[f].timeAvg;
      tempHits.adcX2=primaryHits_B[b].adcX;
      tempHits.adcX1=primaryHits_F[f].adcX;
      tempHits.adcY2=primaryHits_B[b].adcY;
      tempHits.adcY1=primaryHits_F[f].adcY;
        tempHits.moduleX1 = primaryHits_F[f].geoX;
        tempHits.moduleX2 = primaryHits_B[b].geoX;
        tempHits.moduleY1 = primaryHits_F[f].geoY;
        tempHits.moduleY2 = primaryHits_B[b].geoY;
        tempHits.stripX1 = primaryHits_F[f].stripX;
        tempHits.stripX2 = primaryHits_B[b].stripX;
        tempHits.stripY1 = primaryHits_F[f].stripY;
        tempHits.stripY2 = primaryHits_B[b].stripY;
      combTrackHits.push_back(tempHits);
}
    }
	else {    
	if (fabs(primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg)<fFronttoBackTimingCut){
	//cout<<"Reconstructed Hits Converted to CTB: "<<moduletoCTB(primaryHits_F[f].geoX, primaryHits_F[f].geoY)<<','<<moduletoCTB(primaryHits_B[b].geoX, primaryHits_B[b].geoY)<<endl;
      combHits tempHits;
      tempHits.hitPositionX1 = primaryHits_F[f].hitPositionX;
      tempHits.hitPositionY1 = primaryHits_F[f].hitPositionY;
      tempHits.hitPositionZ1 = primaryHits_F[f].hitPositionZ;
      tempHits.hitPositionX2 = primaryHits_B[b].hitPositionX;
      tempHits.hitPositionY2 = primaryHits_B[b].hitPositionY;
      tempHits.hitPositionZ2 = primaryHits_B[b].hitPositionZ;
      tempHits.timeDiff=primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg;
      tempHits.t0=primaryHits_F[f].timeAvg;
      tempHits.adcX2=primaryHits_B[b].adcX;
      tempHits.adcX1=primaryHits_F[f].adcX;
      tempHits.adcY2=primaryHits_B[b].adcY;
      tempHits.adcY1=primaryHits_F[f].adcY;
        tempHits.moduleX1 = primaryHits_F[f].geoX;
        tempHits.moduleX2 = primaryHits_B[b].geoX;
        tempHits.moduleY1 = primaryHits_F[f].geoY;
        tempHits.moduleY2 = primaryHits_B[b].geoY;
        tempHits.stripX1 = primaryHits_F[f].stripX;
        tempHits.stripX2 = primaryHits_B[b].stripX;
        tempHits.stripY1 = primaryHits_F[f].stripY;
        tempHits.stripY2 = primaryHits_B[b].stripY;
      combTrackHits.push_back(tempHits);
}
}
}
  }
 
  // Reconstruciton information
  art::Handle < vector < recob::Track > > trackListHandle;
  vector < art::Ptr < recob::Track > > trackList;
  art::Handle< std::vector<recob::PFParticle> > PFPListHandle; 
  vector<art::Ptr<recob::PFParticle> > pfplist;


  if (event.getByLabel(fTrackModuleLabel, trackListHandle)) {
    art::fill_ptr_vector(trackList, trackListHandle);
  }


  if(event.getByLabel("pandora",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);

  art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, event ,"pandora");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,event,"pandoraTrack");
  art::FindManyP<recob::Hit> trackHits(trackListHandle, event, "pandoraTrack");
  int nTracksReco = trackList.size();
  //cout<<"Number of Potential CRT Reconstructed Through-Going Muon: "<<combTrackHits.size()<<endl;
    art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits
    std::vector<art::Ptr<recob::Hit>> hitlist;
    if(event.getByLabel("hitpdune", hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
  art::FindManyP < recob::Hit > hitsFromTrack(trackListHandle, event, fTrackModuleLabel);
  int tempId = 0;
  allTracksPair.clear();
  for (int iRecoTrack = 0; iRecoTrack < nTracksReco; ++iRecoTrack) {
    if (combTrackHits.size()<1) break;
    std::vector< art::Ptr<recob::Hit> > allHits =  hitsFromTrack.at(iRecoTrack);

      art::Ptr<recob::Track> ptrack(trackListHandle, iRecoTrack);

     
	std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(iRecoTrack);
	if(!pfps.size()) continue;
	std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
	if(t0s.size()){ 
	  auto t0=t0s.at(0);
	  int t_zero=t0->Time();
	  cout<<"Pandora T0: "<<t_zero<<endl;
   	}
    



    int nTrajectoryPoints = trackList[iRecoTrack] -> NumberTrajectoryPoints();
    int lastPoint = nTrajectoryPoints;
// Get track positions and find angles

    double trackStartPositionZ_noSCE = trackList[iRecoTrack]->Vertex().Z();
    double trackEndPositionZ_noSCE = trackList[iRecoTrack] -> End().Z();

    double trackStartPositionX_notCorrected = trackList[iRecoTrack]->Vertex().X();
    double trackStartPositionY_noSCE = trackList[iRecoTrack]->Vertex().Y();


    double trackEndPositionX_notCorrected = trackList[iRecoTrack] -> End().X();
    double trackEndPositionY_noSCE = trackList[iRecoTrack] -> End().Y();

    int firstHit=0;
    int lastHit=allHits.size()-2;
    if (trackStartPositionZ_noSCE>trackEndPositionZ_noSCE){
    trackEndPositionZ_noSCE = trackList[iRecoTrack]->Vertex().Z();
    trackStartPositionZ_noSCE = trackList[iRecoTrack] -> End().Z();
    trackEndPositionX_notCorrected = trackList[iRecoTrack]->Vertex().X();
    trackEndPositionY_noSCE = trackList[iRecoTrack]->Vertex().Y();


    trackStartPositionX_notCorrected=trackList[iRecoTrack] -> End().X();
    trackStartPositionY_noSCE = trackList[iRecoTrack] -> End().Y();
    firstHit=lastHit;
    lastHit=0;
    }




   



    if ((trackEndPositionZ_noSCE > 660 && trackStartPositionZ_noSCE < 50) || (trackStartPositionZ_noSCE > 660 && trackEndPositionZ_noSCE < 50)) {

/*
if(fMCCSwitch){

      art::ServiceHandle < cheat::BackTrackerService > backTracker;
      art::ServiceHandle < cheat::ParticleInventoryService > partInventory;
      art::Ptr<recob::Track> ptrack(trackListHandle, iRecoTrack);

      std::vector<art::Ptr<recob::Hit>> allHits=trackHits.at(iRecoTrack); //storing hits for ith track
 
      int trackid=-1;
      std::map<int,double> trkide;
      std::map<int,double> trknumelec;


      for(size_t h=0; h<allHits.size();h++){
	art::Ptr<recob::Hit> hit=allHits[h];
	std::vector<sim::TrackIDE> eveIDs = backTracker->HitToTrackIDEs(hit);
	for(size_t e=0;e<eveIDs.size(); ++e){
	  trkide[eveIDs[e].trackID] += eveIDs[e].energy;
	}
      }
      double  maxe = -1;
      double tote = 0;
      for(std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
	tote += ii->second;
	if((ii->second)>maxe){
	  maxe = ii->second;
	  trackid = ii->first;
	 
	}
      }
	const simb::MCParticle *particle = partInventory->TrackIdToParticle_P(trackid);

	
	int nTrajectory=particle->NumberTrajectoryPoints();

	int approxExit=0;
	// Find Exit
	for (int i=0; i<nTrajectory-2; i++){
	approxExit=i;
	if (particle->Position(i).Z()>1077) break;
	}
	int beamLeft=0;
	for (int i=0; i<nTrajectory-2; i++){
	beamLeft=i;
	if (particle->Position(i).Z()>-254) break;
	}

	int beamRight=0;
	for (int i=0; i<nTrajectory-2; i++){
	beamRight=i;
	if (particle->Position(i).Z()>-1000) break;
	}
	cout<<beamRight<<beamLeft<<approxExit<<endl;
	if (((particle->Position(beamLeft).Y()<620 && particle->Position(beamLeft).Y()>-43 && (particle->Position(beamLeft).X())<133 && (particle->Position(beamLeft).X())>-196) || (particle->Position(beamRight).Y()<620 && particle->Position(beamRight).Y()>-43 && (particle->Position(beamRight).X())<401  && (particle->Position(beamRight).X())>229) )  && (particle->Position(beamRight).Y()<548 && particle->Position(beamRight).Y()>-147 && (particle->Position(beamRight).X())<336  && (particle->Position(beamRight).X())>-333) ){
truthEnergy=particle->E();
mccTrackId=trackid;
fMCCMuon->Fill();

}

}*/

      
      for (unsigned int iCombinatorialTrack = 0; iCombinatorialTrack < combTrackHits.size(); iCombinatorialTrack++) {
		double xOffset=0;
    

		double trackStartPositionX_noSCE=trackStartPositionX_notCorrected;
		double trackEndPositionX_noSCE=trackEndPositionX_notCorrected;
		if (t0s.empty()){
		int RDOffset=0;
		if (!fMCCSwitch) RDOffset=111;
		double ticksOffset=0;
		cout<<(combTrackHits[iCombinatorialTrack].t0+RDOffset)<<endl;
		cout<<detectorPropertiesService->GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat)<<endl;
		if (!fMCCSwitch) ticksOffset = (combTrackHits[iCombinatorialTrack].t0+RDOffset)/25.f+detectorPropertiesService->GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);

		else if (fMCCSwitch) ticksOffset = (combTrackHits[iCombinatorialTrack].t0/500.f)+detectorPropertiesService->GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		
	       xOffset=detectorPropertiesService->ConvertTicksToX(ticksOffset,allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		//double xOffset=.08*ticksOffset
		
	 trackStartPositionX_noSCE=trackStartPositionX_notCorrected-xOffset;
         trackEndPositionX_noSCE=trackEndPositionX_notCorrected-xOffset;
	}

   double trackStartPositionX=trackStartPositionX_noSCE;
   double trackStartPositionY=trackStartPositionY_noSCE;
   double trackStartPositionZ=trackStartPositionZ_noSCE;

   double trackEndPositionX=trackEndPositionX_noSCE;
   double trackEndPositionY=trackEndPositionY_noSCE;
   double trackEndPositionZ=trackEndPositionZ_noSCE;


    cout<<fSCECorrection<<endl;
    if (fSCECorrection){
     trackStartPositionX=trackStartPositionX_noSCE-SCE->GetPosOffsets(geo::Point_t(trackStartPositionX_noSCE, trackStartPositionY_noSCE, trackStartPositionZ_noSCE)).X();
     trackStartPositionY=trackStartPositionY_noSCE+SCE->GetPosOffsets(geo::Point_t(trackStartPositionX_noSCE, trackStartPositionY_noSCE, trackStartPositionZ_noSCE)).Y();
     trackStartPositionZ=trackStartPositionZ_noSCE+SCE->GetPosOffsets(geo::Point_t(trackStartPositionX_noSCE, trackStartPositionY_noSCE, trackStartPositionZ_noSCE)).Z();


     trackEndPositionX=trackEndPositionX_noSCE-SCE->GetPosOffsets(geo::Point_t(trackEndPositionX_noSCE, trackEndPositionY_noSCE, trackEndPositionZ_noSCE)).X();
     trackEndPositionY=trackEndPositionY_noSCE+SCE->GetPosOffsets(geo::Point_t(trackEndPositionX_noSCE, trackEndPositionY_noSCE, trackEndPositionZ_noSCE)).Y();
     trackEndPositionZ=trackEndPositionZ_noSCE+SCE->GetPosOffsets(geo::Point_t(trackEndPositionX_noSCE, trackEndPositionY_noSCE, trackEndPositionZ_noSCE)).Z();
	}
        double X1 = combTrackHits[iCombinatorialTrack].hitPositionX1;
        double X2 = combTrackHits[iCombinatorialTrack].hitPositionX2;
        double Y1 = combTrackHits[iCombinatorialTrack].hitPositionY1;
        double Y2 = combTrackHits[iCombinatorialTrack].hitPositionY2;
        double Z1 = combTrackHits[iCombinatorialTrack].hitPositionZ1;
        double Z2 = combTrackHits[iCombinatorialTrack].hitPositionZ2;
	double timeCRTDelta= combTrackHits[iCombinatorialTrack].timeDiff;



        double averageSignedDistance = 0;
        double averageSignedDistanceYZ = 0;
        double averageSignedDistanceXZ = 0;
        double averageSignedDistanceXY = 0;

        for (int trackpoint = 0; trackpoint < lastPoint; trackpoint++) {
	  double trackPosX=trackList[iRecoTrack] -> LocationAtPoint(trackpoint).X()+xOffset;
	  double trackPosY=trackList[iRecoTrack] -> LocationAtPoint(trackpoint).Y();
	  double trackPosZ=trackList[iRecoTrack] -> LocationAtPoint(trackpoint).Z();
	   TVector3 trackPos(trackPosX, trackPosY, trackPosZ);
			double distanceYZ = signedPointToLineDistance( Y1,Z1, Y2,Z2, trackPos.Y(), trackPos.Z() ); //only the Y and Z of trackpos will be used
			double distanceXZ = signedPointToLineDistance( X1,Z1, X2,Z2, trackPos.X(), trackPos.Z() );

double distanceXY = signedPointToLineDistance( X1,Y1, X2,Y2, trackPos.X(), trackPos.Y());
			double distance=signed3dDistance(X1, Y1, Z1, X2, Y2, Z2, trackPos);
			
                                                                                      
			
			//calculate average distance (signed)
			averageSignedDistance += distance/(lastPoint+1);	
		        averageSignedDistanceYZ += distanceYZ/(lastPoint+1);
			averageSignedDistanceXZ += distanceXZ/(lastPoint+1);
averageSignedDistanceXY += distanceXY/(lastPoint+1);
	
		    }

        double crtAngleXZ = setAngle(atan2(X2 - X1, Z2 - Z1));
        double crtAngleYZ = setAngle(atan2(Y2 - Y1, Z2 - Z1));
        double crtAngleXY = setAngle(atan2(X2 - X1, Y2 - Y1));
        double recoAngleYZ = setAngle(atan2(trackEndPositionY - trackStartPositionY, trackEndPositionZ - trackStartPositionZ));
        double recoAngleXZ = setAngle(atan2(trackEndPositionX - trackStartPositionX, trackEndPositionZ - trackStartPositionZ));

        double recoAngleXY = setAngle(atan2(trackEndPositionX - trackStartPositionX, trackEndPositionY - trackStartPositionY));
	// Make metrics for a CRT pair to compare later
	TVector3 trackStart(trackStartPositionX, trackStartPositionY, trackStartPositionZ);
	TVector3 trackEnd(trackEndPositionX, trackEndPositionY, trackEndPositionZ);
	TVector3 v1(X1,Y1,Z1);
	TVector3 v2(X2, Y2, Z2);


            TVector3 v4(trackStartPositionX,
                        trackStartPositionY,
                        trackStartPositionZ);
            TVector3 v5(trackEndPositionX,
                        trackEndPositionY,
                        trackEndPositionZ);
	TVector3 trackVector = (v5-v4).Unit();
	TVector3 hitVector=(v2-v1).Unit();

              double predictedHitPositionY1 = (v1.Z()-v5.Z())/(v4.Z()-v5.Z())*(v4.Y()-v5.Y())+v5.Y();
              double predictedHitPositionY2 = (v2.Z()-v5.Z())/(v4.Z()-v5.Z())*(v4.Y()-v5.Y())+v5.Y();

              double predictedHitPositionX1 = (v1.Z()-v5.Z())/(v4.Z()-v5.Z())*(v4.X()-v5.X())+v5.X();
              double predictedHitPositionX2 = (v2.Z()-v5.Z())/(v4.Z()-v5.Z())*(v4.X()-v5.X())+v5.X();

	double dotProductCos=trackVector*hitVector;

        double deltaX1 = (predictedHitPositionX1-X1);
        double deltaX2 = (predictedHitPositionX2-X2);

	double deltaX=fabs(deltaX2)+fabs(deltaX1);

        double deltaY1 = (predictedHitPositionY1-Y1);

        double deltaY2 = (predictedHitPositionY2-Y2);
	double deltaY=fabs(deltaY2)+fabs(deltaY1);
        tracksPair tPair;
        tPair.tempId = tempId;
        tPair.CRTTrackId = iCombinatorialTrack;
        tPair.recoId = iRecoTrack;
        tPair.deltaX_F = deltaX1;
 
        tPair.deltaX_B = deltaX2;
        tPair.deltaY_F = deltaY1;
        tPair.deltaY_B = deltaY2;
	tPair.timeDiff=timeCRTDelta;
	tPair.deltaX=deltaX;
        tPair.deltaY=deltaY;
        tPair.dotProductCos=dotProductCos;
        tPair.deltaAngleXY = (recoAngleXY - crtAngleXY);
        tPair.deltaAngleYZ = (recoAngleYZ - crtAngleYZ);
        tPair.deltaAngleXZ = (recoAngleXZ - crtAngleXZ);
        tPair.averageSignedDistanceXY = averageSignedDistanceXY;
        tPair.averageSignedDistanceYZ = averageSignedDistanceYZ;
        tPair.averageSignedDistanceXZ = averageSignedDistanceXZ;
        tPair.averageSignedDistance = averageSignedDistance;
        tPair.moduleX1 = combTrackHits[iCombinatorialTrack].moduleX1;
        tPair.moduleX2 = combTrackHits[iCombinatorialTrack].moduleX2;
        tPair.moduleY1 = combTrackHits[iCombinatorialTrack].moduleY1;
        tPair.moduleY2 = combTrackHits[iCombinatorialTrack].moduleY2;
      tPair.adcX2=combTrackHits[iCombinatorialTrack].adcX2;
      tPair.adcX1=combTrackHits[iCombinatorialTrack].adcX1;
      tPair.adcY2=combTrackHits[iCombinatorialTrack].adcY2;
      tPair.adcY1=combTrackHits[iCombinatorialTrack].adcY1;

        tPair.stripX1 = combTrackHits[iCombinatorialTrack].stripX1;
        tPair.stripX2 = combTrackHits[iCombinatorialTrack].stripX2;
        tPair.stripY1 = combTrackHits[iCombinatorialTrack].stripY1;
        tPair.stripY2 = combTrackHits[iCombinatorialTrack].stripY2;
        tPair.X1 = X1;
        tPair.Y1 = Y1;
        tPair.Z1 = Z1;
        tPair.X2 = X2;
        tPair.Y2 = Y2;
        tPair.Z2 = Z2;
	tPair.xOffset=xOffset;
	tPair.t0=combTrackHits[iCombinatorialTrack].t0;
        tPair.trackStartPosition=trackStart;
	tPair.trackEndPosition=trackEnd;
	if (t0s.empty()) tPair.pandoraT0Check=0;
	else tPair.pandoraT0Check=1;
        allTracksPair.push_back(tPair);


      }

      tempId++;
    } //iRecoTrack
    }

     //Sort pair by ascending order of absolute distance
    sort(allTracksPair.begin(), allTracksPair.end(), sortPair());

    // Compare, sort, and eliminate CRT hits for just the best one


    vector < tracksPair > allUniqueTracksPair;
    while (allTracksPair.size()) {
      allUniqueTracksPair.push_back(allTracksPair.front());
      allTracksPair.erase(remove_if(allTracksPair.begin(), allTracksPair.end(), removePairIndex(allTracksPair.front())),
        allTracksPair.end());
    }

	cout<<"Number of reco and CRT pairs: "<<allUniqueTracksPair.size()<<endl;
// For the best one, add the validation metrics to a tree
    if (allUniqueTracksPair.size() > 0) {
      for (unsigned int u = 0; u < allUniqueTracksPair.size(); u++) {
	  trackID = allUniqueTracksPair[u].recoId;
         averageSignedDistanceYZ = allUniqueTracksPair[u].averageSignedDistanceYZ;

        averageSignedDistanceXZ = allUniqueTracksPair[u].averageSignedDistanceXZ;
        averageSignedDistance = allUniqueTracksPair[u].averageSignedDistance;

	displAngleXY = allUniqueTracksPair[u].deltaAngleXY;

	displAngleYZ = allUniqueTracksPair[u].deltaAngleYZ;

	displAngleXZ = allUniqueTracksPair[u].deltaAngleXZ;

	deltaX_F=allUniqueTracksPair[u].deltaX_F;
	deltaX_B=allUniqueTracksPair[u].deltaX_B;
	deltaY_F=allUniqueTracksPair[u].deltaY_F;
	deltaY_B=allUniqueTracksPair[u].deltaY_B;
	dotCos=allUniqueTracksPair[u].dotProductCos;
	trackX1=allUniqueTracksPair[u].trackStartPosition.X();
	trackY1=allUniqueTracksPair[u].trackStartPosition.Y();
	trackZ1=allUniqueTracksPair[u].trackStartPosition.Z();

	trackX2=allUniqueTracksPair[u].trackEndPosition.X();
	trackY2=allUniqueTracksPair[u].trackEndPosition.Y();
	trackZ2=allUniqueTracksPair[u].trackEndPosition.Z();

	moduleX_F=allUniqueTracksPair[u].moduleX1;
	moduleX_B=allUniqueTracksPair[u].moduleX2;
	moduleY_F=allUniqueTracksPair[u].moduleY1;
	moduleY_B=allUniqueTracksPair[u].moduleY2;

	adcX_F=allUniqueTracksPair[u].adcX1;
	adcY_F=allUniqueTracksPair[u].adcY1;
	adcX_B=allUniqueTracksPair[u].adcX2;
	adcY_B=allUniqueTracksPair[u].adcY2;

	stripX_F=allUniqueTracksPair[u].stripX1;
	stripX_B=allUniqueTracksPair[u].stripX2;
	stripY_F=allUniqueTracksPair[u].stripY1;
	stripY_B=allUniqueTracksPair[u].stripY2;

	X_F=allUniqueTracksPair[u].X1;
	Y_F=allUniqueTracksPair[u].Y1;
	Z_F=allUniqueTracksPair[u].Z1;

	X_B=allUniqueTracksPair[u].X2;
	Y_B=allUniqueTracksPair[u].Y2;
	Z_B=allUniqueTracksPair[u].Z2;
	recoPandoraT0Check=allUniqueTracksPair[u].pandoraT0Check;

	measuredT0=allUniqueTracksPair[u].t0;
	measuredXOffset=allUniqueTracksPair[u].xOffset;
        cout<<fabs(allUniqueTracksPair[u].dotProductCos)<<endl;
	CRT_TOF=allUniqueTracksPair[u].timeDiff;	
        if (fabs(allUniqueTracksPair[u].dotProductCos)>0.99) {
        //cout<<allUniqueTracksPair[u].timeDiff<<endl;
	//cout<<fabs(allUniqueTracksPair[u].dotProductCos)<<endl;

          cout << "Displacement: " << averageSignedDistance << ',' << averageSignedDistanceYZ << ',' << averageSignedDistanceXZ << endl;
	  cout<< "Delta X and Y: "<<deltaX_F<<','<<deltaY_F<<','<<deltaX_B<<','<<deltaY_B<<endl;
	  cout<< "Predicted X and Y: "<< deltaX_F+X_F<<','<<deltaY_F+Y_F<<','<<deltaX_B+X_B<<','<<deltaY_B+Y_B<<endl;
	  cout<< "Detected X and Y: "<< X_F<<','<<Y_F<<','<<X_B<<','<<Y_B<<endl;
	  cout<<moduleX_F<<','<<stripX_F<<endl;
	  cout<<moduleY_F<<','<<stripY_F<<endl;
	  cout<<moduleX_B<<','<<stripX_B<<endl;
	  cout<<moduleY_B<<','<<stripY_B<<endl;
	  cout<<"ADC Values: "<<adcX_F<<','<<adcY_F<<','<<adcX_B<<','<<adcY_B<<endl;
	fCRTTree->Fill();
	   

        }
      }
    }

 }


// Setup CRT 
void CRT::TwoCRTMatching::beginJob() {
	art::ServiceHandle<art::TFileService> fileServiceHandle;
       fCRTTree = fileServiceHandle->make<TTree>("Displacement", "event by event info");
        fMCCMuon= fileServiceHandle->make<TTree>("MCCTruths", "event by event info");
	fCRTTree->Branch("nEvents", &nEvents, "fnEvents/I");
	fCRTTree->Branch("hRun", &run, "run/I");
	fCRTTree->Branch("hSubRun", &subRun, "subRun/I");
	fCRTTree->Branch("htrackID", &trackID, "trackID/I");
        fCRTTree->Branch("hmccTrackId", &mccTrackId, "mccTrackId/I");
	fMCCMuon->Branch("nEvents", &nEvents, "nEvents/I");
	fMCCMuon->Branch("truthEnergy", &truthEnergy, "truthEnergy/D");
	fMCCMuon->Branch("hRun", &run, "run/I");
	fMCCMuon->Branch("hSubRun", &subRun, "subRun/I");
        fMCCMuon->Branch("hmccTrackId", &mccTrackId, "mccTrackId/I");
	fCRTTree->Branch("nHitsPerEvent", &nHitsPerEvent, "fnHitsPerEvent/I");
	fCRTTree->Branch("haverageSignedDistanceYZ", &averageSignedDistanceYZ, "averageSignedDistanceYZ/D");
	fCRTTree->Branch("haverageSignedDistanceXZ", &averageSignedDistanceXZ, "averageSignedDistanceXZ/D");
	fCRTTree->Branch("haverageSignedDistance", &averageSignedDistance, "averageSignedDistance/D");
	fCRTTree->Branch("hdisplAngleXY", &displAngleXY, "displAngleXY/D");
	fCRTTree->Branch("hdisplAngleYZ", &displAngleYZ, "displAngleYZ/D");
	fCRTTree->Branch("hdisplAngleXZ", &displAngleXZ, "displAngleXZ/D");

	fCRTTree->Branch("hdeltaX_F", &deltaX_F, "deltaX_F/D");
	fCRTTree->Branch("hdeltaX_B", &deltaX_B, "deltaX_B/D");
	fCRTTree->Branch("hdeltaY_F", &deltaY_F, "deltaY_F/D");
	fCRTTree->Branch("hdeltaY_B", &deltaY_B, "deltaY_B/D");
	fCRTTree->Branch("hdotProductCos", &dotCos, "dotCos/D");


	fCRTTree->Branch("hX_F", &X_F, "X_F/D");
	fCRTTree->Branch("hX_B", &X_B, "X_B/D");
	fCRTTree->Branch("hY_F", &Y_F, "Y_F/D");
	fCRTTree->Branch("hY_B", &Y_B, "Y_B/D");
	fCRTTree->Branch("hZ_F", &Z_F, "Z_F/D");
	fCRTTree->Branch("hZ_B", &Z_B, "Z_B/D");

	fCRTTree->Branch("htrackStartX", &trackX1, "trackX1/D");
	fCRTTree->Branch("htrackStartY", &trackY1, "trackY1/D");
	fCRTTree->Branch("htrackStartZ", &trackZ1, "trackZ1/D");

	fCRTTree->Branch("htrackEndX", &trackX2, "trackX2/D");
	fCRTTree->Branch("htrackEndY", &trackY2, "trackY2/D");
	fCRTTree->Branch("htrackEndZ", &trackZ2, "trackZ2/D");
	fCRTTree->Branch("CRT_TOF", &CRT_TOF, "CRT_TOF/D");

	fCRTTree->Branch("hmoduleX_F", &moduleX_F, "moduleX_F/I");
	fCRTTree->Branch("hmoduleX_B", &moduleX_B, "moduleX_B/I");
	fCRTTree->Branch("hmoduleY_F", &moduleY_F, "moduleY_F/I");
	fCRTTree->Branch("hmoduleY_B", &moduleY_B, "moduleY_B/I");

	fCRTTree->Branch("hstripX_F", &stripX_F, "stripX_F/I");
	fCRTTree->Branch("hstripX_B", &stripX_B, "stripX_B/I");
	fCRTTree->Branch("hstripY_F", &moduleY_F, "stripY_F/I");
	fCRTTree->Branch("hstripY_B", &stripY_B, "stripY_B/I");

	fCRTTree->Branch("hadcX_F", &adcX_F, "adcX_F/I");
	fCRTTree->Branch("hadcY_F", &adcY_F, "adcY_F/I");
	fCRTTree->Branch("hadcX_B", &adcX_B, "adcX_B/I");
	fCRTTree->Branch("hadcY_B", &adcY_B, "adcY_B/I");

	fCRTTree->Branch("hmeasuredT0", &measuredT0, "measuredT0/D");
	fCRTTree->Branch("hmeasuredXOffset", &measuredXOffset, "measuredXOffset/D");
	fCRTTree->Branch("hrecoPandoraT0Check", &recoPandoraT0Check, "recoPandoraT0Check/B");


}
// Endjob actions
void CRT::TwoCRTMatching::endJob() 
{



}














DEFINE_ART_MODULE(CRT::TwoCRTMatching)
