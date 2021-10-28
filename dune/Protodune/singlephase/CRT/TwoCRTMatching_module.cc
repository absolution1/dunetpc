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
#include "nug4/ParticleNavigation/ParticleList.h"


#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"


#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"


//Local includes
#include "dune/Protodune/singlephase/CRT/data/CRTTrigger.h"



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
#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
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
  bool moduleMatcher(int module1, int module2);
  double signedPointToLineDistance(double firstPoint1, double firstPoint2,   double secondPoint1, double secondPoint2, double trackPoint1, double trackPoint2);
  double signed3dDistance(double firstPoint1, double firstPoint2, double firstPoint3, double secondPoint1, double secondPoint2, double secondPoint3, TVector3 trackPos);
  void beginJob() override;
  void endJob() override;
  void createPNG(TH1D * histo);
  double setAngle(double angle);
  int nEvents = 0;
  int nHaloMuons = 0;
  int nHitsPerEvent=0;
  int nMCCMuons=0;
  private: ofstream logFile;

  //Parameters for reading in CRT::Triggers and associated AuxDetSimChannels.
  art::InputTag fCRTLabel; //Label for the module that analyzed 
  art::InputTag fCTBLabel;
  TTree * fCRTTree;
  //TTree * fTimingTree;
  TTree * fCRTdQTree;
  TTree * fMCCMuon;
  TTree * fTrackInfo;
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
    double mccT0;
    double SCECorrectX_F, SCECorrectY_F, SCECorrectZ_F;
    double SCECorrectX_B, SCECorrectY_B, SCECorrectZ_B;
    int endTPC;
    double trkPosX, trkPosY, trkPosZ, displXZ, displYZ;
    double mccCRTStartX, mccCRTStartY, mccCRTStartZ;
    double mccCRTEndX, mccCRTEndY, mccCRTEndZ;
    int candidateCRT;
    double trkhitx, trkhity,trkhitz, trkhitt0, crtt0, trkdqdx;
    int trkhitIntegral;
    //int tmoduleX, tDiff, tmoduleY;
    //double xPos, yPos, zPos;
    double  sigmaHit;
    int TPCID, WireID;
     int sumADC, rangeTime; 
   long long timeStamp;
  typedef struct // Structures for arrays to move hits from raw to reco to validation
  {

    int channel;
    int module;
    int channelGeo;
    int adc;
    int triggerTime;
    int triggerNumber;
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
    int trigNumberX;
    int trigNumberY;
  }
  recoHits;



  typedef struct // These are displacement metrics for track and hit reco
  {
    int tempId;
    int recoId;
    int adcX1;
    int adcY1;
    int adcX2;
    int adcY2;
    double deltaX_F;
    double deltaX_B;
    double deltaY_F;
    double deltaY_B;

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
    double mccT0;
    int endTPC;
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

	//return (fabs(pair1.dotProductCos)>fabs(pair2.dotProductCos));
        return (fabs(pair1.deltaX_F)+fabs(pair1.deltaY_F)+fabs(pair1.deltaX_B)+fabs(pair1.deltaY_B)<fabs(pair2.deltaX_F)+fabs(pair2.deltaY_F)+fabs(pair2.deltaX_B)+fabs(pair2.deltaY_B));
	//return ((fabs(pair1.dotProductCos)>.998 && pair1.deltaY<pair2.deltaY && pair1.deltaX<pair2.deltaX));

  }
  };

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
bool CRT::TwoCRTMatching::moduleMatcher(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 6 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 26 || module2 == 27)) || (module1 == 31 && (module2 == 26 || module2 == 27)) || (module1 == 7 && (module2 == 12 || module2 == 13)) || (module1 == 15 && (module2 == 12 || module2 == 13)) || (module1 == 18 && (module2 == 24 || module2 == 25)) || (module1 == 30 && (module2 == 24 || module2 == 25)) || (module1 == 1 && (module2 == 4 || module2 == 5)) || (module1 == 9 && (module2 == 4 || module2 == 5)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 0 && (module2 == 2 || module2 == 3)) || (module1 == 8 && (module2 == 2 || module2 == 3)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23))) return 1;
  else return 0;

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
    fADCThreshold=800;
    fModuletoModuleTimingCut=40;
    fFronttoBackTimingCut=120;
    
}
  else {
    fModuleSwitch=0;
    fADCThreshold=20;
    fModuletoModuleTimingCut=5;
    fFronttoBackTimingCut=8;
    art::ValidHandle<std::vector<raw::RDTimeStamp>> timingHandle = event.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");
    timeStamp=timingHandle->at(0).GetTimeStamp();

} 


     auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();


   if(!fMCCSwitch){
   //const auto& pdspctbs = event.getValidHandle<std::vector<raw::ctb::pdspctb>>(fCTB_tag);
	art::ValidHandle<std::vector<raw::RDTimeStamp>> timingHandle = event.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");

	const raw::RDTimeStamp& timeStamp = timingHandle->at(0);
	if (fCTBTriggerOnly){
	if(timeStamp.GetFlags()!= 13) return;}
  }
  int nHits = 0;




	//Detector properties service
  primaryHits_F.clear();
  primaryHits_B.clear();
  allTracksPair.clear();
  tempHits_F.clear();
  tempHits_B.clear(); // Arrays to compile hits and move them through


  //Get triggers
  cout << "Getting triggers" << endl;
  const auto & triggers = event.getValidHandle < std::vector < CRT::Trigger >> (fCRTLabel);

  art::FindManyP < sim::AuxDetSimChannel > trigToSim(triggers, event, fCRTLabel);

  //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
  art::ServiceHandle < geo::Geometry > geom;



  //Mapping from channel to trigger
  std::unordered_map < size_t, double > prevTimes;
  int hitID = 0;
  cout << "Looking for hits in Triggers" << endl;

  int trigID=0;
  for (const auto & trigger: * triggers) {
    const auto & hits = trigger.Hits();
    for (const auto & hit: hits) { // Collect hits on all modules
	//cout<<hits.size()<<','<<hit.ADC()<<endl;
      if (hit.ADC() > fADCThreshold) { // Keep if they are above threshold

        tempHits tHits;
	if (!fMCCSwitch){

        tHits.module = trigger.Channel(); // Values to add to array
	tHits.channel=hit.Channel();
        tHits.adc = hit.ADC();
	tHits.triggerTime=trigger.Timestamp()-timeStamp;
	}
	else{
        tHits.module = trigger.Channel(); // Values to add to array
	tHits.channel=hit.Channel();
        tHits.adc = hit.ADC();
	tHits.triggerTime=trigger.Timestamp();
	}
	 //cout<<trigger.Channel()<<','<<hit.Channel()<<','<<hit.ADC()<<endl;
        nHits++;
	tHits.triggerNumber=trigID;
        const auto & trigGeo = geom -> AuxDet(trigger.Channel()); // Get geo  
        const auto & csens = trigGeo.SensitiveVolume(hit.Channel());
        const auto center = csens.GetCenter();
        if (center.Z() < 100) tempHits_F.push_back(tHits); // Sort F/B from Z
        else tempHits_B.push_back(tHits);
        hitID++;
      }
    }
    trigID++;
  }
  nHitsPerEvent=nHits;
  cout << "Hits compiled for event: " << nEvents << endl;
  cout << "Number of Hits above Threshold:  " << hitID << endl;

  for (unsigned int f = 0; f < tempHits_F.size(); f++) {
    for (unsigned int f_test = 0; f_test < tempHits_F.size(); f_test++) {
       if (fabs(tempHits_F[f_test].triggerTime-tempHits_F[f].triggerTime)>fModuletoModuleTimingCut) continue;
      const auto & trigGeo = geom -> AuxDet(tempHits_F[f].module);
      const auto & trigGeo2 = geom -> AuxDet(tempHits_F[f_test].module);

      const auto & hit1Geo = trigGeo.SensitiveVolume(tempHits_F[f].channel);
      const auto hit1Center = hit1Geo.GetCenter();
      // Create 2D hits from geo of the Y and X modules
       const auto & hit2Geo = trigGeo2.SensitiveVolume(tempHits_F[f_test].channel);
      const auto hit2Center = hit2Geo.GetCenter();
      bool moduleMatched;
      moduleMatched=moduleMatcher(tempHits_F[f_test].module, tempHits_F[f].module);
      if (moduleMatched) {
        // Get the center of the hits (CRT_Res=2.5 cm)
        double hitX = hit1Center.X();
	for (unsigned int a = 0; a < tempHits_F.size(); a++)
	{
	if(tempHits_F[a].module==tempHits_F[f].module && (tempHits_F[a].channel-1)==tempHits_F[f].channel) hitX=hit1Center.X()+1.25;
	}
	double hitYPrelim=hit2Center.Y();
	for (unsigned int a = 0; a < tempHits_F.size(); a++)
	{
	if(tempHits_F[a].module==tempHits_F[f_test].module && (tempHits_F[a].channel-1)==tempHits_F[f_test].channel) hitYPrelim=hit2Center.Y()+1.25;
	}
	

	
	double hitY=hitYPrelim;
        double hitZ = (hit1Center.Z() + hit2Center.Z()) / 2.f;

        recoHits rHits;
        rHits.hitPositionX = hitX;
        rHits.hitPositionY = hitY;
        rHits.hitPositionZ = hitZ;
	rHits.geoX=tempHits_F[f].module;
	rHits.geoY=tempHits_F[f_test].module;
	rHits.stripX=tempHits_F[f].channel;
	rHits.stripY=tempHits_F[f_test].channel;

	rHits.adcX=tempHits_F[f].adc;
	rHits.adcY=tempHits_F[f_test].adc;
	rHits.trigNumberX=tempHits_F[f].triggerNumber;
	rHits.trigNumberY=tempHits_F[f_test].triggerNumber;
	rHits.timeAvg = (tempHits_F[f_test].triggerTime+tempHits_F[f].triggerTime)/2.0;
       primaryHits_F.push_back(rHits); // Add array
    }
    }
  }
  for (unsigned int f = 0; f < tempHits_B.size(); f++) {
    for (unsigned int f_test = 0; f_test < tempHits_B.size(); f_test++) { // Same as above but for back CRT
       if (fabs(tempHits_B[f_test].triggerTime-tempHits_B[f].triggerTime)>fModuletoModuleTimingCut) continue;

      const auto & trigGeo = geom -> AuxDet(tempHits_B[f].module);
      const auto & trigGeo2 = geom -> AuxDet(tempHits_B[f_test].module);
      const auto & hit1Geo = trigGeo.SensitiveVolume(tempHits_B[f].channel);
      const auto hit1Center = hit1Geo.GetCenter();

      const auto & hit2Geo = trigGeo2.SensitiveVolume(tempHits_B[f_test].channel);
      const auto hit2Center = hit2Geo.GetCenter();
      bool moduleMatched;
      moduleMatched=moduleMatcher(tempHits_B[f_test].module, tempHits_B[f].module);

      if (moduleMatched) {
        double hitX = hit1Center.X();
	
	
	for (unsigned int a = 0; a < tempHits_B.size(); a++)
	{
	if(tempHits_B[a].module==tempHits_B[f].module && (tempHits_B[a].channel-1)==tempHits_B[f].channel) hitX=hit1Center.X()+1.25;
	}
	
        double hitYPrelim = hit2Center.Y();
	
	for (unsigned int a = 0; a < tempHits_B.size(); a++)
	{
	if(tempHits_B[a].module==tempHits_B[f_test].module && (tempHits_B[a].channel-1)==tempHits_B[f_test].channel) hitYPrelim=hit2Center.Y()+1.25;
	}
	double hitY=hitYPrelim;

	
        double hitZ = (hit1Center.Z() + hit2Center.Z()) / 2.f;

        recoHits rHits;
        rHits.hitPositionX = hitX;
        rHits.hitPositionY = hitY;
        rHits.hitPositionZ = hitZ;
	rHits.geoX=tempHits_B[f].module;
	rHits.geoY=tempHits_B[f_test].module;
	rHits.adcX=tempHits_B[f].adc;
	rHits.adcY=tempHits_B[f_test].adc;
	rHits.stripX=tempHits_B[f].channel;
	rHits.stripY=tempHits_B[f_test].channel;
	rHits.trigNumberX=tempHits_B[f].triggerNumber;
	rHits.trigNumberY=tempHits_B[f_test].triggerNumber;
	rHits.timeAvg = (tempHits_B[f_test].triggerTime+tempHits_B[f].triggerTime)/2.0;
	 primaryHits_B.push_back(rHits); 

	 }
    }
  }

	std::cout<<primaryHits_F.size()<<','<<primaryHits_B.size()<<std::endl;
  // Reconstruciton information
  
  vector < art::Ptr < recob::Track > > trackList;
  auto trackListHandle = event.getHandle < vector < recob::Track > >(fTrackModuleLabel);
  if (trackListHandle) {
    art::fill_ptr_vector(trackList, trackListHandle);
  }

  vector<art::Ptr<recob::PFParticle> > pfplist;
  auto PFPListHandle = event.getHandle< std::vector<recob::PFParticle> >("pandora");
  if (PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);
     
   if (pfplist.size()<1) return;    
art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, event ,"pandora");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,event,"pandoraTrack");
   
 art::FindManyP<recob::Hit> trackHits(trackListHandle, event, "pandoraTrack");
  int nTracksReco = trackList.size();
  //cout<<"Number of Potential CRT Reconstructed Through-Going Muon: "<<combTrackHits.size()<<endl;
    std::vector<art::Ptr<recob::Hit>> hitlist; // to get information about the hits
    auto hitListHandle = event.getHandle< std::vector<recob::Hit> >("hitpdune");
    if (hitListHandle)
      art::fill_ptr_vector(hitlist, hitListHandle);
  art::FindManyP < recob::Hit > hitsFromTrack(trackListHandle, event, fTrackModuleLabel);


  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);

    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, event, "pandoracalo");
  int tempId = 0;
  allTracksPair.clear();
  for (int iRecoTrack = 0; iRecoTrack < nTracksReco; ++iRecoTrack) {
    if (pfplist.size()<1) break;  
    std::vector< art::Ptr<recob::Hit> > allHits =  hitsFromTrack.at(iRecoTrack);

      art::Ptr<recob::Track> ptrack(trackListHandle, iRecoTrack);

     
	std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(iRecoTrack);
	if(!pfps.size()) continue;
	std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
	/*if(t0s.size()){ 
	  auto t0=t0s.at(0);
	  int t_zero=t0->Time();
	  cout<<"Pandora T0: "<<t_zero<<endl;
   	}*/
    





    double trackStartPositionZ_noSCE = trackList[iRecoTrack]->Vertex().Z();
    double trackEndPositionZ_noSCE = trackList[iRecoTrack] -> End().Z();

    double trackStartPositionX_notCorrected = trackList[iRecoTrack]->Vertex().X();
    double trackStartPositionY_noSCE = trackList[iRecoTrack]->Vertex().Y();


    double trackEndPositionX_notCorrected = trackList[iRecoTrack] -> End().X();
    double trackEndPositionY_noSCE = trackList[iRecoTrack] -> End().Y();

    int firstHit=0;
    //int lastHit=allHits.size()-2;
    if (trackStartPositionZ_noSCE>trackEndPositionZ_noSCE){
    trackEndPositionZ_noSCE = trackList[iRecoTrack]->Vertex().Z();
    trackStartPositionZ_noSCE = trackList[iRecoTrack] -> End().Z();
    trackEndPositionX_notCorrected = trackList[iRecoTrack]->Vertex().X();
    trackEndPositionY_noSCE = trackList[iRecoTrack]->Vertex().Y();


    trackStartPositionX_notCorrected=trackList[iRecoTrack] -> End().X();
    trackStartPositionY_noSCE = trackList[iRecoTrack] -> End().Y();
    
    }




   



    if ((trackEndPositionZ_noSCE > 660 && trackStartPositionZ_noSCE < 50) || (trackStartPositionZ_noSCE > 660 && trackEndPositionZ_noSCE < 50)) {


      double min_delta = DBL_MAX;
      double best_XF = DBL_MAX;
      double best_YF = DBL_MAX;
      double best_ZF = DBL_MAX;
      double best_XB = DBL_MAX;
      double best_YB = DBL_MAX;
      double best_ZB = DBL_MAX;
      double best_dotProductCos = DBL_MAX;
      double best_deltaXF = DBL_MAX;
      double best_deltaYF = DBL_MAX;
      double best_deltaXB = DBL_MAX;
      double best_deltaYB = DBL_MAX;
      double best_T=DBL_MAX;
      int bestHitIndex_F=0;
      int bestHitIndex_B=0;

      /*
      int best_trigXF=0;
      int best_trigYF=0;
      int best_trigXB=0;
      int best_trigYB=0;
       */
    for (unsigned int f = 0; f < primaryHits_F.size(); f++) {
        for (unsigned int b = 0; b < primaryHits_B.size(); b++) {

      double X1 = primaryHits_F[f].hitPositionX;
      double Y1 = primaryHits_F[f].hitPositionY;
      double Z1 = primaryHits_F[f].hitPositionZ;
      double X2 = primaryHits_B[b].hitPositionX;
      double Y2 = primaryHits_B[b].hitPositionY;
      double Z2= primaryHits_B[b].hitPositionZ;

     

	        if (fabs(primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg)>fFronttoBackTimingCut) continue;
		
		//std::cout<<"FOUND A COMBO"<<std::endl;
		double t0=(primaryHits_F[f].timeAvg+primaryHits_B[b].timeAvg)/2.f;

	        int tempId = 0;
		double xOffset=0;
    

		double trackStartPositionX_noSCE=trackStartPositionX_notCorrected;
		double trackEndPositionX_noSCE=trackEndPositionX_notCorrected;
		if (t0s.empty()){
		
		int RDOffset=0;
		if (!fMCCSwitch) RDOffset=111;
		double ticksOffset=0;

                if (!fMCCSwitch) ticksOffset = (t0+RDOffset)/25.f+detProp.GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);

                else if (fMCCSwitch) ticksOffset = (t0/500.f)+detProp.GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		
               xOffset=detProp.ConvertTicksToX(ticksOffset,allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
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
   if (fSCECorrection && SCE->EnableCalSpatialSCE()){
if(geom->PositionToTPCID(geo::Point_t(trackEndPositionX, trackEndPositionY, trackEndPositionZ)).deepestIndex()<13 && geom->PositionToTPCID(geo::Point_t(trackStartPositionX, trackStartPositionY, trackStartPositionZ)).deepestIndex()<13){ 
            auto const & posOffsets_F = SCE->GetCalPosOffsets(geo::Point_t(trackStartPositionX, trackStartPositionY, trackStartPositionZ), geom->PositionToTPCID(geo::Point_t(trackStartPositionX, trackStartPositionY, trackStartPositionZ)).deepestIndex());
            trackStartPositionX -= posOffsets_F.X();
            trackStartPositionY += posOffsets_F.Y();
            trackStartPositionZ += posOffsets_F.Z();
            auto const & posOffsets_B = SCE->GetCalPosOffsets(geo::Point_t(trackEndPositionX, trackEndPositionY, trackEndPositionZ), geom->PositionToTPCID(geo::Point_t(trackEndPositionX, trackEndPositionY, trackEndPositionZ)).deepestIndex());
            trackEndPositionX -= posOffsets_B.X();
            trackEndPositionY += posOffsets_B.Y();
            trackEndPositionZ += posOffsets_B.Z();
	}
    }





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

        double deltaY1 = (predictedHitPositionY1-Y1);

        double deltaY2 = (predictedHitPositionY2-Y2);



      tempId++;
     //iRecoTrack
      if (min_delta > std::abs(deltaX1)+std::abs(deltaX2) + std::abs(deltaY1)+std::abs(deltaY2) ){

	   min_delta=std::abs(deltaX1)+std::abs(deltaX2) + std::abs(deltaY1)+std::abs(deltaY2);
            best_XF = X1;
            best_YF = Y1;
            best_ZF = Z1;
            best_XB = X2;
            best_YB = Y2;
            best_ZB = Z2;
            best_dotProductCos = dotProductCos;

            best_deltaXF = deltaX1;
            best_deltaYF = deltaY1;
            best_deltaXB = deltaX2;
            best_deltaYB = deltaY2;
	    /*
	    best_trigXF=primaryHits_F[f].trigNumberX;
	    best_trigYF=primaryHits_F[f].trigNumberY;
	    best_trigXB=primaryHits_B[b].trigNumberX;
	    best_trigYB=primaryHits_B[b].trigNumberY;
	    */
	    best_T = t0;
	    bestHitIndex_F=f;
	    bestHitIndex_B=b;
	    if (!fMCCSwitch) best_T=(111.f+t0)*20.f;
	    // Added 111 tick CRT-CTB offset
          }
        }
      }
	X_F=best_XF; Y_F=best_YF; Z_F=best_ZF; X_B=best_XB; Y_B=best_YB; Z_B=best_ZB; int  f=bestHitIndex_F; int b=bestHitIndex_B;
 double t0=best_T;
 deltaX_F=best_deltaXF; deltaY_F=best_deltaYF; deltaX_B=best_deltaXB; deltaY_B=best_deltaYB;
	std::cout<<deltaX_F<<','<<deltaX_B<<','<<deltaY_F<<','<<deltaY_B<<','<<t0<<std::endl;
	if (deltaX_F==DBL_MAX && deltaY_F==DBL_MAX && deltaX_B==DBL_MAX && deltaY_B==DBL_MAX)	continue;


		double trackStartPositionX_noSCE=trackStartPositionX_notCorrected;
		double trackEndPositionX_noSCE=trackEndPositionX_notCorrected;
		double xOffset=0;
		if (t0s.empty()){
		double ticksOffset=0;


                ticksOffset=t0/500.f+detProp.GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		
               xOffset=detProp.ConvertTicksToX(ticksOffset,allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		
	 trackStartPositionX_noSCE=trackStartPositionX_notCorrected-xOffset;
         trackEndPositionX_noSCE=trackEndPositionX_notCorrected-xOffset;
	}


   double trackStartPositionX=trackStartPositionX_noSCE;
   double trackStartPositionY=trackStartPositionY_noSCE;
   double trackStartPositionZ=trackStartPositionZ_noSCE;

   double trackEndPositionX=trackEndPositionX_noSCE;
   double trackEndPositionY=trackEndPositionY_noSCE;
   double trackEndPositionZ=trackEndPositionZ_noSCE;


	TVector3 trackStart(trackStartPositionX, trackStartPositionY, trackStartPositionZ);
	TVector3 trackEnd(trackEndPositionX, trackEndPositionY, trackEndPositionZ);

        tracksPair tPair;
        tPair.tempId = tempId;
        tPair.recoId = iRecoTrack;
        tPair.deltaX_F = deltaX_F;
 
        tPair.deltaX_B = deltaX_B;
        tPair.deltaY_F = deltaY_F;
        tPair.deltaY_B = deltaY_B;

        tPair.dotProductCos=best_dotProductCos;

        tPair.moduleX1 = primaryHits_F[f].geoX;
        tPair.moduleX2 = primaryHits_B[b].geoX;
        tPair.moduleY1 = primaryHits_F[f].geoY;
        tPair.moduleY2 = primaryHits_B[b].geoY;
	tPair.timeDiff=primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg;
      tPair.adcX2=primaryHits_B[b].adcX;
      tPair.adcX1=primaryHits_F[f].adcX;
      tPair.adcY2=primaryHits_B[b].adcY;
      tPair.adcY1=primaryHits_F[f].adcY;
        tPair.stripX1 = primaryHits_F[f].stripX;
        tPair.stripX2 = primaryHits_B[b].stripX;
        tPair.stripY1 = primaryHits_F[f].stripY;
        tPair.stripY2 = primaryHits_B[b].stripY;
        tPair.X1 = X_F;
        tPair.Y1 = Y_F;
        tPair.Z1 = Z_F;
        tPair.X2 = X_B;
        tPair.Y2 = Y_B;
        tPair.Z2 = Z_B;

        tPair.endTPC=allHits[0]->WireID().TPC;
	
	tPair.xOffset=xOffset;
	tPair.t0=t0;
        tPair.trackStartPosition=trackStart;
	tPair.trackEndPosition=trackEnd;

	if (t0s.empty()) tPair.pandoraT0Check=0;
	else tPair.pandoraT0Check=1;
        allTracksPair.push_back(tPair);


      

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

	endTPC=allUniqueTracksPair[u].endTPC;

	recoPandoraT0Check=allUniqueTracksPair[u].pandoraT0Check;
	mccT0=allUniqueTracksPair[u].mccT0;
	measuredT0=allUniqueTracksPair[u].t0;
	measuredXOffset=allUniqueTracksPair[u].xOffset;
	CRT_TOF=allUniqueTracksPair[u].timeDiff;	
        if (fabs(allUniqueTracksPair[u].dotProductCos)>0.99 && fabs(deltaX_F)+fabs(deltaX_B)<40 && fabs(deltaY_F)+fabs(deltaY_B)<40) {
        //cout<<allUniqueTracksPair[u].timeDiff<<endl;
	//cout<<fabs(allUniqueTracksPair[u].dotProductCos)<<endl;

	  cout<< "Delta X and Y: "<<deltaX_F<<','<<deltaY_F<<','<<deltaX_B<<','<<deltaY_B<<endl;
	  cout<< "Predicted X and Y: "<< deltaX_F+X_F<<','<<deltaY_F+Y_F<<','<<deltaX_B+X_B<<','<<deltaY_B+Y_B<<endl;
	  cout<< "Detected X and Y: "<< X_F<<','<<Y_F<<','<<X_B<<','<<Y_B<<endl;
	  cout<<moduleX_F<<','<<stripX_F<<endl;
	  cout<<moduleY_F<<','<<stripY_F<<endl;
	  cout<<moduleX_B<<','<<stripX_B<<endl;
	  cout<<moduleY_B<<','<<stripY_B<<endl;
	  cout<<"ADC Values: "<<adcX_F<<','<<adcY_F<<','<<adcX_B<<','<<adcY_B<<endl;

 
	int iRecoTrack=trackID;
         averageSignedDistanceXZ=0;
	 averageSignedDistanceYZ=0;
         averageSignedDistance=0; 		
	const size_t lastPoint=trackList[iRecoTrack]->NumberTrajectoryPoints();
        for (size_t trackpoint = 0; trackpoint < lastPoint; ++trackpoint) {
	  double trackPosX=trackList[iRecoTrack] -> LocationAtPoint(trackpoint).X()-measuredXOffset;
	  double trackPosY=trackList[iRecoTrack] -> LocationAtPoint(trackpoint).Y();
	  double trackPosZ=trackList[iRecoTrack] -> LocationAtPoint(trackpoint).Z();

	   if (trackPosY==-999) continue;
	   if (trackPosX==-999) continue;
	TVector3 trackPos(trackPosX, trackPosY, trackPosZ);
			double distanceYZ = signedPointToLineDistance( Y_F,Z_F, Y_B,Z_B, trackPos.Y(), trackPos.Z() ); //only the Y and Z of trackpos will be used
			double distanceXZ = signedPointToLineDistance( X_F,Z_F, X_B,Z_B, trackPos.X(), trackPos.Z() );


				double distance=signed3dDistance(X_F,Y_F,Z_F,X_B,Y_B,Z_B, trackPos);

	
           trkPosX=trackPosX;
	  trkPosY=trackPosY;
	  trkPosZ=trackPosZ;    
	  
	   displXZ=distanceXZ;
	 
	displYZ	=distanceYZ; 

			averageSignedDistance += distance/(lastPoint+1);	
		        averageSignedDistanceYZ += distanceYZ/(lastPoint+1);
			averageSignedDistanceXZ += distanceXZ/(lastPoint+1);
//averageSignedDistanceXY += distanceXY/(lastPoint+1);   
	fTrackInfo->Fill();                                                               
	
		    }

	
      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(iRecoTrack);
      //std::vector< art::Ptr<recob::Hit> > allHits =  hitsFromTrack.at(iRecoTrack);
      for(size_t ical = 0; ical<calos.size(); ++ical){
	if (calos[ical]->PlaneID().Plane!=2) continue;
	
	const size_t NHits=calos[ical]->dEdx().size();
	for(size_t iHit = 0; iHit < NHits; ++iHit){
	   
	  const auto& TrkPos = (calos[ical] -> XYZ()[iHit]);
	
          trkdqdx=(calos[ical]->dQdx()[iHit]);
	  size_t tpIndex=(calos[ical]->TpIndices()[iHit]);
	  if (hitlist[tpIndex]->Multiplicity()>1) continue;
          trkhitt0=hitlist[tpIndex]->PeakTime()*500.f;
	  trkhitIntegral=hitlist[tpIndex]->Integral();
	  crtt0=allUniqueTracksPair[u].t0;
	  trkhitx=TrkPos.X()-measuredXOffset;
	  trkhity=TrkPos.Y();
	  trkhitz=TrkPos.Z();
	  WireID=hitlist[tpIndex]->WireID().Wire;
	  TPCID=hitlist[tpIndex]->WireID().TPC;
	  sumADC=hitlist[tpIndex]->SummedADC();
	  sigmaHit=hitlist[tpIndex]->SigmaIntegral();
	  rangeTime=hitlist[tpIndex]->EndTick()-hitlist[tpIndex]->StartTick();
	  fCRTdQTree->Fill();
	}
 
      } 
	
	
	fCRTTree->Fill();
	   

        }
      }
    }

 }


// Setup CRT 
void CRT::TwoCRTMatching::beginJob() {
	art::ServiceHandle<art::TFileService> fileServiceHandle;
       fCRTTree = fileServiceHandle->make<TTree>("Displacement", "track by track info");
        fMCCMuon= fileServiceHandle->make<TTree>("MCCTruths", "event by event info");

        fTrackInfo= fileServiceHandle->make<TTree>("TrackInfo", "track by track info");
	fCRTdQTree=fileServiceHandle->make<TTree>("CRTdQ", "track by track info");
	fCRTTree->Branch("nEvents", &nEvents, "fnEvents/I");
	fCRTTree->Branch("hRun", &run, "run/I");
	fCRTTree->Branch("hSubRun", &subRun, "subRun/I");
	fCRTTree->Branch("htrackID", &trackID, "trackID/I");
        fCRTTree->Branch("hmccTrackId", &mccTrackId, "mccTrackId/I");

	fMCCMuon->Branch("nEvents", &nEvents, "nEvents/I");
	fMCCMuon->Branch("mccCRTStartX",&mccCRTStartX,"mccCRTStartX/D");
	fMCCMuon->Branch("mccCRTStartY",&mccCRTStartY,"mccCRTStartY/D");
	fMCCMuon->Branch("mccCRTStartZ",&mccCRTStartZ,"mccCRTStartZ/D");
	fMCCMuon->Branch("mccCRTEndX",&mccCRTEndX,"mccCRTEndX/D");
	fMCCMuon->Branch("mccCRTEndY",&mccCRTEndY,"mccCRTEndY/D");
	fMCCMuon->Branch("mccCRTEndZ",&mccCRTEndZ,"mccCRTEndZ/D");
	fMCCMuon->Branch("truthEnergy", &truthEnergy, "truthEnergy/D");
	fMCCMuon->Branch("hRun", &run, "run/I");
	fMCCMuon->Branch("hSubRun", &subRun, "subRun/I");
        fMCCMuon->Branch("hmccTrackId", &mccTrackId, "mccTrackId/I");
        fMCCMuon->Branch("hcandidateCRT", &candidateCRT, "candidateCRT/I");
	fMCCMuon->Branch("hmccT0", &mccT0, "mccT0/D");


	fCRTdQTree->Branch("trkhitx",&trkhitx,"trkhitx/D");
	fCRTdQTree->Branch("trkhity",&trkhity,"trkhity/D");
	fCRTdQTree->Branch("trkhitz",&trkhitz,"trkhitz/D");
	fCRTdQTree->Branch("trkdqdx",&trkdqdx,"trkdqdx/D");
	fCRTdQTree->Branch("trkhitt0",&trkhitt0,"trkhitt0/D");
	fCRTdQTree->Branch("trkhitIntegral",&trkhitIntegral,"trkhitIntegral/I");
	fCRTdQTree->Branch("crtt0",&crtt0,"crtt0/D");
	fCRTdQTree->Branch("hX_F", &X_F, "X_F/D");
	fCRTdQTree->Branch("hX_B", &X_B, "X_B/D");
	fCRTdQTree->Branch("hY_F", &Y_F, "Y_F/D");
	fCRTdQTree->Branch("hY_B", &Y_B, "Y_B/D");
	fCRTdQTree->Branch("hZ_F", &Z_F, "Z_F/D");
	fCRTdQTree->Branch("hZ_B", &Z_B, "Z_B/D");

	fCRTdQTree->Branch("hTPCID", &TPCID, "TPCID/I");
	fCRTdQTree->Branch("hWireID", &WireID, "WireID/I");
	fCRTdQTree->Branch("hsumADC",&sumADC,"sumADC/I");
	fCRTdQTree->Branch("hrangeTime",&rangeTime,"rangeTime/I");
	fCRTdQTree->Branch("hsigmaHit",&sigmaHit,"sigmaHit/D");

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

	fCRTTree->Branch("hSCECorrectX_F", &SCECorrectX_F, "SCECorrectX_F/D");
	fCRTTree->Branch("hSCECorrectY_F", &SCECorrectY_F, "SCECorrectY_F/D");
	fCRTTree->Branch("hSCECorrectZ_F", &SCECorrectZ_F, "SCECorrectZ_F/D");
	fCRTTree->Branch("hSCECorrectX_B", &SCECorrectX_B, "SCECorrectX_B/D");
	fCRTTree->Branch("hSCECorrectY_B", &SCECorrectY_B, "SCECorrectY_B/D");
	fCRTTree->Branch("hSCECorrectZ_B", &SCECorrectZ_B, "SCECorrectZ_B/D");
	fCRTTree->Branch("hendTPC", &endTPC, "endTPC/I");

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
	fCRTTree->Branch("hstripY_F", &stripY_F, "stripY_F/I");
	fCRTTree->Branch("hstripY_B", &stripY_B, "stripY_B/I");

	fCRTTree->Branch("hadcX_F", &adcX_F, "adcX_F/I");
	fCRTTree->Branch("hadcY_F", &adcY_F, "adcY_F/I");
	fCRTTree->Branch("hadcX_B", &adcX_B, "adcX_B/I");
	fCRTTree->Branch("hadcY_B", &adcY_B, "adcY_B/I");

	fCRTTree->Branch("hmeasuredT0", &measuredT0, "measuredT0/D");
	fCRTTree->Branch("hmeasuredXOffset", &measuredXOffset, "measuredXOffset/D");
	fCRTTree->Branch("hmccT0", &mccT0, "mccT0/D");
	fCRTTree->Branch("hrecoPandoraT0Check", &recoPandoraT0Check, "recoPandoraT0Check/B");



	fTrackInfo->Branch("nEvents", &nEvents, "nEvents/I");
	fTrackInfo->Branch("hdotCos", &dotCos, "dotCos/D");


	fTrackInfo->Branch("htrkPosX", &trkPosX, "trkPosX/D");
	fTrackInfo->Branch("htrkPosY", &trkPosY, "trkPosY/D");
	fTrackInfo->Branch("htrkPosZ", &trkPosZ, "trkPosZ/D");
	fTrackInfo->Branch("hdisplXZ", &displXZ, "displXZ/D");
	fTrackInfo->Branch("hdisplYZ", &displYZ, "displYZ/D");
        fTrackInfo->Branch("hmeasuredXOffset", &measuredXOffset, "measuredXOffset/D");

	fTrackInfo->Branch("hX_F", &X_F, "X_F/D");
	fTrackInfo->Branch("hX_B", &X_B, "X_B/D");
	fTrackInfo->Branch("hY_F", &Y_F, "Y_F/D");
	fTrackInfo->Branch("hY_B", &Y_B, "Y_B/D");
	fTrackInfo->Branch("hZ_F", &Z_F, "Z_F/D");
	fTrackInfo->Branch("hZ_B", &Z_B, "Z_B/D");




}
// Endjob actions
void CRT::TwoCRTMatching::endJob() 
{



}














DEFINE_ART_MODULE(CRT::TwoCRTMatching)
