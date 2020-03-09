////////////////////////////////////////////////////////////////////////
// Class:       TwoCRTReco
// Plugin Type: Analyzer (art v2_11_02)
// File:        TwoCRTReco_module.cc
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
  class TwoCRTReco;
}


class CRT::TwoCRTReco: public art::EDAnalyzer {
  public: // Setup functions and variables
    explicit TwoCRTReco(fhicl::ParameterSet
      const & p);
  std::string fTrackModuleLabel = "pandoraTrack";
  //std::string fTrackModuleLabel = "pmtrack";
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  TwoCRTReco(TwoCRTReco
    const & ) = delete;
  TwoCRTReco(TwoCRTReco && ) = delete;
  TwoCRTReco & operator = (TwoCRTReco
    const & ) = delete;
  TwoCRTReco & operator = (TwoCRTReco && ) = delete;
  void analyze(art::Event const & e) override;
  //bool filter(art::Event const & e) override;
// Declare functions and variables for validation
  bool moduleMatcher(int module1, int module2);
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

    int run, subRun;
    bool fMCCSwitch;
    bool fCTBTriggerOnly;
    bool fSCECorrection;
    bool fModuleSwitch;
    int fADCThreshold;
    int fModuletoModuleTimingCut;
    int fFronttoBackTimingCut;


  //double averageSignedDistanceXY;
    
    int adcX_F, adcX_B, adcY_F, adcY_B;
    double X_F, X_B, Y_F, Y_B, Z_F, Z_B;
    double trackX1, trackX2, trackY1, trackY2, trackZ1, trackZ2;
    int moduleX_F, moduleX_B, moduleY_F, moduleY_B;
    int stripX_F, stripX_B, stripY_F, stripY_B;
    double measuredT0;
    double CRT_TOF;
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
  }
  tracksPair;





  std::vector < recoHits > primaryHits_F;
  std::vector < recoHits > primaryHits_B;

  std::vector < tempHits > tempHits_F;
  std::vector < tempHits > tempHits_B;
  std::vector < tracksPair > allTracksPair;

};

CRT::TwoCRTReco::TwoCRTReco(fhicl::ParameterSet
    const & p):
  EDAnalyzer(p), fCRTLabel(p.get < art::InputTag > ("CRTLabel")),  fCTBLabel(p.get<art::InputTag>("CTBLabel")) {
    consumes < std::vector < CRT::Trigger >> (fCRTLabel);
    consumes < std::vector < art::Assns < sim::AuxDetSimChannel, CRT::Trigger >>> (fCRTLabel); 
  fMCCSwitch=(p.get<bool>("MCC"));
  fCTBTriggerOnly=(p.get<bool>("CTBOnly"));
  fSCECorrection=(p.get<bool>("SCECorrection"));
  }




// v6 Geo Channel Map
bool CRT::TwoCRTReco::moduleMatcher(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 6 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 26 || module2 == 27)) || (module1 == 31 && (module2 == 26 || module2 == 27)) || (module1 == 7 && (module2 == 12 || module2 == 13)) || (module1 == 15 && (module2 == 12 || module2 == 13)) || (module1 == 18 && (module2 == 24 || module2 == 25)) || (module1 == 30 && (module2 == 24 || module2 == 25)) || (module1 == 1 && (module2 == 4 || module2 == 5)) || (module1 == 9 && (module2 == 4 || module2 == 5)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 0 && (module2 == 2 || module2 == 3)) || (module1 == 8 && (module2 == 2 || module2 == 3)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23))) return 1;
  else return 0;

}


int CRT::TwoCRTReco::moduletoCTB(int module2, int module1){
  if (module1 == 15 && module2 == 12 ) return 15;
  else if (module1 == 15 &&  module2 == 13) return 10;
  else if (module1 == 7 &&  module2 == 12) return 8;
  else if (module1 == 7 &&  module2 == 13) return 9;
  else if (module1 == 19 &&  module2 == 26) return 4;
  else if (module1 == 19 &&  module2 == 27) return 13;
  else if (module1 == 31 &&  module2 == 26) return 3;
  else if (module1 == 31 &&  module2 == 27) return 2;
  else if (module1 == 30 &&  module2 == 24) return 1;
  else if (module1 == 30 &&  module2 == 25) return 0;
  else if (module1 == 18 &&  module2 == 24) return 12;
  else if (module1 == 18 &&  module2 == 25) return 11;
  else if (module1 == 6  &&  module2 == 11) return 7;
  else if (module1 == 6 &&  module2 == 10) return 6;
  else if (module1 == 14  &&  module2 == 11) return 14;
  else if (module1 == 14 &&  module2 == 10) return 5;
  else if (module1 == 0 &&  module2 == 3) return 25;
  else if (module1 == 0 &&  module2 == 2) return 24;
  else if (module1 == 8 &&  module2 == 3) return 26;
  else if (module1 == 8 &&  module2 == 2) return 30;
  else if (module1 == 17 &&  module2 == 23) return 27;
  else if (module1 == 17 &&  module2 == 22) return 28;
  else if (module1 == 29 &&  module2 == 23) return 16;
  else if (module1 == 29 &&  module2 == 22) return 17;
  else if (module1 == 28 &&  module2 == 21) return 18;
  else if (module1 == 28 &&  module2 == 20) return 19;
  else if (module1 == 16 &&  module2 == 21) return 29;
  else if (module1 == 16 &&  module2 == 20) return 20;
  else if (module1 == 9 &&  module2 == 5) return 31;
  else if (module1 == 9 &&  module2 == 4) return 21;
  else if (module1 == 1 &&  module2 == 5) return 23;
  else if (module1 == 1 &&  module2 == 4) return 22;
  else return -1;
}


void CRT::TwoCRTReco::createPNG(TH1D * histo) {
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
double CRT::TwoCRTReco::setAngle(double angle) {
  if (angle < 0) {
    angle += 3.14159265359;
  }
  angle *= 180.00 / 3.14159265359;
  return angle;
}

//bool CRT::TwoCRTReco::filter(art::Event const & event)
void CRT::TwoCRTReco::analyze(art::Event
  const & event) // Analysis module

{

    nEvents++;	
    run=event.run();
    subRun=event.subRun();
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

   // Find CTB pixels
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
	//cout<<nEvents<<" TJYang Pixels: "<<pixel0<<","<<pixel1<<endl;
	}
	else if (fCTBTriggerOnly) return;
	      }
	  }
	}

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

	std::cout<<"Number of Hits: "<<primaryHits_F.size()<<','<<primaryHits_B.size()<<std::endl;

	std::cout<<"CTB Pixes: "<<pixel0<<','<<pixel1<<std::endl;
    int tempId = 0;
    adcX_F=0;
    adcX_B=0;
    adcY_F=0;
    adcY_B=0;
    int crtPixel0=-1;
    int crtPixel1=-1;
    for (unsigned int f = 0; f < primaryHits_F.size(); f++) {
      if (pixel0==-1 || pixel1==-1) break; 
      for (unsigned int b = 0; b < primaryHits_B.size(); b++) {
      if (fabs(primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg)>fFronttoBackTimingCut) continue;
      crtPixel0=moduletoCTB(primaryHits_F[f].geoX,primaryHits_F[f].geoY);
      crtPixel1=moduletoCTB(primaryHits_B[b].geoX, primaryHits_B[b].geoY);
      if (crtPixel0!=pixel0 || crtPixel1!=pixel1) continue;	
      //std::cout<<"HEY"<<std::endl;
      if (adcX_F<primaryHits_F[f].adcX && adcY_F<primaryHits_F[f].adcY && adcX_B<primaryHits_B[b].adcX  && adcY_B<primaryHits_B[b].adcY){
      X_F = primaryHits_F[f].hitPositionX;
      Y_F = primaryHits_F[f].hitPositionY;
      Z_F = primaryHits_F[f].hitPositionZ;
      X_B = primaryHits_B[b].hitPositionX;
      Y_B = primaryHits_B[b].hitPositionY;
      Z_B= primaryHits_B[b].hitPositionZ;
      adcX_F=primaryHits_F[f].adcX;
      adcX_B=primaryHits_B[b].adcX;
      adcY_F=primaryHits_F[f].adcY;
      adcY_B=primaryHits_B[b].adcY;
      moduleX_F=primaryHits_F[f].geoX;
      moduleY_F=primaryHits_F[f].geoY;
      moduleX_B=primaryHits_B[b].geoX;
      moduleY_B=primaryHits_B[b].geoY;

      stripX_F=primaryHits_F[f].stripX;
      stripY_F=primaryHits_F[f].stripY;
      stripX_B=primaryHits_B[b].stripX;
      stripY_B=primaryHits_B[b].stripY;
      //std::cout<<"FOUND A COMBO"<<std::endl;
      measuredT0=(primaryHits_F[f].timeAvg+primaryHits_B[b].timeAvg)/2.f;
	CRT_TOF=(primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg);
	tempId++;

	}	
        }
      }
    if (adcX_F>0) fCRTTree->Fill();
    // Filter return if (adcX_F==0) return false;
    // Filter return
    // return true;
 }


// Setup CRT 
void CRT::TwoCRTReco::beginJob() {
	art::ServiceHandle<art::TFileService> fileServiceHandle;
       fCRTTree = fileServiceHandle->make<TTree>("CRTCandidate", "track by track info");

	fCRTTree->Branch("nEvents", &nEvents, "fnEvents/I");
	fCRTTree->Branch("hRun", &run, "run/I");



	fCRTTree->Branch("hX_F", &X_F, "X_F/D");
	fCRTTree->Branch("hX_B", &X_B, "X_B/D");
	fCRTTree->Branch("hY_F", &Y_F, "Y_F/D");
	fCRTTree->Branch("hY_B", &Y_B, "Y_B/D");
	fCRTTree->Branch("hZ_F", &Z_F, "Z_F/D");
	fCRTTree->Branch("hZ_B", &Z_B, "Z_B/D");


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





}
// Endjob actions
void CRT::TwoCRTReco::endJob() 
{



}














DEFINE_ART_MODULE(CRT::TwoCRTReco)
