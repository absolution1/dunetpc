////////////////////////////////////////////////////////////////////////
// Class:       CRTTimingValidation
// Plugin Type: analyzer (art v2_11_02)
// File:        CRTTimingValidation_module.cc
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

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "lardataobj/RawData/RDTimeStamp.h"


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
using namespace std;   // Namespaces established to make life easier
using namespace ROOT::Math;
namespace CRT {
  class CRTTimingValidation;

}



class CRT::CRTTimingValidation: public art::EDAnalyzer {
  public: // Setup functions and variables
    explicit CRTTimingValidation(fhicl::ParameterSet
      const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  CRTTimingValidation(CRTTimingValidation
    const & ) = delete;
  CRTTimingValidation(CRTTimingValidation && ) = delete;
  CRTTimingValidation & operator = (CRTTimingValidation
    const & ) = delete;
  CRTTimingValidation & operator = (CRTTimingValidation && ) = delete;
  void analyze(art::Event
    const & e) override;
// Declare functions and variables for validation
  bool moduleCheckX(int module);
  int moduletoCTB(int module2, int module1);
  void beginJob() override;
  void endJob() override;

  //Parameters for reading in CRT::Triggers and associated AuxDetSimChannels.
  art::InputTag fCRTLabel; //Label for the module that produced 
  art::InputTag fCTBLabel; 

  TTree * fCRTTree;
  TTree * fCRTTreeF;
  TTree * fCRTTreeB;
    bool fMCCSwitch;
    bool fModuleSwitch;
    int fADCThreshold;
	int nEvents=0;
	int matchedCTBtoCRT=0;
	int ctb_F, ctb_B;

    double CRT_TOF;
    int T0Offset_F;
    int T0Offset_B;
    double T_F, T_B;
    int moduleX_F, moduleX_B, moduleY_F, moduleY_B;
    int stripX_F, stripX_B, stripY_F, stripY_B;
    double RDminDeltaT_F, RDminDeltaT_B;

  typedef struct // Structures for arrays to move hits from raw to reco to validation
  {

    int module;
    int64_t triggerTime;
  }
  tempTrigger;

  typedef struct // Structures for arrays to move hits from raw to reco to validation
  {

    int moduleX;
    int moduleY;
    int triggerDiff;
    double triggerTimeAvg;
    int RDDeltaT;
  }
  tempHits;

  typedef struct // Structures for arrays to move hits from raw to reco to validation
  {

    int module1;
    int module2;
    int ctbTimestamp;
  }
  ctbHits;


  std::vector < tempTrigger > trigger_F_X;
  std::vector < tempTrigger > trigger_F_Y;
  std::vector < tempTrigger > trigger_B_X;
  std::vector < tempTrigger > trigger_B_Y;

  std::vector < tempHits > hits_F;
  std::vector < tempHits > hits_B;
  std::vector <ctbHits> ctbTriggers;



};

CRT::CRTTimingValidation::CRTTimingValidation(fhicl::ParameterSet
    const & p):
  EDAnalyzer(p), fCRTLabel(p.get < art::InputTag > ("CRTLabel")), fCTBLabel(p.get<art::InputTag>("CTBLabel")) {
    consumes < std::vector < CRT::Trigger >> (fCRTLabel);
    consumes < std::vector < art::Assns < sim::AuxDetSimChannel, CRT::Trigger >>> (fCRTLabel); // CRT art consumables
  fMCCSwitch=(p.get<bool>("MCC"));

  }




bool CRT::CRTTimingValidation::moduleCheckX(int module){
if (module==4 || module==5 || module==6 || module==7 || module==8 || module==9 || module==10 || module==11 || module==20 || module==21 || module==22 || module==23 || module==24 || module==25 || module==26 || module==27) return 1;
else return 0;
}


int CRT::CRTTimingValidation::moduletoCTB(int module2, int module1){
  if (module1 == 13 && module2 == 6 ) return 15;
  else if (module1 == 13 &&  module2 == 7) return 10;
  else if (module1 == 1 &&  module2 == 6) return 8;
  else if (module1 == 1 &&  module2 == 7) return 9;
  else if (module1 == 16 &&  module2 == 20) return 4;
  else if (module1 == 16 &&  module2 == 21) return 13;
  else if (module1 == 28 &&  module2 == 20) return 3;
  else if (module1 == 28 &&  module2 == 21) return 2;
  else if (module1 == 29 &&  module2 == 22) return 1;
  else if (module1 == 29 &&  module2 == 23) return 0;
  else if (module1 == 17 &&  module2 == 22) return 12;
  else if (module1 == 17 &&  module2 == 23) return 11;
  else if (module1 == 0  &&  module2 == 5) return 7;
  else if (module1 == 0 &&  module2 == 4) return 6;
  else if (module1 == 12  &&  module2 == 5) return 14;
  else if (module1 == 12 &&  module2 == 4) return 5;
  else if (module1 == 3 &&  module2 == 8) return 25;
  else if (module1 == 3 &&  module2 == 9) return 24;
  else if (module1 == 15 &&  module2 == 8) return 26;
  else if (module1 == 15 &&  module2 == 9) return 31;
  else if (module1 == 18 &&  module2 == 26) return 27;
  else if (module1 == 18 &&  module2 == 27) return 28;
  else if (module1 == 30 &&  module2 == 26) return 16;
  else if (module1 == 30 &&  module2 == 27) return 17;
  else if (module1 == 31 &&  module2 == 24) return 18;
  else if (module1 == 31 &&  module2 == 25) return 19;
  else if (module1 == 19 &&  module2 == 24) return 29;
  else if (module1 == 19 &&  module2 == 25) return 20;
  else if (module1 == 14 &&  module2 == 10) return 30;
  else if (module1 == 14 &&  module2 == 11) return 21;
  else if (module1 == 2 &&  module2 == 10) return 23;
  else if (module1 == 2 &&  module2 == 11) return 22;
  else return -1;
}





void CRT::CRTTimingValidation::analyze(art::Event
  const & event) // Analysis module
{
//cout<<"Setting everything up"<<endl;
   art::ValidHandle<std::vector<raw::RDTimeStamp>> timingHandle = event.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");

 using timestamp_t = int64_t;


	const raw::RDTimeStamp& timeStamp = timingHandle->at(0);
	if(timeStamp.GetFlags()!= 13) return;

  if (fMCCSwitch){
    fModuleSwitch=1;
    fADCThreshold=800;
    
}
  else {
    fModuleSwitch=0;
    fADCThreshold=50;


}

 


  hits_F.clear();
  hits_B.clear();
  trigger_F_X.clear();

  trigger_F_Y.clear();

  trigger_B_Y.clear();

  trigger_B_X.clear();
  ctbTriggers.clear();
 // Arrays to compile hits and move them through

    //allTracksPair.clear();

  //Get triggers
  //cout << "Getting triggers" << endl;
  const auto & triggers = event.getValidHandle < std::vector < CRT::Trigger >> (fCRTLabel);

  art::FindManyP < sim::AuxDetSimChannel > trigToSim(triggers, event, fCRTLabel);

  //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
  art::ServiceHandle < geo::Geometry > geom;

  //Mapping from channel to trigger
  std::unordered_map < size_t, double > prevTimes;
  int hitID = 0;
  //cout << "Looking for hits in Triggers" << endl;

  for (const auto & trigger: * triggers) {

        tempTrigger tTrigger;

        tTrigger.module = trigger.Channel(); // Values to add to array
	tTrigger.triggerTime=trigger.Timestamp();

	//cout<<ctbPixels[0]<<endl;

        const auto & trigGeo = geom -> AuxDet(trigger.Channel()); // Get geo  
        const auto & csens = trigGeo.SensitiveVolume(0);
        const auto center = csens.GetCenter();
        if (center.Z() < 100 and moduleCheckX(trigger.Channel())==1) trigger_F_X.push_back(tTrigger);
        else if (moduleCheckX(trigger.Channel())==1) trigger_B_X.push_back(tTrigger);
        else if (center.Z() < 100 and moduleCheckX(trigger.Channel())!=1) trigger_F_Y.push_back(tTrigger);
        else trigger_B_Y.push_back(tTrigger);
        hitID++;
      }
for (unsigned int i=0; i < trigger_F_X.size(); i++){
for (unsigned int j=0; j < trigger_F_Y.size(); j++){
T0Offset_F=trigger_F_X[i].triggerTime-trigger_F_Y[j].triggerTime;
moduleX_F=trigger_F_X[i].module;
moduleY_F=trigger_F_Y[j].module;



        tempHits tHits;

        tHits.moduleX = trigger_F_X[i].module; // Values to add to array
	tHits.moduleY=trigger_F_Y[j].module;
	tHits.triggerDiff=T0Offset_F;
	tHits.triggerTimeAvg=(trigger_F_X[i].triggerTime+trigger_F_Y[j].triggerTime)/2.;
	int minDeltaT=9999;      
	for(const auto& time: *timingHandle)
      {
	
        const timestamp_t& rawTime = time.GetTimeStamp();
        const auto deltaT = rawTime - (trigger_F_X[i].triggerTime);
	//cout<<rawTime<<','<<deltaT<<endl;
        if(fabs(deltaT) < fabs(minDeltaT)) minDeltaT = deltaT;
      }
	tHits.RDDeltaT=minDeltaT;
	RDminDeltaT_F=minDeltaT;
if (fabs(T0Offset_F)<5) {fCRTTreeF->Fill(); hits_F.push_back(tHits);}

}
}

for (unsigned int i=0; i<trigger_B_X.size(); i++){
for (unsigned int j=0; j<trigger_B_Y.size(); j++){
T0Offset_B=trigger_B_X[i].triggerTime-trigger_B_Y[j].triggerTime;
moduleX_B=trigger_B_X[i].module;
moduleY_B=trigger_B_Y[j].module;



        tempHits tHits;

        tHits.moduleX = trigger_B_X[i].module; // Values to add to array
	tHits.moduleY=trigger_B_Y[j].module;
	tHits.triggerDiff=T0Offset_B;
	tHits.triggerTimeAvg=(trigger_B_X[i].triggerTime+trigger_B_Y[j].triggerTime)/2.;
	int minDeltaT=9999;  
	//cout<<"Timing offset:"<<trigger_B_X[i].triggerTime-trigger_F_X[i].triggerTime<<endl;   
	for(const auto& time: *timingHandle)
      {
        const timestamp_t& rawTime = time.GetTimeStamp();
        const auto deltaT = rawTime - (trigger_B_X[i].triggerTime);
        if(fabs(deltaT) < fabs(minDeltaT)) minDeltaT=deltaT;
      }
	
	tHits.RDDeltaT=minDeltaT;
	RDminDeltaT_B=minDeltaT;
if (fabs(T0Offset_B)<5) {fCRTTreeB->Fill(); hits_B.push_back(tHits);}
}
}

       const auto& pdspctbs = *event.getValidHandle<std::vector<raw::ctb::pdspctb>>(fCTBLabel);
    std::vector<int> uS, dS;

//cout<<pdspctbs.size()<<endl;




	const size_t npdspctbs = pdspctbs.size();
	for(size_t j=0;j<npdspctbs;++j)
	  {
	    const std::vector<raw::ctb::Trigger> HLTriggers = pdspctbs[j].GetHLTriggers();
	    const std::vector<raw::ctb::ChStatus> chs = pdspctbs[j].GetChStatusAfterHLTs();
for (size_t k=0; k<HLTriggers.size(); ++k)
	      { 



	

	
		dS.clear(); uS.clear();
		//cout<<chs[k].timestamp<<endl;
		int num = chs[k].crt;
		//cout<<num<<endl;
	
		const std::string binary = std::bitset<32>(num).to_string();	

	//std::vector<CRT::Trigger> inWindow(triggers->size());
	
       //const auto crtMask = fGeom->toCTBMask(inWindow);
	//cout<<crtMask<<','<<chs[k].crt<<endl;
      //constexpr size_t nBits = 32;
      //std::bitset<nBits> diff(crtMask ^ status.crt), crtOnly(crtMask & (~status.crt)), ctbOnly(status.crt & (~crtMask));
	const auto crtmask=chs[k].crt;
        int pixel0 = -1;
        int pixel1 = -1;
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
        if (pixel0!=-1 && pixel1!=1) {



	cout<<nEvents<<" TJYang Pixels: "<<pixel0<<","<<pixel1<<endl;
	ctb_F=pixel0;
	ctb_B=pixel1;


	}
	else return;

		

	

	      }
	  }
int loopTimer=0;
for (unsigned int i=0; i<hits_F.size(); i++){
for (unsigned int j=0; j<hits_B.size(); j++){
moduleX_B=hits_B[j].moduleX;
moduleY_B=hits_B[j].moduleY;

moduleX_F=hits_F[i].moduleX;
moduleY_F=hits_F[i].moduleY;
T_F=hits_F[i].triggerTimeAvg;
T_B=hits_B[j].triggerTimeAvg;
CRT_TOF=hits_B[j].triggerTimeAvg-hits_F[i].triggerTimeAvg;


if (fabs(CRT_TOF)<5 && ctb_F==moduletoCTB(hits_F[i].moduleX,hits_F[i].moduleY) && ctb_B==moduletoCTB(hits_B[j].moduleX,hits_B[j].moduleY)) {cout<<nEvents<<" CRT to CTB Front: "<<hits_F[i].moduleX<<','<<hits_F[i].moduleY<<','<<moduletoCTB(hits_F[i].moduleX,hits_F[i].moduleY)<<endl; 
cout<<nEvents<<" CRT to CTB Back: "<<hits_B[j].moduleX<<','<<hits_B[j].moduleY<<','<<moduletoCTB(hits_B[j].moduleX,hits_B[j].moduleY)<<endl;  fCRTTree->Fill();
matchedCTBtoCRT++;

cout<<matchedCTBtoCRT<<endl;
}
loopTimer++;
}
}

nEvents++;
 }


// Setup CRT 
void CRT::CRTTimingValidation::beginJob() {
	art::ServiceHandle<art::TFileService> fileServiceHandle;
       fCRTTreeF = fileServiceHandle->make<TTree>("T_F", "event by event info");
       fCRTTreeB = fileServiceHandle->make<TTree>("T_B", "event by event info");

       fCRTTree = fileServiceHandle->make<TTree>("Matching TOF", "event by event info");


	fCRTTreeF->Branch("nEvents", &nEvents, "nEvents/I");
	fCRTTreeB->Branch("nEvents", &nEvents, "nEvents/I");
	fCRTTree->Branch("nEvents", &nEvents, "nEvents/I");
	fCRTTreeF->Branch("RDminDeltaT_F", &RDminDeltaT_F, "RDminDeltaT_F/D");
	fCRTTreeB->Branch("minDeltaT_B", &RDminDeltaT_B, "RDminDeltaT_B/D");

	fCRTTreeF->Branch("T0Offset_F", &T0Offset_F, "T0Offset_F/I");
	fCRTTreeB->Branch("T0Offset_B", &T0Offset_B, "T0Offset_B/I");


	fCRTTreeF->Branch("hmoduleX_F", &moduleX_F, "moduleX_F/I");
	fCRTTreeB->Branch("hmoduleX_B", &moduleX_B, "moduleX_B/I");

	fCRTTreeF->Branch("hmoduleY_F", &moduleY_F, "moduleY_F/I");
	fCRTTreeB->Branch("hmoduleY_B", &moduleY_B, "moduleY_B/I");

	fCRTTree->Branch("CRT_TOF", &CRT_TOF, "CRT_TOF/D");
fCRTTree->Branch("T_B", &T_B, "T_B/D");
fCRTTree->Branch("T_F", &T_F, "T_F/D");


	fCRTTree->Branch("hmoduleX_F", &moduleX_F, "moduleX_F/I");
	fCRTTree->Branch("hmoduleX_B", &moduleX_B, "moduleX_B/I");
	fCRTTree->Branch("hmoduleY_F", &moduleY_F, "moduleY_F/I");
	fCRTTree->Branch("hmoduleY_B", &moduleY_B, "moduleY_B/I");


}
// Endjob actions
void CRT::CRTTimingValidation::endJob() 
{

cout<<matchedCTBtoCRT<<endl;

}














DEFINE_ART_MODULE(CRT::CRTTimingValidation)
