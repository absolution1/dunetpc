////////////////////////////////////////////////////////////////////////
// Class:       SingleCRTMatchingProducer
// Plugin Type: producer (art v2_10_03)
// File:        SingleCRTMatchingProducer_module.cc
//
// Generated at Wed Jun 27 04:09:39 2018 by Andrew Olivier using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "canvas/Persistency/Common/Ptr.h"
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
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

//LArSoft includes

#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nug4/ParticleNavigation/ParticleList.h"


#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"



#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

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
  class SingleCRTMatchingProducer;
}


class CRT::SingleCRTMatchingProducer : public art::EDProducer {
public:
  explicit SingleCRTMatchingProducer(fhicl::ParameterSet const & p);

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SingleCRTMatchingProducer(SingleCRTMatchingProducer const &) = delete;
  SingleCRTMatchingProducer(SingleCRTMatchingProducer &&) = delete;
  SingleCRTMatchingProducer & operator = (SingleCRTMatchingProducer const &) = delete;
  SingleCRTMatchingProducer & operator = (SingleCRTMatchingProducer &&) = delete;
  
  bool moduleMatcher(int module1, int module2);
  void produce(art::Event & event) override;
  std::string fTrackModuleLabel = "pandoraTrack";

  int nEvents = 0;
  int nHaloMuons = 0;
  private: ofstream logFile;
  const std::string reco="reco";
  //Parameters for reading in CRT::Triggers and associated AuxDetSimChannels.
  art::InputTag fCRTLabel; //Label for the module that produced 
    bool fMCCSwitch;
    bool fModuleSwitch;
    int fADCThreshold;
    bool fSCECorrection;
    int fModuletoModuleTimingCut;
    int fFronttoBackTimingCut;
    int fOpCRTTDiffCut;


    double dotCos;
    int adcX, adcY;
    double X_CRT, Y_CRT, Z_CRT;
    double trackX1, trackX2, trackY1, trackY2, trackZ1, trackZ2;
    int moduleX, moduleY; 
    int stripX, stripY;
    double deltaX, deltaY;
    double CRTT0;
    double flashTime;
    double opCRTTDiff;
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
    int stripX, stripY, moduleX, moduleY;
    double hitPositionX;
    double hitPositionY;
    double hitPositionZ;
    double timeAvg;
    int trigNumberX, trigNumberY;
  }
  recoHits;

  typedef struct // These are displacement metrics for track and hit reco
  {
    int tempId;
    int CRTTrackId;
    int recoId;
    int adcX1;
    int adcY1;
    double deltaX;
    double deltaY;
    double dotProductCos;
    double X1;
    double Y1;
    double Z1;
    int trackID;
    TVector3 trackStartPosition;
    TVector3 trackEndPosition;
    int moduleX1, moduleY1;
    int stripX1, stripY1;
    double flashTDiff;
    double timeAvg;
    int trigNumberX, trigNumberY; 
  }
  tracksPair;

  struct removePairIndex // Struct to remove tracks after sorting
  {
    const tracksPair tracksPair1;
    removePairIndex(const tracksPair & tracksPair0): tracksPair1(tracksPair0) {}

    bool operator()(const tracksPair & tracksPair2) {
      return (tracksPair1.recoId == tracksPair2.recoId || tracksPair1.CRTTrackId == tracksPair2.CRTTrackId || (tracksPair1.stripX1==tracksPair2.stripX1 && tracksPair1.moduleX1==tracksPair2.moduleX1) || (tracksPair1.stripY1==tracksPair2.stripY1 && tracksPair1.moduleY1==tracksPair2.moduleY1));
    }
  };

  struct sortPair // Struct to sort to find best CRT track for TPC track
  {
    bool operator()(const tracksPair & pair1,
      const tracksPair & pair2) {
     
	//return (fabs(pair1.dotProductCos)>fabs(pair2.dotProductCos));
	return (fabs(pair1.deltaX)+fabs(pair1.deltaY)<fabs(pair2.deltaX)+fabs(pair2.deltaY));

  }
  };
  std::vector < recoHits > primaryHits_F;
  std::vector < recoHits > primaryHits_B;

  std::vector < tempHits > tempHits_F;
  std::vector < tempHits > tempHits_B;
  std::vector < tracksPair > tracksPair_F;
  std::vector < tracksPair > tracksPair_B;
};

CRT::SingleCRTMatchingProducer::SingleCRTMatchingProducer(fhicl::ParameterSet
    const & p):
  EDProducer{p}, fCRTLabel(p.get < art::InputTag > ("CRTLabel")) {
    consumes < std::vector < CRT::Trigger >> (fCRTLabel);
    consumes < std::vector < art::Assns < sim::AuxDetSimChannel, CRT::Trigger >>> (fCRTLabel); // CRT art consumables
  fMCCSwitch=(p.get<bool>("MCC"));
  produces< std::vector<anab::T0> >();
  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Track, anab::T0> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();
  produces< art::Assns<anab::CosmicTag, anab::T0> >();
  produces< art::Assns<CRT::Trigger, anab::CosmicTag> >();
  fSCECorrection=(p.get<bool>("SCECorrection"));
  }



// v6 Geo Channel Map
bool CRT::SingleCRTMatchingProducer::moduleMatcher(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 6 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 26 || module2 == 27)) || (module1 == 31 && (module2 == 26 || module2 == 27)) || (module1 == 7 && (module2 == 12 || module2 == 13)) || (module1 == 15 && (module2 == 12 || module2 == 13)) || (module1 == 18 && (module2 == 24 || module2 == 25)) || (module1 == 30 && (module2 == 24 || module2 == 25)) || (module1 == 1 && (module2 == 4 || module2 == 5)) || (module1 == 9 && (module2 == 4 || module2 == 5)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 0 && (module2 == 2 || module2 == 3)) || (module1 == 8 && (module2 == 2 || module2 == 3)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23))) return 1;
  else return 0;

}


//Turn sim::AuxDetSimChannels into CRT::Hits. 
void CRT::SingleCRTMatchingProducer::produce(art::Event & event)
{

  std::unique_ptr< std::vector<anab::T0> > T0col( new std::vector<anab::T0>);
 
  auto CRTTrack=std::make_unique< std::vector< anab::CosmicTag > > ();

 std::unique_ptr< art::Assns<anab::CosmicTag, anab::T0> > CRTT0assn( new art::Assns<anab::CosmicTag, anab::T0>);

 std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag> > TPCCRTassn( new art::Assns<recob::Track, anab::CosmicTag>);
 std::unique_ptr< art::Assns<recob::Track, anab::T0> > TPCT0assn( new art::Assns<recob::Track, anab::T0>);

 std::unique_ptr< art::Assns<CRT::Trigger, anab::CosmicTag>> CRTTriggerassn( new art::Assns<CRT::Trigger, anab::CosmicTag>);

  if (fMCCSwitch){
    fModuleSwitch=1;
    fADCThreshold=800;
    fModuletoModuleTimingCut=4;
    fFronttoBackTimingCut=100;
    fOpCRTTDiffCut=200;
    
}
  else {
    fModuleSwitch=0;
    fADCThreshold=10;
    fModuletoModuleTimingCut=5;
    fFronttoBackTimingCut=8;
    fOpCRTTDiffCut=1000000;


}

/*   if(!fMCCSwitch){
   art::ValidHandle<std::vector<raw::RDTimeStamp>> timingHandle = event.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");
   //const auto& pdspctbs = event.getValidHandle<std::vector<raw::ctb::pdspctb>>(fCTB_tag);

	//const raw::RDTimeStamp& timeStamp = timingHandle->at(0);
	//if(timeStamp.GetFlags()!= 13) {event.put(std::move(CRTTrack)); event.put(std::move(T0col)); event.put(std::move(TPCCRTassn)); event.put(std::move(CRTT0assn));  event.put(std::move(TPCT0assn)); return;}
 }*/

  int nHits = 0;
  art::ServiceHandle < cheat::BackTrackerService > backTracker;
  art::ServiceHandle < cheat::ParticleInventoryService > partInventory;

  primaryHits_F.clear();
  primaryHits_B.clear();
  tracksPair_F.clear();
  tracksPair_B.clear();
  tempHits_F.clear();
  tempHits_B.clear(); // Arrays to compile hits and move them through
  primaryHits_F.clear();
    //allTracksPair.clear();
  logFile.open("ProtoDUNE.log"); // Logfile I don't use right now

  //Get triggers
  //cout << "Getting triggers" << endl;
  const auto & triggers = event.getValidHandle < std::vector < CRT::Trigger >> (fCRTLabel);

  art::FindManyP < sim::AuxDetSimChannel > trigToSim(triggers, event, fCRTLabel);

  vector < art::Ptr < CRT::Trigger > > crtList;
  auto crtListHandle = event.getHandle < vector < CRT::Trigger > >(fCRTLabel);
  if (crtListHandle) {
    art::fill_ptr_vector(crtList, crtListHandle);
  }

  //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
  art::ServiceHandle < geo::Geometry > geom;
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  //Mapping from channel to trigger
  std::unordered_map < size_t, double > prevTimes;
  int hitID = 0;
  int trigID=0;
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
        tHits.triggerNumber=trigID;
	}
	else{
        tHits.module = trigger.Channel(); // Values to add to array
        tHits.channelGeo = hit.Channel();
	tHits.channel=hit.Channel();
        tHits.adc = hit.ADC();
	tHits.triggerTime=trigger.Timestamp();
        tHits.triggerNumber=trigID;
	}
	 //cout<<trigger.Channel()<<','<<hit.Channel()<<','<<hit.ADC()<<endl;
        nHits++;

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

  //cout << "Hits compiled for event: " << nEvents << endl;
  //cout << "Number of Hits above Threshold:  " << hitID << endl;

  for (unsigned int f = 0; f < tempHits_F.size(); f++) {
    for (unsigned int f_test = 0; f_test < tempHits_F.size(); f_test++) {
       if (fabs(tempHits_F[f_test].triggerTime-tempHits_F[f].triggerTime)>fModuletoModuleTimingCut) continue;
      const auto & trigGeo = geom -> AuxDet(tempHits_F[f].module);
      const auto & trigGeo2 = geom -> AuxDet(tempHits_F[f_test].module);

      const auto & hit1Geo = trigGeo.SensitiveVolume(tempHits_F[f].channelGeo);
      const auto hit1Center = hit1Geo.GetCenter();
      // Create 2D hits from geo of the Y and X modules
       const auto & hit2Geo = trigGeo2.SensitiveVolume(tempHits_F[f_test].channelGeo);
      const auto hit2Center = hit2Geo.GetCenter();
      bool moduleMatched;
      moduleMatched=moduleMatcher(tempHits_F[f_test].module, tempHits_F[f].module);

      if (moduleMatched) {
        // Get the center of the hits (CRT_Res=2.5 cm)
        double hitX = hit1Center.X();
	for (unsigned int a = 0; a < tempHits_F.size(); a++)
	{
	if(tempHits_F[a].module==tempHits_F[f].module && (tempHits_F[a].channelGeo-1)==tempHits_F[f].channelGeo) hitX=hit1Center.X()+1.25;
	}
	double hitYPrelim=hit2Center.Y();
	for (unsigned int a = 0; a < tempHits_F.size(); a++)
	{
	if(tempHits_F[a].module==tempHits_F[f_test].module && (tempHits_F[a].channelGeo-1)==tempHits_F[f_test].channelGeo) hitYPrelim=hit2Center.Y()+1.25;
	}
	

	
	double hitY=hitYPrelim;
        double hitZ = (hit1Center.Z() + hit2Center.Z()) / 2.f;

        recoHits rHits;
        rHits.adcX=tempHits_F[f].adc;
	rHits.adcY=tempHits_F[f_test].adc;
        rHits.hitPositionX = hitX;
        rHits.hitPositionY = hitY;
        rHits.hitPositionZ = hitZ;
	rHits.moduleX=tempHits_F[f].module;
        rHits.moduleY=tempHits_F[f_test].module;
        rHits.trigNumberX=tempHits_F[f].triggerNumber;
	rHits.trigNumberY=tempHits_F[f_test].triggerNumber;
	rHits.stripX=tempHits_F[f].channel;
	rHits.stripY=tempHits_F[f_test].channel;
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
      const auto & hit1Geo = trigGeo.SensitiveVolume(tempHits_B[f].channelGeo);
      const auto hit1Center = hit1Geo.GetCenter();

      const auto & hit2Geo = trigGeo2.SensitiveVolume(tempHits_B[f_test].channelGeo);
      const auto hit2Center = hit2Geo.GetCenter();
      bool moduleMatched;
      moduleMatched=moduleMatcher(tempHits_B[f_test].module, tempHits_B[f].module);


      if (moduleMatched) {
        double hitX = hit1Center.X();
	
	
	for (unsigned int a = 0; a < tempHits_B.size(); a++)
	{
	if(tempHits_B[a].module==tempHits_B[f].module && (tempHits_B[a].channelGeo-1)==tempHits_B[f].channelGeo) hitX=hit1Center.X()+1.25;
	}
	
        double hitYPrelim = hit2Center.Y();
	
	for (unsigned int a = 0; a < tempHits_B.size(); a++)
	{
	if(tempHits_B[a].module==tempHits_B[f_test].module && (tempHits_B[a].channel-1)==tempHits_B[f_test].channel) hitYPrelim=hit2Center.Y()+1.25;
	}
	double hitY=hitYPrelim;

	
        double hitZ = (hit1Center.Z() + hit2Center.Z()) / 2.f;

        recoHits rHits;
        rHits.adcX=tempHits_B[f].adc;
	rHits.adcY=tempHits_B[f_test].adc;
        rHits.hitPositionX = hitX;
        rHits.hitPositionY = hitY;
        rHits.hitPositionZ = hitZ;
        rHits.moduleX=tempHits_B[f].module;
        rHits.moduleY=tempHits_B[f_test].module;
	rHits.stripX=tempHits_B[f].channel;
	rHits.trigNumberX=tempHits_B[f].triggerNumber;
	rHits.trigNumberY=tempHits_B[f_test].triggerNumber;
	rHits.stripY=tempHits_B[f_test].channel;
	rHits.timeAvg = (tempHits_B[f_test].triggerTime+tempHits_B[f].triggerTime)/2.0;
	 primaryHits_B.push_back(rHits); 

	 }
    }
  }

     auto const t0CandPtr = art::PtrMaker<anab::T0>(event);
     auto const crtPtr = art::PtrMaker<anab::CosmicTag>(event);	
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

  int nTracksReco = trackList.size();
  art::FindManyP < recob::Hit > hitsFromTrack(trackListHandle, event, fTrackModuleLabel);
  int tempId = 0;

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);

  for (int iRecoTrack = 0; iRecoTrack < nTracksReco; ++iRecoTrack) {
    if (primaryHits_F.size()+primaryHits_B.size()<1) break;
   std::vector< art::Ptr<recob::Hit> > allHits =  hitsFromTrack.at(iRecoTrack);

      art::Ptr<recob::Track> ptrack(trackListHandle, iRecoTrack);
     
	std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(iRecoTrack);
	if(!pfps.size()) continue;
	std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
	if(t0s.size()){ 
	  //auto t0=t0s.at(0);
	  //int t_zero=t0->Time();
	  //cout<<"Pandora T0: "<<t_zero<<endl;
   	}
    int firstHit=0;
    int lastHit=allHits.size()-2;


// Get track positions and find angles
 //  if (fMCCSwith==1){
    double trackStartPositionZ_noSCE = trackList[iRecoTrack]->Vertex().Z();
    double trackEndPositionZ_noSCE = trackList[iRecoTrack] -> End().Z();

    double trackStartPositionX_noSCE = trackList[iRecoTrack]->Vertex().X();
    double trackStartPositionY_noSCE = trackList[iRecoTrack]->Vertex().Y();


    double trackEndPositionX_noSCE = trackList[iRecoTrack] -> End().X();
    double trackEndPositionY_noSCE = trackList[iRecoTrack] -> End().Y();

   if (trackStartPositionZ_noSCE>trackEndPositionZ_noSCE){
    trackEndPositionZ_noSCE = trackList[iRecoTrack]->Vertex().Z();
    trackStartPositionZ_noSCE = trackList[iRecoTrack] -> End().Z();
    trackEndPositionX_noSCE = trackList[iRecoTrack]->Vertex().X();
    trackEndPositionY_noSCE = trackList[iRecoTrack]->Vertex().Y();


    trackStartPositionX_noSCE = trackList[iRecoTrack] -> End().X();
    trackStartPositionY_noSCE = trackList[iRecoTrack] -> End().Y();
    firstHit=lastHit;
    lastHit=0;
    }
 if ((trackEndPositionZ_noSCE>90 && trackEndPositionZ_noSCE < 660 && trackStartPositionZ_noSCE <50 && trackStartPositionZ_noSCE<660) || (trackStartPositionZ_noSCE>90 && trackStartPositionZ_noSCE < 660 && trackEndPositionZ_noSCE <50 && trackEndPositionZ_noSCE<660)) {

      double min_delta = DBL_MAX;

      double best_dotProductCos = DBL_MAX;
      double best_deltaXF = DBL_MAX;
      double best_deltaYF = DBL_MAX;
      int bestHitIndex_F=0;
      double best_trackX1=DBL_MAX;
      double best_trackX2=DBL_MAX;
      for (unsigned int iHit_F = 0; iHit_F < primaryHits_F.size(); iHit_F++) {
   double xOffset=0;

		double trackStartPositionX_notCorrected=trackStartPositionX_noSCE;
		double trackEndPositionX_notCorrected=trackEndPositionX_noSCE;
		if (!t0s.empty()){
		if (event.isRealData() && fabs(t0s.at(0)->Time()-(primaryHits_F[iHit_F].timeAvg*20.f))>100000) continue;
		if (!event.isRealData() && fabs(t0s.at(0)->Time()-primaryHits_F[iHit_F].timeAvg)>100000) continue;
		}
		if (t0s.empty()){
		int RDOffset=0;
		if (!fMCCSwitch) RDOffset=111;
		double ticksOffset=0;
		//cout<<(primaryHits_F[iHit_F].timeAvg+RDOffset)<<endl;
                //cout<<detProp.GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat)<<endl;
                if (!fMCCSwitch) ticksOffset = (primaryHits_F[iHit_F].timeAvg+RDOffset)/25.f+detProp.GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);

                else if (fMCCSwitch) ticksOffset = (primaryHits_F[iHit_F].timeAvg/500.f)+detProp.GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		
               xOffset=detProp.ConvertTicksToX(ticksOffset,allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		//double xOffset=.08*ticksOffset
		
	 trackStartPositionX_noSCE=trackStartPositionX_notCorrected-xOffset;
         trackEndPositionX_noSCE=trackEndPositionX_notCorrected-xOffset;
	if (fabs(xOffset)>300 || ((trackStartPositionX_notCorrected<0 && trackStartPositionX_noSCE>0) || (trackEndPositionX_notCorrected<0 && trackEndPositionX_noSCE>0)) || ((trackStartPositionX_notCorrected>0 && trackStartPositionX_noSCE<0) || (trackEndPositionX_notCorrected>0 && trackEndPositionX_noSCE<0)) ) continue;
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

        double X1 = primaryHits_F[iHit_F].hitPositionX;

        double Y1 = primaryHits_F[iHit_F].hitPositionY;

        double Z1 = primaryHits_F[iHit_F].hitPositionZ;

	// Make metrics for a CRT pair to compare later
	TVector3 trackStart(trackStartPositionX, trackStartPositionY, trackStartPositionZ);
	TVector3 trackEnd(trackEndPositionX, trackEndPositionY, trackEndPositionZ);
	TVector3 v1(X1,Y1,Z1);
	TVector3 v2(trackStartPositionX, trackStartPositionY, trackStartPositionZ);

            TVector3 v4(trackStartPositionX,
                        trackStartPositionY,
                        trackStartPositionZ);
            TVector3 v5(trackEndPositionX,
                        trackEndPositionY,
                        trackEndPositionZ);
	TVector3 trackVector = (v5-v4).Unit();
	TVector3 hitVector=(v2-v1).Unit();




              double predictedHitPositionY1 = (v1.Z()-v5.Z())/(v4.Z()-v5.Z())*(v4.Y()-v5.Y())+v5.Y();


              double predictedHitPositionX1 = (v1.Z()-v5.Z())/(v4.Z()-v5.Z())*(v4.X()-v5.X())+v5.X();


	if (predictedHitPositionX1<-200 || predictedHitPositionX1>580 || predictedHitPositionY1<-50 || predictedHitPositionY1>620) continue;

	double dotProductCos=trackVector*hitVector;

        double deltaX1 = (predictedHitPositionX1-X1);


        double deltaY1 = (predictedHitPositionY1-Y1);


   if (min_delta > std::abs(deltaX1) + std::abs(deltaY1) ){

	   min_delta=std::abs(deltaX1)+std::abs(deltaY1);

            best_dotProductCos = dotProductCos;
	    best_trackX1=trackStartPositionX;
	    best_trackX2=trackEndPositionX;
            best_deltaXF = deltaX1;
            best_deltaYF = deltaY1;
	    bestHitIndex_F=iHit_F;
          }

	}
   int  f=bestHitIndex_F; 
 
   double deltaX=best_deltaXF; double deltaY=best_deltaYF;
   double dotProductCos=best_dotProductCos;

   double trackStartPositionX=best_trackX1;
   double trackStartPositionY=trackStartPositionY_noSCE;
   double trackStartPositionZ=trackStartPositionZ_noSCE;

   double trackEndPositionX=best_trackX2;
   double trackEndPositionY=trackEndPositionY_noSCE;
   double trackEndPositionZ=trackEndPositionZ_noSCE;


	TVector3 trackStart(trackStartPositionX, trackStartPositionY, trackStartPositionZ);
	TVector3 trackEnd(trackEndPositionX, trackEndPositionY, trackEndPositionZ);

         tracksPair tPair;
        tPair.tempId = tempId;
        tPair.CRTTrackId = f;
        tPair.recoId = iRecoTrack;

	tPair.deltaX=deltaX;
        tPair.deltaY=deltaY;
        tPair.dotProductCos=dotProductCos;

        tPair.moduleX1 = primaryHits_F[f].moduleX;
        tPair.moduleY1 = primaryHits_F[f].moduleY;

        tPair.adcX1=primaryHits_F[f].adcX;
        tPair.adcY1=primaryHits_F[f].adcY;

        tPair.stripX1 = primaryHits_F[f].stripX;
        tPair.stripY1 = primaryHits_F[f].stripY;
        tPair.trigNumberX = primaryHits_F[f].trigNumberX;
        tPair.trigNumberY = primaryHits_F[f].trigNumberY;
        tPair.X1 = primaryHits_F[f].hitPositionX;
        tPair.Y1 = primaryHits_F[f].hitPositionY;
        tPair.Z1 = primaryHits_F[f].hitPositionZ;
	tPair.timeAvg=primaryHits_F[f].timeAvg;
        tPair.trackStartPosition=trackStart;
	tPair.trackEndPosition=trackEnd;
        tracksPair_B.push_back(tPair);

      }
	
 if ( (trackStartPositionZ_noSCE<620 && trackEndPositionZ_noSCE > 660 && trackStartPositionZ_noSCE > 50 && trackEndPositionZ_noSCE > 50) || (trackStartPositionZ_noSCE>660 && trackEndPositionZ_noSCE < 620 && trackStartPositionZ_noSCE > 50 && trackEndPositionZ_noSCE > 50)) {
      double min_delta = DBL_MAX;

      double best_dotProductCos = DBL_MAX;
      double best_deltaXF = DBL_MAX;
      double best_deltaYF = DBL_MAX;     
      int bestHitIndex_B=0;
      double best_trackX1=DBL_MAX;
      double best_trackX2=DBL_MAX;
      for (unsigned int iHit_B = 0; iHit_B < primaryHits_B.size(); iHit_B++) {
double xOffset=0;
    

		double trackStartPositionX_notCorrected=trackStartPositionX_noSCE;
		double trackEndPositionX_notCorrected=trackEndPositionX_noSCE;
		if (!t0s.empty()){
		if (event.isRealData() && fabs(t0s.at(0)->Time()-(primaryHits_B[iHit_B].timeAvg*20.f))>100000) continue;
		if (!event.isRealData() && fabs(t0s.at(0)->Time()-primaryHits_B[iHit_B].timeAvg)>100000) continue;
	}
		if (t0s.empty()){
		int RDOffset=0;
		if (!fMCCSwitch) RDOffset=111;
		double ticksOffset=0;
                if (!fMCCSwitch) ticksOffset = (primaryHits_B[iHit_B].timeAvg+RDOffset)/25.f+detProp.GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);

                else if (fMCCSwitch) ticksOffset = (primaryHits_B[iHit_B].timeAvg/500.f)+detProp.GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		
               xOffset=detProp.ConvertTicksToX(ticksOffset,allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
		//double xOffset=.08*ticksOffset
		
	 trackStartPositionX_noSCE=trackStartPositionX_notCorrected-xOffset;
         trackEndPositionX_noSCE=trackEndPositionX_notCorrected-xOffset;

	if (fabs(xOffset)>300 || ((trackStartPositionX_notCorrected<0 && trackStartPositionX_noSCE>0) || (trackEndPositionX_notCorrected<0 && trackEndPositionX_noSCE>0)) || ((trackStartPositionX_notCorrected>0 && trackStartPositionX_noSCE<0) || (trackEndPositionX_notCorrected>0 && trackEndPositionX_noSCE<0)) ) continue;
	}

   double trackStartPositionX=trackStartPositionX_noSCE;
   double trackStartPositionY=trackStartPositionY_noSCE;
   double trackStartPositionZ=trackStartPositionZ_noSCE;

   double trackEndPositionX=trackEndPositionX_noSCE;
   double trackEndPositionY=trackEndPositionY_noSCE;
   double trackEndPositionZ=trackEndPositionZ_noSCE;
        double X1 = primaryHits_B[iHit_B].hitPositionX;

        double Y1 = primaryHits_B[iHit_B].hitPositionY;

        double Z1 = primaryHits_B[iHit_B].hitPositionZ;


 
	// Make metrics for a CRT pair to compare later
	TVector3 trackStart(trackStartPositionX, trackStartPositionY, trackStartPositionZ);
	TVector3 trackEnd(trackEndPositionX, trackEndPositionY, trackEndPositionZ);
	TVector3 v1(X1,Y1,Z1);
	TVector3 v2(trackStartPositionX, trackStartPositionY, trackStartPositionZ);

            TVector3 v4(trackStartPositionX,
                        trackStartPositionY,
                        trackStartPositionZ);
            TVector3 v5(trackEndPositionX,
                        trackEndPositionY,
                        trackEndPositionZ);
	TVector3 trackVector = (v5-v4).Unit();
	TVector3 hitVector=(v2-v1).Unit();




              double predictedHitPositionY1 = (v1.Z()-v5.Z())/(v4.Z()-v5.Z())*(v4.Y()-v5.Y())+v5.Y();


              double predictedHitPositionX1 = (v1.Z()-v5.Z())/(v4.Z()-v5.Z())*(v4.X()-v5.X())+v5.X();

	if (abs(predictedHitPositionX1)>340 || predictedHitPositionY1<-160 || predictedHitPositionY1>560) continue;

	double dotProductCos=trackVector*hitVector;

        double deltaX1 = (predictedHitPositionX1-X1);

        double deltaY1 = (predictedHitPositionY1-Y1);

      if (min_delta > std::abs(deltaX1) + std::abs(deltaY1) ){

	   min_delta=std::abs(deltaX1)+std::abs(deltaY1);

            best_dotProductCos = dotProductCos;
	    best_trackX1=trackStartPositionX;
	    best_trackX2=trackEndPositionX;
            best_deltaXF = deltaX1;
            best_deltaYF = deltaY1;
	    bestHitIndex_B=iHit_B; 
          }

	}
   int  f=bestHitIndex_B; 
 
   double deltaX=best_deltaXF; double deltaY=best_deltaYF;
   double dotProductCos=best_dotProductCos;


   double trackStartPositionX=best_trackX1;
   double trackStartPositionY=trackStartPositionY_noSCE;
   double trackStartPositionZ=trackStartPositionZ_noSCE;

   double trackEndPositionX=best_trackX2;
   double trackEndPositionY=trackEndPositionY_noSCE;
   double trackEndPositionZ=trackEndPositionZ_noSCE;

	TVector3 trackStart(trackStartPositionX, trackStartPositionY, trackStartPositionZ);
	TVector3 trackEnd(trackEndPositionX, trackEndPositionY, trackEndPositionZ);


        tracksPair tPair;
        tPair.tempId = tempId;
        tPair.CRTTrackId = f;
        tPair.recoId = iRecoTrack;

	tPair.deltaX=deltaX;
        tPair.deltaY=deltaY;
        tPair.dotProductCos=dotProductCos;

        tPair.moduleX1 = primaryHits_B[f].moduleX;
        tPair.moduleY1 = primaryHits_B[f].moduleY;

        tPair.adcX1=primaryHits_B[f].adcX;
        tPair.adcY1=primaryHits_B[f].adcY;

        tPair.stripX1 = primaryHits_B[f].stripX;
        tPair.stripY1 = primaryHits_B[f].stripY;
       tPair.trigNumberX = primaryHits_B[f].trigNumberX;
        tPair.trigNumberY = primaryHits_B[f].trigNumberY;
        tPair.X1 = primaryHits_B[f].hitPositionX;
        tPair.Y1 = primaryHits_B[f].hitPositionY;
        tPair.Z1 = primaryHits_B[f].hitPositionZ;
	tPair.timeAvg=primaryHits_B[f].timeAvg;
        tPair.trackStartPosition=trackStart;
	tPair.trackEndPosition=trackEnd;
        tracksPair_B.push_back(tPair);


      tempId++;
    } //iRecoTrack
    }

     //Sort pair by ascending order of absolute distance
    sort(tracksPair_F.begin(), tracksPair_F.end(), sortPair());
    sort(tracksPair_B.begin(), tracksPair_B.end(), sortPair());


    vector < tracksPair > allUniqueTracksPair;
    while (tracksPair_F.size()) {
      allUniqueTracksPair.push_back(tracksPair_F.front());
      tracksPair_F.erase(remove_if(tracksPair_F.begin(), tracksPair_F.end(), removePairIndex(tracksPair_F.front())),
        tracksPair_F.end());
    }

    while (tracksPair_B.size()) {
      allUniqueTracksPair.push_back(tracksPair_B.front());
      tracksPair_B.erase(remove_if(tracksPair_B.begin(), tracksPair_B.end(), removePairIndex(tracksPair_B.front())),
        tracksPair_B.end());
    }

	//cout<<"Number of reco and CRT pairs: "<<allUniqueTracksPair.size()<<endl;
    if (allUniqueTracksPair.size() > 0) {
      for (unsigned int u = 0; u < allUniqueTracksPair.size(); u++) {

	int CRTTrackId=allUniqueTracksPair[u].CRTTrackId;
	int TPCTrackId=allUniqueTracksPair[u].recoId;
	deltaX=allUniqueTracksPair[u].deltaX;

	deltaY=allUniqueTracksPair[u].deltaY;
	opCRTTDiff=allUniqueTracksPair[u].flashTDiff;

	dotCos=fabs(allUniqueTracksPair[u].dotProductCos);
	trackX1=allUniqueTracksPair[u].trackStartPosition.X();
	trackY1=allUniqueTracksPair[u].trackStartPosition.Y();
	trackZ1=allUniqueTracksPair[u].trackStartPosition.Z();

	trackX2=allUniqueTracksPair[u].trackEndPosition.X();
	trackY2=allUniqueTracksPair[u].trackEndPosition.Y();
	trackZ2=allUniqueTracksPair[u].trackEndPosition.Z();

	moduleX=allUniqueTracksPair[u].moduleX1;
	moduleY=allUniqueTracksPair[u].moduleY1;

	adcX=allUniqueTracksPair[u].adcX1;
	adcY=allUniqueTracksPair[u].adcY1;

	CRTT0=allUniqueTracksPair[u].timeAvg;

	if (!fMCCSwitch) CRTT0=(111.f+allUniqueTracksPair[u].timeAvg)*20.f;
	// Added 111 tick CRT-CTB Timing Offset
	stripX=allUniqueTracksPair[u].stripX1;
	stripY=allUniqueTracksPair[u].stripY1;

	X_CRT=allUniqueTracksPair[u].X1;
	Y_CRT=allUniqueTracksPair[u].Y1;
	Z_CRT=allUniqueTracksPair[u].Z1;

       	flashTime=-1*opCRTTDiff-CRTT0;
        if ( fabs(trackX1)<400 &&  fabs(trackX2)<400 && fabs(deltaX)<60 &&  fabs(deltaY)<60 && dotCos>0.9995 ) {
	//cout<<"Found Matched Single CRT Tag with CRT*TPC: "<<fabs(allUniqueTracksPair[u].dotProductCos)<<endl;
	//cout<<"Displacement of match:"<<deltaX<<','<<deltaY<<endl;

	std::vector<float> hitF;
	std::vector<float> hitB;
	hitF.push_back(X_CRT); hitF.push_back(Y_CRT); hitF.push_back(Z_CRT);
	hitB.push_back(trackX1); hitB.push_back(trackY1); hitB.push_back(trackZ1);
	CRTTrack->push_back(anab::CosmicTag(hitF,hitB, fabs(allUniqueTracksPair[u].dotProductCos),anab::CosmicTagID_t::kUnknown));
	if (Z_CRT<100) T0col->push_back(anab::T0(CRTT0, 1, 1,CRTTrackId,fabs(allUniqueTracksPair[u].dotProductCos)));
        else T0col->push_back(anab::T0(CRTT0, 2, 1,CRTTrackId,fabs(allUniqueTracksPair[u].dotProductCos) ));
	auto const crtTrackPtr = crtPtr(CRTTrack->size()-1);
	auto const t0CP = t0CandPtr(CRTTrackId);
	CRTT0assn->addSingle(crtTrackPtr,t0CP);
	
	util::CreateAssn(*this, event, *T0col, trackList[TPCTrackId], *TPCT0assn);
	util::CreateAssn(*this, event, *CRTTrack, trackList[TPCTrackId], *TPCCRTassn);
	util::CreateAssn(*this, event, *CRTTrack, crtList[allUniqueTracksPair[u].trigNumberX], *CRTTriggerassn);
	util::CreateAssn(*this, event, *CRTTrack, crtList[allUniqueTracksPair[u].trigNumberY], *CRTTriggerassn);

        }
      }
    }
  nEvents++;
event.put(std::move(T0col)); event.put(std::move(CRTTrack)); event.put(std::move(TPCCRTassn));  event.put(std::move(TPCT0assn)); event.put(std::move(CRTT0assn)); event.put(std::move(CRTTriggerassn));

}




DEFINE_ART_MODULE(CRT::SingleCRTMatchingProducer)
