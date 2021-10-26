////////////////////////////////////////////////////////////////////////
// Class:       TwoCRTMatchingProducer
// Plugin Type: producer (art v2_10_03)
// File:        TwoCRTMatchingProducer_module.cc
// Author: Richie Diurba (rdiurba@fnal.gov)
//	   Tingjun Yang (tjyang@fnal.gov)
//
// Generated at Wed Jun 27 04:09:39 2018 by Andrew Olivier using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "lardata/Utilities/AssociationUtil.h"

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
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"


#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

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
using namespace std;   // Namespaces established to make life easier
using namespace ROOT::Math;

namespace CRT {
  class TwoCRTMatchingProducer;
}

class CRT::TwoCRTMatchingProducer : public art::EDProducer {
public:
  explicit TwoCRTMatchingProducer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TwoCRTMatchingProducer(TwoCRTMatchingProducer const &) = delete;
  TwoCRTMatchingProducer(TwoCRTMatchingProducer &&) = delete;
  TwoCRTMatchingProducer & operator = (TwoCRTMatchingProducer const &) = delete;
  TwoCRTMatchingProducer & operator = (TwoCRTMatchingProducer &&) = delete;
  bool moduleMatcher(int module1, int module2);
  // Required functions.

  void produce(art::Event & event) override;

  int nEvents = 0;
  int nHitsPerEvent=0;
  std::string fTrackModuleLabel = "pandoraTrack";


private:
    art::InputTag fCRTLabel; //Label for the module that analyzed 
    bool fSCECorrection;
    bool fModuleSwitch;
    int fADCThreshold;
    int fModuletoModuleTimingCut;
    int fFronttoBackTimingCut;
    long long timeStamp;


typedef struct // Structures for arrays to move hits from raw to reco to validation
  {
    int triggerNumber;
    int channel;
    int module;
    int adc;
    int triggerTime;
  }
  tempHits;

  typedef struct {
    int tempId;	
    int trigNumberX;
    int trigNumberY;
    double hitPositionX;
    double hitPositionY;
    double hitPositionZ;
    double timeAvg;
  }
  recoHits;


  std::vector < recoHits > primaryHits_F;
  std::vector < recoHits > primaryHits_B;

  std::vector < tempHits > tempHits_F;
  std::vector < tempHits > tempHits_B;



};


CRT::TwoCRTMatchingProducer::TwoCRTMatchingProducer(fhicl::ParameterSet const & p): EDProducer{p}, fCRTLabel(p.get < art::InputTag > ("CRTLabel"))
{
  consumes < std::vector < CRT::Trigger >> (fCRTLabel);
  consumes < std::vector < art::Assns < sim::AuxDetSimChannel, CRT::Trigger >>> (fCRTLabel);
  produces< std::vector<anab::T0> >();

  produces< std::vector<anab::CosmicTag> >();  
  produces< art::Assns<recob::Track, anab::T0> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();
  produces< art::Assns<CRT::Trigger, anab::CosmicTag> >();

  fSCECorrection=(p.get<bool>("SCECorrection"));
}

// v6 Geo Channel Map
bool CRT::TwoCRTMatchingProducer::moduleMatcher(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 6 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 26 || module2 == 27)) || (module1 == 31 && (module2 == 26 || module2 == 27)) || (module1 == 7 && (module2 == 12 || module2 == 13)) || (module1 == 15 && (module2 == 12 || module2 == 13)) || (module1 == 18 && (module2 == 24 || module2 == 25)) || (module1 == 30 && (module2 == 24 || module2 == 25)) || (module1 == 1 && (module2 == 4 || module2 == 5)) || (module1 == 9 && (module2 == 4 || module2 == 5)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 0 && (module2 == 2 || module2 == 3)) || (module1 == 8 && (module2 == 2 || module2 == 3)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23))) return 1;
  else return 0;

}


//Turn sim::AuxDetSimChannels into CRT::Hits. 
void CRT::TwoCRTMatchingProducer::produce(art::Event & event)
{
    bool fMCCSwitch=!event.isRealData();
    nEvents++;	
  std::unique_ptr< std::vector<anab::T0> > T0col( new std::vector<anab::T0>);
  auto CRTTrack=std::make_unique< std::vector< anab::CosmicTag > > (); 


 std::unique_ptr< art::Assns<CRT::Trigger, anab::CosmicTag>> CRTTriggerassn( new art::Assns<CRT::Trigger, anab::CosmicTag>);

 std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag> > TPCCRTassn( new art::Assns<recob::Track, anab::CosmicTag>);
 std::unique_ptr< art::Assns<recob::Track, anab::T0> > TPCT0assn( new art::Assns<recob::Track, anab::T0>);
   

  if (fMCCSwitch){
    fModuleSwitch=1;
    fADCThreshold=800;
    fModuletoModuleTimingCut=2;
    fFronttoBackTimingCut=100;
    
}
  else {
    fModuleSwitch=0;
    fADCThreshold=10;
    fModuletoModuleTimingCut=5;
    fFronttoBackTimingCut=8;
    art::ValidHandle<std::vector<raw::RDTimeStamp>> timingHandle = event.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");
    timeStamp=timingHandle->at(0).GetTimeStamp();
} 
  int nHits = 0;

	//Detector properties service
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);

  primaryHits_F.clear();
  primaryHits_B.clear();
  tempHits_F.clear();
  tempHits_B.clear(); // Arrays to compile hits and move them through


  //Get triggers
  //cout << "Getting triggers" << endl;
  vector < art::Ptr < CRT::Trigger > > crtList;
  auto crtListHandle = event.getHandle < vector < CRT::Trigger > >(fCRTLabel);
  if (crtListHandle) {
    art::fill_ptr_vector(crtList, crtListHandle);
  }
  const auto & triggers = event.getValidHandle < std::vector < CRT::Trigger >> (fCRTLabel);

  art::FindManyP < sim::AuxDetSimChannel > trigToSim(triggers, event, fCRTLabel);

  //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
  art::ServiceHandle < geo::Geometry > geom;

  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  //Mapping from channel to trigger
  std::unordered_map < size_t, double > prevTimes;
  int hitID = 0;
  //cout << "Looking for hits in Triggers" << endl;
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
  //cout << "Hits compiled for event: " << nEvents << endl;
  //cout << "Number of Hits above Threshold:  " << hitID << endl;

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
	rHits.trigNumberX=tempHits_B[f].triggerNumber;
	rHits.trigNumberY=tempHits_B[f_test].triggerNumber;
	rHits.timeAvg = (tempHits_B[f_test].triggerTime+tempHits_B[f].triggerTime)/2.0;
	 primaryHits_B.push_back(rHits); 

	 }
    }
  }
  vector < art::Ptr < recob::Track > > trackList;
  auto trackListHandle = event.getHandle < vector < recob::Track > >(fTrackModuleLabel);
  if (trackListHandle) {
    art::fill_ptr_vector(trackList, trackListHandle);
  }
  else{
    event.put(std::move(T0col)); event.put(std::move(CRTTrack)); event.put(std::move(TPCCRTassn));  event.put(std::move(TPCT0assn)); event.put(std::move(CRTTriggerassn));
    return;
  }

  vector<art::Ptr<recob::PFParticle> > pfplist;
  auto PFPListHandle = event.getHandle< std::vector<recob::PFParticle> >("pandora");
  if (PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);
  if(pfplist.size()<1) return;
  art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, event ,"pandora");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,event,"pandoraTrack");
  int nTracksReco = trackList.size();
  art::FindManyP < recob::Hit > hitsFromTrack(trackListHandle, event, fTrackModuleLabel);

  for (int iRecoTrack = 0; iRecoTrack < nTracksReco; ++iRecoTrack) {
    std::vector< art::Ptr<recob::Hit> > allHits =  hitsFromTrack.at(iRecoTrack);

      art::Ptr<recob::Track> ptrack(trackListHandle, iRecoTrack);
     
	std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(iRecoTrack);
	if(!pfps.size()) continue;
	std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
    

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
      int best_trigXF=0;
      int best_trigYF=0;
      int best_trigXB=0;
      int best_trigYB=0;

    for (unsigned int f = 0; f < primaryHits_F.size(); f++) {
        for (unsigned int b = 0; b < primaryHits_B.size(); b++) {

      double X1 = primaryHits_F[f].hitPositionX;
      double Y1 = primaryHits_F[f].hitPositionY;
      double Z1 = primaryHits_F[f].hitPositionZ;
      double X2 = primaryHits_B[b].hitPositionX;
      double Y2 = primaryHits_B[b].hitPositionY;
      double Z2= primaryHits_B[b].hitPositionZ;
     

	        if (fabs(primaryHits_F[f].timeAvg-primaryHits_B[b].timeAvg)>fFronttoBackTimingCut) continue;
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
	    best_trigXF=primaryHits_F[f].trigNumberX;
	    best_trigYF=primaryHits_F[f].trigNumberY;
	    best_trigXB=primaryHits_B[b].trigNumberX;
	    best_trigYB=primaryHits_B[b].trigNumberY;
	    best_T = t0;
	    if (!fMCCSwitch) best_T=(111.f+best_T)*20.f;
	    // Added 111 tick CRT-CTB offset
          }
        }
      }
      if (std::abs(best_dotProductCos)>0.99 && std::abs(best_deltaXF)+std::abs(best_deltaXB)<40 && std::abs(best_deltaYF)+std::abs(best_deltaYB)<40 ) {
        //std::cout<<"Found match with TPC*CRT "<<best_dotProductCos<<std::endl;

//std::cout<<"Displacement in front and back "<<best_deltaXF<<","<<best_deltaYF<<","<<best_deltaXB<<","<<best_deltaYB<<std::endl;
        std::vector<float> hitF;
	std::vector<float> hitB;
	hitF.push_back(best_XF); hitF.push_back(best_YF); hitF.push_back(best_ZF);
	hitB.push_back(best_XB); hitB.push_back(best_YB); hitB.push_back(best_ZB);

        CRTTrack->push_back(anab::CosmicTag(hitF,hitB, std::abs(best_dotProductCos),anab::CosmicTagID_t::kNotTagged));

        T0col->push_back(anab::T0(best_T, 13,2,iRecoTrack,best_dotProductCos));
        util::CreateAssn(*this, event, *CRTTrack, trackList[iRecoTrack], *TPCCRTassn);
 	util::CreateAssn(*this, event, *T0col, trackList[iRecoTrack], *TPCT0assn);
        util::CreateAssn(*this, event, *CRTTrack, crtList[best_trigXF], *CRTTriggerassn);
        util::CreateAssn(*this, event, *CRTTrack, crtList[best_trigYF], *CRTTriggerassn);

        util::CreateAssn(*this, event, *CRTTrack, crtList[best_trigXB], *CRTTriggerassn);

        util::CreateAssn(*this, event, *CRTTrack, crtList[best_trigYB], *CRTTriggerassn);
	
        }
      }
    }

 event.put(std::move(T0col)); event.put(std::move(CRTTrack)); event.put(std::move(TPCCRTassn));  event.put(std::move(TPCT0assn)); event.put(std::move(CRTTriggerassn));
}




DEFINE_ART_MODULE(CRT::TwoCRTMatchingProducer)
