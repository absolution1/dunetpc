////////////////////////////////////////////////////////////////////////
// Class:       TwoCRTMatchingProducer
// Plugin Type: producer (art v3_02_06)
// File:        TwoCRTMatchingProducer_module.cc
//
// author: diurb001@umn.edu
//         tjyang@fnal.gov
//
// Generated at Thu Jul 25 22:15:30 2019 by Tingjun Yang using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/NumericUtils.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"


#include <memory>
#include <map>
#include <vector>
#include <iostream>

namespace CRT {
  class TwoCRTMatchingProducer;
}


class CRT::TwoCRTMatchingProducer : public art::EDProducer {
public:
  explicit TwoCRTMatchingProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TwoCRTMatchingProducer(TwoCRTMatchingProducer const&) = delete;
  TwoCRTMatchingProducer(TwoCRTMatchingProducer&&) = delete;
  TwoCRTMatchingProducer& operator=(TwoCRTMatchingProducer const&) = delete;
  TwoCRTMatchingProducer& operator=(TwoCRTMatchingProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fCRTLabel; //Label for the module that analyzed 
  art::InputTag fCTBLabel;
  art::InputTag fTrackModuleLabel;
  bool fMCCSwitch;
  bool fCTBTriggerOnly;
  bool fSCECorrection;
  int fADCThreshold;
  int fModuletoModuleTimingCut;
  int fFronttoBackTimingCut;
  unsigned short fMaxHitsinGroup;

  typedef struct {
    double x;
    double y;
    double z;
    unsigned long long t;
    int trigX;
    int trigY;
  }
    recoHits;

  bool moduleMatcher(int module1, int module2);
};


CRT::TwoCRTMatchingProducer::TwoCRTMatchingProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fCRTLabel(p.get < art::InputTag > ("CRTLabel")),  
  fCTBLabel(p.get<art::InputTag>("CTBLabel")),
  fTrackModuleLabel(p.get<art::InputTag>("TrackModuleLabel","pandoraTrack")),
  fMaxHitsinGroup(p.get<unsigned short>("MaxHitsinGroup",2))
  // More initializers here.
{
  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();
  produces< std::vector<anab::T0> >();
  produces< art::Assns<recob::Track, anab::T0> >();
  produces< art::Assns<CRT::Trigger, anab::CosmicTag>>();
}

void CRT::TwoCRTMatchingProducer::produce(art::Event& e)
{
  auto CRTTrack=std::make_unique< std::vector< anab::CosmicTag > > (); 
  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag> > TPCCRTassn( new art::Assns<recob::Track, anab::CosmicTag>);

  auto CRTT0=std::make_unique< std::vector< anab::T0> > (); 

  std::unique_ptr< art::Assns<recob::Track, anab::T0> > TPCT0assn( new art::Assns<recob::Track, anab::T0>);

  std::unique_ptr< art::Assns<CRT::Trigger, anab::CosmicTag>> CRTTrigassn( new art::Assns<CRT::Trigger, anab::CosmicTag>);

  fMCCSwitch = !(e.isRealData());

  //Geometry service
  art::ServiceHandle < geo::Geometry > geom;
  //Detector properties service
  auto const* detectorPropertiesService = lar::providerFrom<detinfo::DetectorPropertiesService>();

  ULong64_t rdtimestamp = 0;

  if (fMCCSwitch){
    fADCThreshold=800;
    fModuletoModuleTimingCut=4;
    fFronttoBackTimingCut=100;
  }
  else {
    fADCThreshold=10;
    fModuletoModuleTimingCut=5;
    fFronttoBackTimingCut=8;
    art::ValidHandle<std::vector<raw::RDTimeStamp>> timingHandle = e.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");
    rdtimestamp = timingHandle->at(0).GetTimeStamp();
  }

  //Get triggers
  art::Handle < std::vector < CRT::Trigger > > crtListHandle;
  std::vector < art::Ptr < CRT::Trigger > > crtList;
  if (e.getByLabel(fCRTLabel, crtListHandle)) {
    art::fill_ptr_vector(crtList, crtListHandle);
  }

  //Group adjacent CRT strips together
  //0 1 2 3 4 5 6 7 8
  //|   | |   |   | |
  //suppose those strips on a module are hit, strip 0 will in group 0, strips 2,3 will be in group 1, strip 5 will be in group 2, strips 7,8 will be in group 3, etc.
  std::vector<std::vector<std::vector<unsigned int>>> crthitsgroup(crtList.size()); //[trigger index][group index][hit index]
  for (size_t i = 0; i<crtList.size(); ++i){
    const auto& trigger = crtList[i];
    std::vector<unsigned int> crthits;
    auto & hits = trigger->Hits();
    //The map sorts the crt::Hit index based on the channel number
    std::map<unsigned short, std::vector<unsigned int>> hitchannelmap; //[hit channel, hit indices]
    for (size_t j = 0; j<hits.size(); ++j){
      if (hits[j].ADC() > fADCThreshold){
        hitchannelmap[hits[j].Channel()].push_back(j);
      }
    }
    for (const auto & [channel, indices] : hitchannelmap){
      for (const auto &index : indices){
        if (crthits.empty()){
          crthits.push_back(index);
        }
        else{
          //reach the limit on the number of hits
          if (crthits.size()==fMaxHitsinGroup){
            crthitsgroup[i].push_back(crthits);
            crthits.clear();
            crthits.push_back(index);
          }
          else{
            //The current crt::Hit is adjacent to the last crt::Hit
            if (channel == hits[crthits.back()].Channel() + 1){
              crthits.push_back(index);
              
            }
            else{
              //Start a new group
              crthitsgroup[i].push_back(crthits);
              crthits.clear();
              crthits.push_back(index);
            }
          }
        }
      }
    }
    if (!crthits.empty()) crthitsgroup[i].push_back(crthits);
  }
  
  /*
  for (size_t i = 0; i<crthitsgroup.size(); ++i){
    std::cout<<"Trigger "<<i<<std::endl;
    for (size_t j = 0; j<crthitsgroup[i].size(); ++j){
      std::cout<<"Group "<<j<<std::endl;
      for (size_t k = 0; k<crthitsgroup[i][j].size(); ++k){
        std::cout<<k<<" "<<crtList[i]->Hits()[crthitsgroup[i][j][k]].Channel()<<std::endl;
      }
    }
  }
  */
  int nmatch = 0;

  //Construct matches between X strips and Y strips
  std::vector< recoHits> hits3d_F;
  std::vector< recoHits> hits3d_B;

  for (size_t itrg = 0; itrg<crtList.size(); ++itrg){
    for (size_t jtrg = 0; jtrg<crtList.size(); ++jtrg){
      auto & xtrig = crtList[itrg];
      auto & ytrig = crtList[jtrg];
      if (!moduleMatcher(ytrig->Channel(), xtrig->Channel())) continue;
      if (int(util::absDiff(xtrig->Timestamp(), ytrig->Timestamp()))>fModuletoModuleTimingCut) continue;
      for (auto const& grpx : crthitsgroup[itrg]){
        for (auto const& grpy : crthitsgroup[jtrg]){
          recoHits hit3d;
          double hitx = 0, hity = 0, hitz0 = 0, hitz1 = 0;
          unsigned long long t0 = 0, t1 = 0;
          //Get the average x value from all crt::Hits in the group
          for (auto const& xstrip : grpx){
            auto const& trigGeo = geom->AuxDet(xtrig->Channel());
            auto const& hitGeo = trigGeo.SensitiveVolume(xtrig->Hits()[xstrip].Channel());
            auto const& hitCenter = hitGeo.GetCenter();
            hitx += hitCenter.X();
            if (!hitz0) hitz0 = hitCenter.Z();
            if (!t0) t0 = xtrig->Timestamp();
          }
          //Get the average y value from all crt::Hits in the group
          for (auto const& ystrip : grpy){
            auto const& trigGeo = geom->AuxDet(ytrig->Channel());
            auto const& hitGeo = trigGeo.SensitiveVolume(ytrig->Hits()[ystrip].Channel());
            auto const& hitCenter = hitGeo.GetCenter();
            hity += hitCenter.Y();
            if (!hitz1) hitz1 = hitCenter.Z();
            if (!t1) t1 = ytrig->Timestamp();
          }
          hit3d.x = hitx/grpx.size();
          hit3d.y = hity/grpy.size();
          hit3d.z = (hitz0+hitz1)/2;
          hit3d.t = (t0+t1)/2;
          hit3d.trigY=jtrg;
	  hit3d.trigX=itrg;
	  if (hit3d.z < 100) hits3d_F.push_back(hit3d);
          else hits3d_B.push_back(hit3d);
        }
      }
      //std::cout<<itrg<<" "<<jtrg<<" "<<crthitsgroup[itrg].size()<<" "<<crthitsgroup[jtrg].size()<<std::endl;
      nmatch += crthitsgroup[itrg].size() * crthitsgroup[jtrg].size();
    }
  }

  //std::cout<<"Total matches: "<<nmatch<<" "<<hits3d_F.size()<<" "<<hits3d_B.size()<<std::endl;

  // Reconstruciton information
  art::Handle < std::vector < recob::Track > > trackListHandle;
  std::vector < art::Ptr < recob::Track > > trackList;
  if (e.getByLabel(fTrackModuleLabel, trackListHandle)) {
    art::fill_ptr_vector(trackList, trackListHandle);
  }
  art::FindManyP < recob::Hit > hitsFromTrack(trackListHandle, e, fTrackModuleLabel);
  for (auto const& track : trackList){
    
    std::vector< art::Ptr<recob::Hit> > allHits =  hitsFromTrack.at(track.key());
    if (allHits.empty()) continue;

    double trackStartPositionX = track->Vertex().X();
    double trackStartPositionY = track->Vertex().Y();
    double trackStartPositionZ = track->Vertex().Z();
    double trackEndPositionX = track->End().X();
    double trackEndPositionY = track->End().Y();
    double trackEndPositionZ = track->End().Z();
    TVector3 trackStart(trackStartPositionX, trackStartPositionY, trackStartPositionZ);
    TVector3 trackEnd(trackEndPositionX, trackEndPositionY, trackEndPositionZ);

    if ((trackEndPositionZ > 660 && trackStartPositionZ < 50) || (trackStartPositionZ > 660 && trackEndPositionZ < 50)) {

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
      for (auto const& fronthit : hits3d_F){
        for (auto const& backhit : hits3d_B){
          if (int(util::absDiff(fronthit.t, backhit.t))>fFronttoBackTimingCut) continue;
          //t0 correction
          double xOffset = 0;
          int RDOffset=0;
          if (!fMCCSwitch) RDOffset=111;
          double ticksOffset = 0;
          if (!fMCCSwitch) ticksOffset = ((fronthit.t+backhit.t)/2.-rdtimestamp+RDOffset)/25.f+detectorPropertiesService->GetXTicksOffset(allHits[0]->WireID());
          else ticksOffset = (fronthit.t+backhit.t)/2./500.f+detectorPropertiesService->GetXTicksOffset(allHits[0]->WireID());
		
          xOffset=detectorPropertiesService->ConvertTicksToX(ticksOffset,allHits[0]->WireID());

          trackStartPositionX=trackStartPositionX-xOffset;
          trackEndPositionX=trackEndPositionX-xOffset;

          TVector3 v1(fronthit.x, fronthit.y, fronthit.z);
          TVector3 v2(backhit.x, backhit.y, backhit.z);

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
          
          double deltaX1 = (predictedHitPositionX1-v1.X());
          double deltaX2 = (predictedHitPositionX2-v2.X());
          double deltaX  = std::abs(deltaX2)+std::abs(deltaX1);

          double deltaY1 = (predictedHitPositionY1-v1.Y());
          double deltaY2 = (predictedHitPositionY2-v2.Y());
          double deltaY  = std::abs(deltaY2)+std::abs(deltaY1); 
	  //if (std::abs(deltaY2)<20) std::cout<<v2.Y()<<std::endl;      
	  if (min_delta > std::abs(deltaX) + std::abs(deltaY)){
            min_delta = std::abs(deltaX) + std::abs(deltaY);
            best_XF = fronthit.x;
            best_YF = fronthit.y;
            best_ZF = fronthit.z;
            best_XB = backhit.x;
            best_YB = backhit.y;
            best_ZB = backhit.z;
            best_dotProductCos = dotProductCos;
            best_deltaXF = deltaX1;
            best_deltaYF = deltaY1;
            best_deltaXB = deltaX2;
            best_deltaYB = deltaY2;
	    best_trigXF=fronthit.trigX;
	    best_trigYF=fronthit.trigY;
	    best_trigXB=backhit.trigX;
	    best_trigYB=backhit.trigY;
	   if (!fMCCSwitch) best_T=((fronthit.t+backhit.t)/2-rdtimestamp+RDOffset)/50.f;
	   else best_T=(fronthit.t+backhit.t)/1000.f;
          }
        }
      }
      if (std::abs(best_dotProductCos)>0.99 && std::abs(best_deltaXF)+std::abs(best_deltaXB)<40 && std::abs(best_deltaYF)+std::abs(best_deltaYB)<40 ) {
        std::cout<<"Track ID = "<<track.key()<<" best_deltaXF =  "<<best_deltaXF<<" best_deltaYF = "<<best_deltaYF<<" best_deltaXB = "<<best_deltaXB<<" best_deltaYB = "<<best_deltaYB<<" best_XF = "<<best_XF<<" best_YF = "<<best_YF<<" best_XB = "<<best_XB<<" best_YB = "<<best_YB<<" best_dotProductCos = "<<best_dotProductCos<<std::endl;
        std::vector<float> hitF;
	std::vector<float> hitB;
	hitF.push_back(best_XF); hitF.push_back(best_YF); hitF.push_back(best_ZF);
	hitB.push_back(best_XB); hitB.push_back(best_YB); hitB.push_back(best_ZB);

        CRTTrack->push_back(anab::CosmicTag(hitF,hitB, std::abs(best_dotProductCos),anab::CosmicTagID_t::kNotTagged));

        CRTT0->push_back(anab::T0(best_T, 13,2,track.key(),best_dotProductCos));
        util::CreateAssn(*this, e, *CRTTrack, track, *TPCCRTassn);
 	util::CreateAssn(*this, e, *CRTT0, track, *TPCT0assn);
        util::CreateAssn(*this, e, *CRTTrack, crtList[best_trigXF], *CRTTrigassn);
        util::CreateAssn(*this, e, *CRTTrack, crtList[best_trigYF], *CRTTrigassn);

        util::CreateAssn(*this, e, *CRTTrack, crtList[best_trigXB], *CRTTrigassn);

        util::CreateAssn(*this, e, *CRTTrack, crtList[best_trigYB], *CRTTrigassn);


      }
    }
  }
  e.put(std::move(CRTTrack)); 
  e.put(std::move(TPCCRTassn));
  e.put(std::move(CRTT0)); 
  e.put(std::move(TPCT0assn));
  e.put(std::move(CRTTrigassn));
}

// v6 Geo Channel Map
bool CRT::TwoCRTMatchingProducer::moduleMatcher(int module1, int module2) {
  // Function checking if two hits could reasonably be matched into a 2D hit
  if ((module1 == 6 && (module2 == 10 || module2 == 11)) || (module1 == 14 && (module2 == 10 || module2 == 11)) || (module1 == 19 && (module2 == 26 || module2 == 27)) || (module1 == 31 && (module2 == 26 || module2 == 27)) || (module1 == 7 && (module2 == 12 || module2 == 13)) || (module1 == 15 && (module2 == 12 || module2 == 13)) || (module1 == 18 && (module2 == 24 || module2 == 25)) || (module1 == 30 && (module2 == 24 || module2 == 25)) || (module1 == 1 && (module2 == 4 || module2 == 5)) || (module1 == 9 && (module2 == 4 || module2 == 5)) || (module1 == 16 && (module2 == 20 || module2 == 21)) || (module1 == 28 && (module2 == 20 || module2 == 21)) || (module1 == 0 && (module2 == 2 || module2 == 3)) || (module1 == 8 && (module2 == 2 || module2 == 3)) || (module1 == 17 && (module2 == 22 || module2 == 23)) || (module1 == 29 && (module2 == 22 || module2 == 23))) return 1;
  else return 0;

}


DEFINE_ART_MODULE(CRT::TwoCRTMatchingProducer)
