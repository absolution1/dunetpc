#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

protoana::ProtoDUNETrackUtils::ProtoDUNETrackUtils(){

}

protoana::ProtoDUNETrackUtils::~ProtoDUNETrackUtils(){

}

std::vector<anab::CosmicTag> protoana::ProtoDUNETrackUtils::GetRecoTrackCosmicTag(const recob::Track &track, art::Event const &evt, std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  unsigned int trackIndex = track.ID();

  // Convert to std::vector<anab::CosmicTag> from std::vector<art::Ptr<anab::CosmicTag>>
  std::vector<anab::CosmicTag> trackTags;

  try{
    const art::FindManyP<anab::CosmicTag> findCosmicTags(recoTracks,evt,trackModule);
    for(unsigned int t = 0; t < findCosmicTags.at(trackIndex).size(); ++t){
      trackTags.push_back((*(findCosmicTags.at(trackIndex)[t])));
    }
  }
  catch(...){
//    std::cerr << "Product not found - returning empty vector" << std::endl;
  }

  return trackTags;    
}

std::vector<anab::T0> protoana::ProtoDUNETrackUtils::GetRecoTrackT0(const recob::Track &track, art::Event const &evt, std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  unsigned int trackIndex = track.ID();

  // Convert to std::vector<anab::T0> from std::vector<art::Ptr<anab::T0>>
  std::vector<anab::T0> trackT0s;
  
  try{
    const art::FindManyP<anab::T0> findTrackT0s(recoTracks,evt,trackModule);
    for(unsigned int t = 0; t < findTrackT0s.at(trackIndex).size(); ++t){
      trackT0s.push_back((*(findTrackT0s.at(trackIndex)[t])));
    }
  }
  catch(...){
//    std::cerr << "Product not found - returning empty vector" << std::endl;
  }
  
  return trackT0s;

}

// Get the Calorimetry(s) from a given reco track
std::vector<anab::Calorimetry> protoana::ProtoDUNETrackUtils::GetRecoTrackCalorimetry(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  std::vector<anab::Calorimetry> caloInfo;
  
  try{
    const art::FindManyP<anab::Calorimetry> findCalorimetry(recoTracks,evt,caloModule);
    std::vector<art::Ptr<anab::Calorimetry>> theseCalos = findCalorimetry.at(track.ID());

    for( auto calo : theseCalos){
      caloInfo.push_back(*calo);
    }
  }
  catch(...){
    std::cerr << "No calorimetry object found... returning empty vector" << std::endl;
  }

  return caloInfo;
}

// Get the hits from a given reco track
const std::vector<const recob::Hit*> protoana::ProtoDUNETrackUtils::GetRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,trackModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(track.ID());

  std::vector<const recob::Hit*> trackHits;

  for(const art::Ptr<recob::Hit> hit : inputHits){

    trackHits.push_back(hit.get());

  }

  return trackHits;  

}

// Get the hits from a given reco track
unsigned int protoana::ProtoDUNETrackUtils::GetNumberRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const{

  return GetRecoTrackHits(track,evt,trackModule).size();

}


