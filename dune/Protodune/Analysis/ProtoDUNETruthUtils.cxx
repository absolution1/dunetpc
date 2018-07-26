#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Event.h"

protoana::ProtoDUNETruthUtils::ProtoDUNETruthUtils(){

}

protoana::ProtoDUNETruthUtils::~ProtoDUNETruthUtils(){

}

// Function to find the best matched true particle to a reconstructed particle. In case of problems, returns a null pointer
const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetMCParticleFromRecoTrack(const recob::Track &track, art::Event const & evt, std::string trackModule) const{

  const simb::MCParticle* mcParticle = 0x0;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return mcParticle;

  // Get the reconstructed tracks
  auto allRecoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  // We need the association between the tracks and the hits
  const art::FindManyP<recob::Hit> findTrackHits(allRecoTracks, evt, trackModule);

  // A track utils object will be useful
  protoana::ProtoDUNETrackUtils utils;
  unsigned int trackIndex = utils.GetTrackIndexNumber(track,evt,trackModule);

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<int, double> trkIDE;
  for (auto const & h : findTrackHits.at(trackIndex))
  {
    for (auto const & ide : bt_serv->HitToTrackIDEs(h)) // loop over std::vector<sim::TrackIDE>
    {
        trkIDE[ide.trackID] += ide.energy; // sum energy contribution by each track ID
    }
  }

  int best_id = 0;
  double tot_e = 0, max_e = 0;
  for (auto const & contrib : trkIDE)
  {
    tot_e += contrib.second;     // sum total energy in these hits
    if (contrib.second > max_e)  // find track ID corresponding to max energy
    {
        max_e = contrib.second;
        best_id = contrib.first;
    }
  }

  if ((max_e > 0) && (tot_e > 0)) // ok, found something reasonable
  {
    if (best_id < 0)            // NOTE: negative ID means this is EM activity
    {                           // caused by track with the same but positive ID
//        best_id = -best_id;     // --> we'll find mother MCParticle of these hits
      return mcParticle;
    }
    mcParticle = pi_serv->TrackIdToParticle_P(best_id); // MCParticle corresponding to track ID
  }

  return mcParticle;
} 

