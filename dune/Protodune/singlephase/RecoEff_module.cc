////////////////////////////////////////////////////////////////////////
// Class:       RecoEff
// Module Type: analyzer
// File:        RecoEff_module.cc
// Author:      Dorota Stefan
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTracker.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "TH1.h"
#include "TEfficiency.h"
#include "TTree.h"

namespace proto
{
	class RecoEff;
}

class proto::RecoEff : public art::EDAnalyzer {
public:

  struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> SimulationLabel {
          Name("SimulationLabel"), Comment("tag of simulation producer")
      };

      fhicl::Atom<art::InputTag> HitModuleLabel {
          Name("HitModuleLabel"), Comment("tag of hits producer")
      };

      fhicl::Atom<art::InputTag> TrackModuleLabel {
          Name("TrackModuleLabel"), Comment("tag of track producer")
      };
    
      fhicl::Atom<size_t> MinHitsPerPlane {
          Name("MinHitsPerPlane"), Comment("min hits per plane for a reconstructable MC particle")
      };
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit RecoEff(Parameters const& config);

  RecoEff(RecoEff const &) = delete;
  RecoEff(RecoEff &&) = delete;
  RecoEff & operator = (RecoEff const &) = delete;
  RecoEff & operator = (RecoEff &&) = delete;

  void analyze(art::Event const & e) override;
  void beginJob() override;
  void endJob() override;

private:
  void ResetVars();
  
  TH1D* fDenominatorHist;
  TH1D* fNominatorHist;
  
  TTree *fTree;
  
  int fRun;
  int fEvent;
  int fNRecoTracks;
  int fReconstructable;
  int fMatched;

  double fEnGen;
  double fEkGen;
  double fT0;
  
  art::InputTag fSimulationLabel;
  art::InputTag fHitModuleLabel;
  art::InputTag fTrackModuleLabel;
  size_t fMinHitsPerPlane;
};

proto::RecoEff::RecoEff(Parameters const& config) : EDAnalyzer(config),
    fSimulationLabel(config().SimulationLabel()),
    fHitModuleLabel(config().HitModuleLabel()),
    fTrackModuleLabel(config().TrackModuleLabel()),
    fMinHitsPerPlane(config().MinHitsPerPlane())
{
}

void proto::RecoEff::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	
	fDenominatorHist = tfs->make<TH1D>("Denominator", "all reconstructable particles", 100., 0., 1500.);
	fNominatorHist = tfs->make<TH1D>("Nominator", "reconstructed and matched tracks", 100., 0., 1500.);

	fTree = tfs->make<TTree>("events", "summary tree");
	fTree->Branch("fRun", &fRun, "fRun/I");
	fTree->Branch("fEvent", &fEvent, "fEvent/I");
	fTree->Branch("fEnGen", &fEnGen, "fEnGen/D");
	fTree->Branch("fEkGen", &fEkGen, "fEkGen/D");
	fTree->Branch("fT0", &fT0, "fT0/D");
	fTree->Branch("fNRecoTracks", &fNRecoTracks, "fNRecoTracks/I");
	fTree->Branch("fReconstructable", &fReconstructable, "fReconstructable/I");
	fTree->Branch("fMatched", &fMatched, "fMatched/I");
}

void proto::RecoEff::endJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    TEfficiency* fEfficiency = tfs->makeAndRegister<TEfficiency>("Efficiency", "tracking efficiency", *fNominatorHist, *fDenominatorHist);
    std::cout << "Efficiency created: " << fEfficiency->GetTitle() << std::endl;
}

void proto::RecoEff::analyze(art::Event const & e)
{
  ResetVars();
  
  fRun = e.run();
  fEvent = e.id().event();

  art::ServiceHandle<cheat::BackTracker> bt;
  
  // we are going to look only for these MC truth particles, which contributed to hits
  // and normalize efficiency to things which indeed generated activity and hits in TPC's:
  // - need a vector of hits associated to these MC particles,
  // - need the energy deposit corresponding to the particle part seen in these hits.
  std::unordered_map<int, std::vector<recob::Hit> > mapTrackIDtoHits;
  std::unordered_map<int, double> mapTrackIDtoHitsEnergy;

  auto const & hitListHandle = *e.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
  for (auto const & h : hitListHandle)
  {
    std::unordered_map<int, double> particleID_E;
    for (auto const & id : bt->HitToTrackID(h)) // loop over std::vector< sim::TrackIDE > contributing to hit h
	{
	    // select only hadronic and muon track, skip EM activity (electron's pdg, negative track id)
        if ((id.trackID > 0) && (abs((bt->TrackIDToParticle(id.trackID))->PdgCode()) != 11))
        {
            particleID_E[id.trackID] += id.energy;
        }
    }

    int best_id = 0;
    double max_e = 0.0;
    for (auto const & contrib : particleID_E)
    {
        if (contrib.second > max_e) // find particle TrackID corresponding to max energy contribution
        {
            max_e = contrib.second;
            best_id = contrib.first;
        }
    }

    if (max_e > 0)
    {
        mapTrackIDtoHits[best_id].push_back(h);
        mapTrackIDtoHitsEnergy[best_id] += max_e;
    }
  }
  // ---------------------------------------------------------------------------------------


  // now lets map particles to their associated hits, but with some filtering of things
  // which have a minimal chance to be reconstructed: so at least some hits in two views
  std::unordered_map<int, std::vector< recob::Hit > > mapTrackIDtoHits_filtered;
  for (auto const & p : mapTrackIDtoHits)
  {
    // *** here more conditions may be applied to select interesting MC particles, eg:
    // - check with BackTracker if this is a beam/cosmic particle:
    //auto origin = bt->TrackIDToMCTruth(p.first)->Origin();
    //if (origin == simb::kCosmicRay) { std::cout << "cosmic" << std::endl; }
    //else if (origin == simb::kSingleParticle) { std::cout << "beam" << std::endl; }
    //else { std::cout << "other" << std::endl; }
    // - check if this the process name is "primary" to select primary particles:
    //std::cout << bt->TrackIDToParticle(p.first)->Process() << std::endl;
    // - check if the mother of this MC particle is primary so you can select secondaries of the primary interaction

    std::unordered_map<geo::View_t, size_t> hit_count;
    for (auto const & h : p.second) { hit_count[h.View()]++; }

    size_t nviews = 0;
    for (auto const & n : hit_count) {  if (n.second > fMinHitsPerPlane) { nviews++; } }
    if (nviews >= 2) { mapTrackIDtoHits_filtered.emplace(p); }
  }
  fReconstructable = mapTrackIDtoHits_filtered.size();

  std::cout << "------------ event: " << fEvent << " has reconstructable: " << fReconstructable << std::endl;
  for (auto const &p : mapTrackIDtoHits_filtered)
  {
    std::cout << " : id " << p.first << " size " << p.second.size() << " en: " << mapTrackIDtoHitsEnergy[p.first] << std::endl;
    fDenominatorHist->Fill(p.second.size());
  }
  // ---------------------------------------------------------------------------------------


  // match reconstructed tracks to MC particles
  const auto trkHandle = e.getValidHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  art::FindManyP< recob::Hit > hitsFromTracks(trkHandle, e, fTrackModuleLabel);

  fNRecoTracks = trkHandle->size();
  for (size_t t = 0; t < trkHandle->size(); ++t)     // loop over tracks
  {
    // *** here you should select if the reconstructed track is interesting, eg:
    // - associated PFParticle has particlular PDG code
    // - or just the recob::Track has interesting ParticleID

    std::unordered_map<int, double> trkID_E;         // map MC particles to their energy contributed to this track t
    const auto & hits = hitsFromTracks.at(t);
    for (const auto & h : hits)                      // loop over hits assigned to track t
    {
        for (auto const & ide : bt->HitToTrackID(h)) // loop over std::vector< sim::TrackIDE >, for a hit h
        {
            if (mapTrackIDtoHits_filtered.find(ide.trackID) != mapTrackIDtoHits_filtered.end())
            { 
                trkID_E[ide.trackID] += ide.energy;  // ...and sum energies contributed to this reco track
            }                                        // by MC particles which we considered as reconstructable
        }
    }

    // find MC particle which cotributed maximum energy to track t
    int best_id = 0;
    double max_e = 0, tot_e = 0;
    for (auto const & entry : trkID_E)
	{
        tot_e += entry.second;    // sum total energy in the hits in track t
        if (entry.second > max_e) // find track ID corresponding to max energy
        {
            max_e = entry.second;
            best_id = entry.first;
        }
    }

    // check if reco track is matching to MC particle:
    if ((max_e > 0.5 * tot_e) &&                         // MC particle has more than 50% energy contribution to the track
        (max_e > 0.5 * mapTrackIDtoHitsEnergy[best_id])) // track covers more than 50% of energy deposited by MC particle in hits
    {
        fNominatorHist->Fill(mapTrackIDtoHits_filtered[best_id].size());
        fMatched++;
    }
  }

  std::cout << "------------ reconstructed: " << fNRecoTracks << " and then matched to MC particles: " << fMatched << std::endl;
  fTree->Fill();
}

void proto::RecoEff::ResetVars()
{
	fRun = 0;
	fEvent = 0;
	fEnGen = 0.0;
	fEkGen = 0.0;
	fT0 = 0.0;
	fMatched = 0;
	fNRecoTracks = 0;
	fReconstructable = 0;
}

DEFINE_ART_MODULE(proto::RecoEff)
