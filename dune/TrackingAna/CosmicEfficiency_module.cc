////////////////////////////////////////////////////////////////////////
// Class:       CosmicEfficiency
// Plugin Type: analyzer (art v2_07_03)
// File:        CosmicEfficiency_module.cc
//
// Generated at Mon Sep  4 06:55:33 2017 by Leigh Whitehead using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "TH1F.h"

namespace cosmic {
  class CosmicEfficiency;
}


class cosmic::CosmicEfficiency : public art::EDAnalyzer {
public:

  // Just make an enum for plotting in a histogram to avoid
  // histograms of the pdg code. These definitions are
  // valid for particle and antiparticle
  // Note the enum starts at one since this is the bin
  // number for the histograms
  enum particleType
  {
    kNotSet = 0,
    kElectron = 1,
    kMuon = 2,
    kPion = 3,
    kProton = 4,
    kKaon = 5,
    kOther = 6
  };

  explicit CosmicEfficiency(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicEfficiency(CosmicEfficiency const &) = delete;
  CosmicEfficiency(CosmicEfficiency &&) = delete;
  CosmicEfficiency & operator = (CosmicEfficiency const &) = delete;
  CosmicEfficiency & operator = (CosmicEfficiency &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  const simb::MCParticle* GetTruthParticle(const std::vector< art::Ptr<recob::Hit> > & hits) const;

  const particleType ConvertPDGCode(int pdg) const;

  // fcl parameters
  art::InputTag fTrackerTag;

  // Useful counters
  unsigned int fTrueCosmicMuon;
  unsigned int fTaggedCosmicMuon;
  unsigned int fTrueOther;
  unsigned int fTaggedOther;
  // Also tag just anything from cosmic or beam origin
  unsigned int fTrueCosmicOrigin;
  unsigned int fTaggedCosmicOrigin;
  unsigned int fTrueOtherOrigin;
  unsigned int fTaggedOtherOrigin;

  // Histograms
  TH1F *hCosmicMuonEff;
  TH1F *hCosmicMuonPur;
  TH1F *hCosmicEff;
  TH1F *hCosmicPur;
  TH1F *hTaggedPDG;
};


cosmic::CosmicEfficiency::CosmicEfficiency(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fTrackerTag(p.get<art::InputTag>("TrackerTag")),
  fTrueCosmicMuon(0), fTaggedCosmicMuon(0),
  fTrueOther(0), fTaggedOther(0),
  fTrueCosmicOrigin(0), fTaggedCosmicOrigin(0),
  fTrueOtherOrigin(0), fTaggedOtherOrigin(0)
{

}

void cosmic::CosmicEfficiency::beginJob()
{

  // Get a file service so we can make some histograms to write out
  art::ServiceHandle<art::TFileService> tfs;

  hCosmicMuonEff = tfs->make<TH1F>("hCosmicMuonEff",";Efficiency",50,0,1.0001);
  hCosmicMuonPur = tfs->make<TH1F>("hCosmicMuonPur",";Purity",50,0,1.0001);
  hCosmicEff = tfs->make<TH1F>("hCosmicEff",";Efficiency",50,0,1.0001);
  hCosmicPur = tfs->make<TH1F>("hCosmicPur",";Purity",50,0,1.0001);
  hTaggedPDG = tfs->make<TH1F>("hTaggedParticleType","",6,0.5,6.5);

  // Let's give the particle type plot useful labels
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kElectron,"Electron");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kMuon,"Muon");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kPion,"Pion");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kProton,"Proton");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kKaon,"Kaon");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kOther,"Other");
}

void cosmic::CosmicEfficiency::analyze(art::Event const & evt)
{

  // We must have MC for this module to make sense
  if(evt.isRealData()) return;

  // Get the reconstructed tracks
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);

  // We need the association between the tracks and the hits
  const art::FindManyP<recob::Hit> findTrackHits(recoTracks, evt, fTrackerTag);

  // Also want the association between tracks and anab::CosmicTag objects
  const art::FindManyP<anab::CosmicTag> findCosmicTags(recoTracks,evt,fTrackerTag);

  // Finally, get the backtracker
  art::ServiceHandle<cheat::BackTracker> bt;

  // Since we will make plots of efficiency per event, we need local versions of the class variables
  unsigned int local_TrueCosmicMuon = 0; 
  unsigned int local_TaggedCosmicMuon = 0;
  unsigned int local_TrueOther = 0;
  unsigned int local_TaggedOther = 0;
  unsigned int local_TrueCosmicOrigin = 0;
  unsigned int local_TaggedCosmicOrigin = 0;
  unsigned int local_TrueOtherOrigin = 0;
  unsigned int local_TaggedOtherOrigin = 0;

  // Loop over the tracks
  for(unsigned int t = 0; t < recoTracks->size(); ++t){
  
    // Get the hits
    auto const& trackHits = findTrackHits.at(t);

    // Get the best match MC Particle for this track
    const simb::MCParticle* trueMatch = GetTruthParticle(trackHits);

    if(trueMatch == 0x0) continue;

    // We need to use the backtracker to find if this true particle is cosmic or not
    bool isCosmic = false;
    bool isCosmicMuon = false;
    const art::Ptr<simb::MCTruth> mc = bt->TrackIDToMCTruth(trueMatch->TrackId());

    // Check if we have a cosmic muon (or anti-muon) and set the true counting variables
    if(mc->Origin() == simb::kCosmicRay){
      if(fabs(trueMatch->PdgCode()) == 13){
        ++local_TrueCosmicMuon;
        isCosmicMuon = true;
      }
      else{
        ++local_TrueOther;
      }
      ++local_TrueCosmicOrigin;
      isCosmic = true;
    }
    else{
      // Definitely not cosmic!
      ++local_TrueOther;
      ++local_TrueOtherOrigin;
    }

    // Now look for a cosmic tag object
    auto const& trackCosmicTag = findCosmicTags.at(t);
    
    // Did we find a cosmic tag for this track?
    if(trackCosmicTag.size() == 0){
      continue;
    }
    else{
      // Specific cosmic muon counters
      if(isCosmicMuon){
        ++local_TaggedCosmicMuon;
      }
      else{
        ++local_TaggedOther;
      }
      // Origin counters
      if(isCosmic){
        ++local_TaggedCosmicOrigin;
      }
      else{
        ++local_TaggedOtherOrigin;
      }
      // Fill the tagged particle hisrogram
      hTaggedPDG->Fill(ConvertPDGCode(trueMatch->PdgCode()));
    }

  } // End loop over reconstructed tracks

  // Fill the plots
  if(local_TrueCosmicMuon > 0){
    hCosmicMuonEff->Fill(local_TaggedCosmicMuon/static_cast<float>(local_TrueCosmicMuon));
  }
  if(local_TaggedCosmicMuon+local_TaggedOther > 0){
    hCosmicMuonPur->Fill(local_TaggedCosmicMuon/static_cast<float>(local_TaggedCosmicMuon + local_TaggedOther));
  }
  if(local_TrueCosmicOrigin > 0){
    hCosmicEff->Fill(local_TaggedCosmicOrigin/static_cast<float>(local_TrueCosmicOrigin));
  }
  if(local_TaggedCosmicOrigin+local_TaggedOtherOrigin > 0){
    hCosmicPur->Fill(local_TaggedCosmicOrigin/static_cast<float>(local_TaggedCosmicOrigin + local_TaggedOtherOrigin));
  } 
  
  // Make sure to increment the main member counters
  fTrueCosmicMuon += local_TrueCosmicMuon;
  fTaggedCosmicMuon += local_TaggedCosmicMuon;
  fTrueOther += local_TrueOther;
  fTaggedOther += local_TaggedOther;

  fTrueCosmicOrigin += local_TrueCosmicOrigin;
  fTaggedCosmicOrigin += local_TaggedCosmicOrigin;
  fTrueOtherOrigin += local_TrueOtherOrigin;
  fTaggedOtherOrigin += local_TaggedOtherOrigin;
}

void cosmic::CosmicEfficiency::endJob()
{
  if(fTrueCosmicMuon != 0 && (fTaggedCosmicMuon+fTaggedOther) !=0){
    std::cout << " == Cosmic Efficiency Summary == " << std::endl;
    std::cout << " - Looking for cosmic muons: " << std::endl;
    std::cout << " Tagged " << fTaggedCosmicMuon << " of " << fTrueCosmicMuon << " :: " << 100.*fTaggedCosmicMuon/static_cast<float>(fTrueCosmicMuon) << "% efficiency" << std::endl;
    std::cout << " Total tags = " << fTaggedCosmicMuon + fTaggedOther << " :: " 
              << 100.*fTaggedCosmicMuon/static_cast<float>(fTaggedCosmicMuon + fTaggedOther) << "% purity" << std::endl;
    std::cout << " - Looking for just cosmic origin: " << std::endl;
    std::cout << " Tagged " << fTaggedCosmicOrigin << " of " << fTrueCosmicOrigin << " :: " << 100.*fTaggedCosmicOrigin/static_cast<float>(fTrueCosmicOrigin) << "% efficiency" << std::endl;
    std::cout << " Total tags = " << fTaggedCosmicOrigin + fTaggedOtherOrigin << " :: " 
              << 100.*fTaggedCosmicOrigin/static_cast<float>(fTaggedCosmicOrigin + fTaggedOtherOrigin) << "% purity" << std::endl;
  }
}

// Function to find the best matched true particle to a reconstructed particle
const simb::MCParticle* cosmic::CosmicEfficiency::GetTruthParticle(const std::vector< art::Ptr<recob::Hit> > & hits) const
{
  const simb::MCParticle* mcParticle = 0;

  art::ServiceHandle<cheat::BackTracker> bt;
  std::unordered_map<int, double> trkIDE;
  for (auto const & h : hits)
  {
    for (auto const & ide : bt->HitToTrackID(h)) // loop over std::vector<sim::TrackIDE>
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
    mcParticle = bt->TrackIDToParticle(best_id); // MCParticle corresponding to track ID
  }

  return mcParticle;
}

// Function to convert the PDG code into the scheme we want to use for our histogram
const cosmic::CosmicEfficiency::particleType cosmic::CosmicEfficiency::ConvertPDGCode(int pdg) const
{
  cosmic::CosmicEfficiency::particleType particle = cosmic::CosmicEfficiency::kNotSet;

  switch(pdg)
  {
    case -11:
    case  11: 
      particle = cosmic::CosmicEfficiency::kElectron;
      break;
    case -13:
    case  13: 
      particle = cosmic::CosmicEfficiency::kMuon;
      break;
    case  211:
    case -211:
      particle = cosmic::CosmicEfficiency::kPion;
      break;
    case  2212:
    case -2212:
      particle = cosmic::CosmicEfficiency::kProton;
      break;
    case  321:
    case -321:
      particle = cosmic::CosmicEfficiency::kKaon;
      break;
    default:
      particle = cosmic::CosmicEfficiency::kOther;
  }

  return particle;
}

DEFINE_ART_MODULE(cosmic::CosmicEfficiency)

