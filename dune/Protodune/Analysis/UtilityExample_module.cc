////////////////////////////////////////////////////////////////////////
// Class:       UtilityExample
// Plugin Type: analyzer (art v2_07_03)
// File:        UtilityExample_module.cc
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

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"

namespace protoana {
  class UtilityExample;
}


class protoana::UtilityExample : public art::EDAnalyzer {
public:

  explicit UtilityExample(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UtilityExample(UtilityExample const &) = delete;
  UtilityExample(UtilityExample &&) = delete;
  UtilityExample & operator = (UtilityExample const &) = delete;
  UtilityExample & operator = (UtilityExample &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // fcl parameters
  std::string fTrackerTag;
  std::string fGeneratorTag;
  bool fVerbose;

};


protoana::UtilityExample::UtilityExample(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose"))
{

}

void protoana::UtilityExample::beginJob()
{

}

void protoana::UtilityExample::analyze(art::Event const & evt)
{

  // We must have MC for this module to make sense
  if(evt.isRealData()) return;

  // Get the reconstructed tracks
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);

  // Bag ourselves a couple of utilities
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNETrackUtils trackUtil;

  // Loop over the tracks
  for(unsigned int t = 0; t < recoTracks->size(); ++t){

    // Match to truth
    const recob::Track thisTrack = (*recoTracks)[t];
    const simb::MCParticle *trueMatch = truthUtil.GetMCParticleFromRecoTrack(thisTrack,evt,fTrackerTag);
    bool hasTruth = (trueMatch != 0x0);
 
    // Check for T0
    std::vector<anab::T0> trackT0 = trackUtil.GetRecoTrackT0(thisTrack,evt,fTrackerTag);
    bool hasT0 = (trackT0.size() != 0);

    // Check for cosmic tag
    std::vector<anab::CosmicTag> trackCosmic = trackUtil.GetRecoTrackCosmicTag(thisTrack,evt,fTrackerTag);
    bool hasTag = (trackCosmic.size() != 0);

    std::cout << "Track properties: ";
    if(hasTruth) std::cout << " Matched true particle of type " << trueMatch->PdgCode() << " :: ";
    if(hasT0)    std::cout << " Track has T0 = " << trackT0[0].fTime << " :: ";
    if(hasTag)   std::cout << " Has cosmic tag ";
    std::cout << std::endl;

  } // End loop over reconstructed tracks

  // Get the generator MCTruth objects and find the GEANT track id of the good particle
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
  if(geantGoodParticle != 0x0){
    std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() << std::endl;
  }

}

void protoana::UtilityExample::endJob()
{

}

DEFINE_ART_MODULE(protoana::UtilityExample)

