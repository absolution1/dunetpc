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
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"

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
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  bool fVerbose;

};


protoana::UtilityExample::UtilityExample(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
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

  // Truth utility
  protoana::ProtoDUNETruthUtils truthUtil;

  // Get the generator MCTruth objects and find the GEANT track id of the good particle
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
  if(geantGoodParticle != 0x0){
    std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() << std::endl;
  }

  // Track utility
  protoana::ProtoDUNETrackUtils trackUtil;

  unsigned int nTracksWithTruth = 0;
  unsigned int nTracksWithT0    = 0;
  unsigned int nTracksWithTag   = 0;

  // Get the reconstructed tracks
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);

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

    if(hasTruth) ++nTracksWithTruth;
    if(hasT0)    ++nTracksWithT0;
    if(hasTag)   ++nTracksWithTag;

  } // End loop over reconstructed tracks

  std::cout << "Found " << recoTracks->size() << " reconstructed tracks:" << std::endl;
  std::cout << " - " << nTracksWithTruth << " successfully associated to the truth information " << std::endl;
  std::cout << " - " << nTracksWithT0    << " have a reconstructed T0" << std::endl;
  std::cout << " - " << nTracksWithTag   << " have a cosmic tag" << std::endl;

  // What about PFParticles?
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  unsigned int nParticlesPrimary   = 0;
  unsigned int nParticlesWithT0    = 0;
  unsigned int nParticlesWithTag   = 0;
  for(unsigned int p = 0; p < recoParticles->size(); ++p){
  
    // Only consider primary particles here
    if(!(recoParticles->at(p).IsPrimary())) continue;
 
    ++nParticlesPrimary;

    // Do we have a T0?
    if(pfpUtil.GetPFParticleT0(recoParticles->at(p),evt,fPFParticleTag).size() != 0) ++nParticlesWithT0;
    if(pfpUtil.GetPFParticleCosmicTag(recoParticles->at(p),evt,fPFParticleTag).size() != 0) ++nParticlesWithTag;

  } 
  std::cout << "Found " << nParticlesPrimary << " reconstructed primary particles:" << std::endl;
  std::cout << " - Total number of particles = " << recoParticles->size() << std::endl;
  std::cout << " - " << nParticlesWithT0    << " have a reconstructed T0" << std::endl;
  std::cout << " - " << nParticlesWithTag   << " have a cosmic tag" << std::endl;

  // We can also build a slice map if we want to
  std::map<unsigned int, std::vector<recob::PFParticle*>> sliceMap;
  sliceMap = pfpUtil.GetPFParticleSliceMap(evt,fPFParticleTag);
  unsigned int nMultiSlice = 0;
  std::cout << "Found " << sliceMap.size() << " slices with PFParticles" << std::endl;
  for(auto slice : sliceMap){
    if(slice.second.size() > 1) ++nMultiSlice;
  }
  std::cout << " - " << nMultiSlice << " have at least one primary PFParticle" << std::endl;

  // Look for the beam particle slice and get the particles if we can
  unsigned short beamSlice = pfpUtil.GetBeamSlice(evt,fPFParticleTag);
  std::cout << "- Beam slice = " << beamSlice << std::endl;
  if(beamSlice != 9999){
    std::vector<recob::PFParticle*> beamSlicePrimaries = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
    std::cout << " - Found the beam slice! " << beamSlicePrimaries.size() << " beam particles in slice " << beamSlice << std::endl;

    for (const recob::PFParticle* prim : beamSlicePrimaries){
      // Vertex of the PFParticle
      const TVector3 vtx = pfpUtil.GetPFParticleVertex(*prim,evt,fPFParticleTag,fTrackerTag);
      std::cout << "Beam particle vertex: " << std::endl;
      vtx.Print();
      // For track-like PFParticles, returns the end of the associated track. This should give
      // the secondary interaction point
      if(pfpUtil.IsPFParticleTracklike(*prim)){
        const TVector3 sec = pfpUtil.GetPFParticleSecondaryVertex(*prim,evt,fTrackerTag,fPFParticleTag);
        std::cout << "Beam particle interaction vertex: " << std::endl;
        sec.Print();
      }
      // If we want to look at the track or shower that makes up the PFParticle:
      /*
      if(pfpUtil.IsPFParticleTracklike(*prim)){
        const recob::Track* trk = pfpUtil.GetPFParticleTrack(*prim,evt,fPFParticleTag,fTrackerTag);
        if(trk != 0x0) std::cout << "Got track" << std::endl;
      }
      else{
        const recob::Shower* shw = pfpUtil.GetPFParticleShower(*prim,evt,fPFParticleTag,fShowerTag);
        if(shw != 0x0) std::cout << "Got shower" << std::endl;
      }
      */

    }
  }

  // See how many clear cosmics we have
  std::cout << "There are " << pfpUtil.GetClearCosmicPFParticles(evt,fPFParticleTag).size() << " clear cosmic muons in this event" << std::endl;

}

void protoana::UtilityExample::endJob()
{

}

DEFINE_ART_MODULE(protoana::UtilityExample)

