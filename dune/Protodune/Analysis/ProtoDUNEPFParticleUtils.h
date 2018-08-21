#ifndef PROTODUNE_PFPARTICLE_UTILS_H
#define PROTODUNE_PFPARTICLE_UTILS_H

///////////////////////////////////////////////////////////////////
// ProtoDUNEPFParticleUtils
//  - Class to help analysers access useful PFParticle information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////////

#include <map>
#include <string>

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "art/Framework/Principal/Event.h"
#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

namespace protoana {

  class ProtoDUNEPFParticleUtils {

  public:

    ProtoDUNEPFParticleUtils();
    ~ProtoDUNEPFParticleUtils();

    /// Get a map of slice index to the PFParticles within it
    std::map<unsigned int,std::vector<recob::PFParticle*>> GetPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const;

    /// Try to get the slice tagged as beam
    unsigned short GetBeamSlice(art::Event const &evt, const std::string particleLabel) const;

    /// Return the pointers for the PFParticles in the beam slice
    std::vector<recob::PFParticle*> GetPFParticlesFromBeamSlice(art::Event const &evt, const std::string particleLabel) const;

    /// Get the cosmic tag(s) from a given PFParticle
    std::vector<anab::CosmicTag> GetPFParticleCosmicTag(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const;
    /// Get the T0(s) from a given PFParticle
    std::vector<anab::T0> GetPFParticleT0(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const;

    /// Access the BDT output used to decide if a slice is beam-like or cosmic-like
    float GetBeamCosmicScore(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Use the pandora metadata to tell us if this is a beam particle or not
    bool IsBeamParticle(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the reconstructed slice associated with a particle
    unsigned short GetPFParticleSliceIndex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the metadata associated to a PFParticle from pandora
    const std::map<std::string,float> GetPFParticleMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

  private:

  };

}

#endif

