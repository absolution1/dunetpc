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

namespace protoana {

  class ProtoDUNEPFParticleUtils {

  public:

    ProtoDUNEPFParticleUtils();
    ~ProtoDUNEPFParticleUtils();

    /// Get a map of slice index to the PFParticles within it
    std::map<unsigned int,std::vector<recob::PFParticle*>> GetPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const;

    /// Get the cosmic tag(s) from a given PFParticle
    std::vector<anab::CosmicTag> GetPFParticleCosmicTag(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const;
    /// Get the T0(s) from a given PFParticle
    std::vector<anab::T0> GetPFParticleT0(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const;

  private:




  };

}

#endif

