#ifndef PROTODUNE_TRUTH_UTILS_H
#define PROTODUNE_TRUTH_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNETruthUtils
//  - Class to help analysers access useful truth information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////


#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Framework/Principal/Event.h"

namespace protoana {

  class ProtoDUNETruthUtils {

  public:

    ProtoDUNETruthUtils();
    ~ProtoDUNETruthUtils();

    const simb::MCParticle* GetMCParticleFromRecoTrack(const recob::Track &track, art::Event const & evt, std::string trackModule) const;

  private:




  };

}

#endif

