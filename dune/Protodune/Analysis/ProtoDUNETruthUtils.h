#ifndef PROTODUNE_TRUTH_UTILS_H
#define PROTODUNE_TRUTH_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNETruthUtils
//  - Class to help analysers access useful truth information
//
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////


#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Framework/Principal/Event.h"

namespace protoana {

  class ProtoDUNETruthUtils {

  public:

    ProtoDUNETruthUtils();
    ~ProtoDUNETruthUtils();

    const simb::MCParticle* GetMCParticleFromRecoTrack(const recob::Track &track, art::Event const & evt, std::string trackModule) const;
    const recob::Track* GetRecoTrackFromMCParticle(const simb::MCParticle &part, art::Event const & evt, std::string trackModule) const;
    const simb::MCParticle* GetMCParticleFromRecoShower(const recob::Shower &shower, art::Event const & evt, std::string showerModule) const;
    const simb::MCParticle* MatchPduneMCtoG4( const simb::MCParticle & pDunePart, const art::Event & evt );
    const simb::MCParticle* GetGeantGoodParticle(const simb::MCTruth &genTruth, const art::Event &evt) const;

    // Converting times in LArSoft can be a bit of a minefield. These functions convert true times in ns
    // to pandora times in ns
    const float ConvertTrueTimeToPandoraTimeNano(const simb::MCParticle &part) const;
    const float ConvertTrueTimeToPandoraTimeNano(const float trueTime) const;
    // Microsecond versions
    const float ConvertTrueTimeToPandoraTimeMicro(const simb::MCParticle &part) const;
    const float ConvertTrueTimeToPandoraTimeMicro(const float trueTime) const;

    // Get interaction process key.
    int GetProcessKey(std::string process);

    // Get the MC truth deposited energy
    double GetDepEnergyMC(const art::Event &evt, geo::GeometryCore const * fGeom, int trackid, int whichview) const;

    // Get first trajectory point in TPC active volume
    int GetFirstTrajectoryPointInTPCActiveVolume(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh);

    // Get MC Particle length in TPC active volume
    double GetMCParticleLengthInTPCActiveVolume(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh);

  private:




  };

}

#endif
