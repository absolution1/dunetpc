#ifndef PROTODUNE_DATA_UTILS_H
#define PROTODUNE_DATA_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNEDataUtils
//  - Class to help analysers access useful beam data 
//    information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include <set>

namespace protoana {

  struct PossibleParticleCands {
    public:
      bool electron = false;
      bool muon = false;
      bool pion = false;
      bool kaon = false;
      bool proton = false;

    inline PossibleParticleCands operator&&(const PossibleParticleCands & b) const {
        return {electron && b.electron, muon && b.muon, pion && b.pion, kaon && b.kaon, proton && b.proton};
    }
    inline PossibleParticleCands operator||(const PossibleParticleCands & b) const {
        return {electron || b.electron, muon || b.muon, pion || b.pion, kaon || b.kaon, proton || b.proton};
    }
  };

  class ProtoDUNEDataUtils {

  public:

    ProtoDUNEDataUtils(fhicl::ParameterSet const &pset);
    ~ProtoDUNEDataUtils();

    void reconfigure(fhicl::ParameterSet const &pset);

    /**
     * Returns true if the ProtoDUNE trigger says this is a beam trigger
     */
    bool IsBeamTrigger(art::Event const & evt) const;

    /// Get number of active fembs in an APA
    int GetNActiveFembsForAPA(art::Event const & evt, int apa) const;

    /**
     * Returns true if the beamline instrumentation has a good trigger
     * that matches the ProtoDUNE trigger.
     */
    bool IsGoodBeamlineTrigger(art::Event const & evt) const;

    /**
     * Uses the beamline momentum and time of flight measurements
     * to estimate the beamline particle mass in GeV/c^2
     * Returns a list of all possible cominations.
     */
    std::vector<double> GetBeamlineMass(art::Event const & evt) const;

    /**
     * Based on the beam energy of the given run and whether each of the Cherenkov
     * detectors fired, sets the particle types this could possibly be to true in the
     * returned PossibleParticleCands object
     *
     */
    PossibleParticleCands GetCherenkovParticleID(art::Event const & evt, const float beamEnergyGeV) const;
    /**
     * Based on the energy of a given run, the beamline momentum, 
     * and time of flight, sets the particle types this could 
     * possibly be to true in the returned PossibleParticleCands 
     * object
     */
    PossibleParticleCands GetTOFParticleID(art::Event const & evt, const float beamEnergyGeV) const;
    /**
     * Combines TOF, momentum, and Cherenkov info to determine a particle ID.
     * Sets the particle types this could possibly be to true in 
     * the returned PossibleParticleCands object
     */
    PossibleParticleCands GetBeamlineParticleID(art::Event const & evt, const float beamEnergyGeV) const;

    /**
     *  Get reconstructed beamline momentum (in GeV/c), tof (in ns), and
     *  flags for the tofChannel and ckov's. Will be < 0 if invalid.
     *
     *  C++ structured binding, so call like:
     *  const auto [momentum, tof, tofChannel,ckov0,ckov1] = dataUtils.GetBeamlineInformation(e);
     *
     *  then you have the normal float momentum, int ckov0 variables, etc. in the current scope.
     */
    const std::tuple<double,double,int,int,int> GetBeamlineVars(art::Event const & evt) const;

    /**
     *  Get reconstructed beamline momentum (in GeV/c), tof (in ns),
     *  flags for the tofChannel and ckov's. Also the timing trigger (12 means beam),
     *  the beam instrumentation trigger (1 means beam trigger), and whether
     *  the timing and beam instrumentation triggers are matched.
     *
     *  All values will be < 0 if invalid.
     *
     *  C++ structured binding, so call like:
     *  const auto [momentum, tof, tofChannel,ckov0,ckov1,timingTrigger,BITrigger,areBIAndTimingMatched] = dataUtils.GetBeamlineInformation(e);
     *
     *  then you have the normal float momentum, int ckov0 variables, etc. in the current scope.
     */
    const std::tuple<double,double,int,int,int,int,int,bool> GetBeamlineVarsAndStatus(art::Event const & evt) const;

  private:

    art::InputTag fTimingTag;
    art::InputTag fBeamEventTag;

    bool fStrictTOF;
    bool fStrictCherenkov;

    float fTOFDistance; // meters
    float fMomentumScaleFactor;
    float fMomentumOffset; // GeV/c
    float fTOFScaleFactor;
    float fTOFOffset; // ns

  };

}

#endif

