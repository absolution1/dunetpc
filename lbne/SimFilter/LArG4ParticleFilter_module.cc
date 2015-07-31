#include <iostream>
#include <utility>
#include <set>


#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "Geometry/Geometry.h"
#include "Simulation/AuxDetSimChannel.h"
#include "SimulationBase/MCTruth.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"

namespace filt{

  class LArG4ParticleFilter : public art::EDFilter {
    public:
      explicit LArG4ParticleFilter(fhicl::ParameterSet const & pset);
      virtual ~LArG4ParticleFilter() {};
      virtual bool filter(art::Event& e);
      void reconfigure(fhicl::ParameterSet const& pset);
      void beginJob() ;

    private:

      bool IsInterestingParticle(const art::Ptr<simb::MCParticle> particle);

      //internal variables
      std::vector<int> fInterestingPDGs;

  };

  LArG4ParticleFilter::LArG4ParticleFilter::LArG4ParticleFilter(fhicl::ParameterSet const & pset)
  {
    this->reconfigure(pset);
  }

  void LArG4ParticleFilter::reconfigure(fhicl::ParameterSet const& pset){
    fInterestingPDGs = pset.get<std::vector<int> >("InterestingPDGs");
    std::cout<<"NInteresting PDGs: " << fInterestingPDGs.size() << std::endl;
    for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
      std::cout<<"-- PDG: " << fInterestingPDGs[i] << std::endl; 
    }
    return;
  }

  bool LArG4ParticleFilter::filter(art::Event & e){

    //art::ServiceHandle<geo::Geometry> geom;


    //Get the vector of particles
    art::Handle<std::vector<simb::MCParticle> > particles;
    e.getByLabel("largeant",particles);
    //Loop through the particles
    for (unsigned int part_i = 0; part_i < particles->size(); part_i++){
      //Get the particle
      const art::Ptr<simb::MCParticle> particle(particles,part_i);
      //Check if the particle is initially worth considering
      if (IsInterestingParticle(particle)){
        return true;
      }
    }

    //Assume that the event is not worth saving
    return false;
  }

  void LArG4ParticleFilter::beginJob() {

    return;
  }

  bool LArG4ParticleFilter::IsInterestingParticle(const art::Ptr<simb::MCParticle> particle){
    //The bulk of the code goes here
    //Each check should probably be its own function, but for now everything can be crunched inline
    bool OK = false;

    //Check the particle PDG first
    int pdg = particle->PdgCode();
    //Loop through the PDG vector and see if we have a match
    for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
      if (pdg == fInterestingPDGs[i]){
        OK = true;
        break;
      }
    }
    if (!OK) return false;

    //Che

    return true;;
  }

  DEFINE_ART_MODULE(LArG4ParticleFilter)
}
