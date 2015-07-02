#include <iostream>
#include <utility>


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

  class LArG4Filter : public art::EDFilter {
    public:
      explicit LArG4Filter(fhicl::ParameterSet const & pset);
      virtual ~LArG4Filter() {};
      virtual bool filter(art::Event& e);
      void reconfigure(fhicl::ParameterSet const& pset);
      void beginJob() ;

    private:

      //Small class to hold the opposing set of counters
      struct CounterSetPair{
        std::vector<unsigned int> setA; //First set
        std::vector<unsigned int> setB; //Second set
        bool IsRequested; //Check if we care about a particular pairwise set of counters (e.g. does the user care about the EW counters etc)
      };
      std::vector<CounterSetPair> fCounterSetPairs;

      //Adjustable parameters
      bool fUseEWCounterPair; //Use the EW counter pair
      bool fUseNupSdownCounterPair; //Use North (up) South (down) counter pair
      bool fUseNdownSupCounterPair; //Use the North (down) South (up) counter pair
      std::vector<int> fInterestingPDGs; //A vector of particle PDGs which want to be filtered on
      double fParticleMinEnergy;  //The minimum energy of a particle to be filtered
      double fParticleMaxEnergy;  //The maximum energy of a particle to be filtered
      double fParticleMinEDep;  //The minimum energy a particle must deposit in an external counter to be filtered
  
      bool IsInterestingParticle(const simb::MCParticle &particle);  //Define whether a particular particle is initially worth saving e.g. is it a muon, pion etc
      bool UsesCounterSetPair(const CounterSetPair &CSP, const std::vector<const sim::AuxDetSimChannel*> &auxDetSimChannels, const simb::MCParticle &particle); //Check if a particle uses a pariwse set of counters
      bool UsesCounterSet(const std::vector<unsigned int> &counterIDs, const std::vector<const sim::AuxDetSimChannel*> &auxDetSimChannels, const simb::MCParticle &particle); //Check if a particle uses one set of the pairwise set of counters

  };

  LArG4Filter::LArG4Filter::LArG4Filter(fhicl::ParameterSet const & pset)
  {
    this->reconfigure(pset);
  }

  void LArG4Filter::reconfigure(fhicl::ParameterSet const& pset){
    fUseEWCounterPair = pset.get<bool>("UseEWCounterPair",1);
    std::cout<<"Use EW counter pair: " << fUseEWCounterPair<<std::endl;
    fUseNupSdownCounterPair = pset.get<bool>("UseNupSdownCounterPair",1);
    std::cout<<"Use N (up) S (down) counter pair: " << fUseNupSdownCounterPair << std::endl;
    fUseNdownSupCounterPair = pset.get<bool>("UseNdownSupCounterPair",1);
    std::cout<<"Use N (down) S (up) counter pair: " << fUseNdownSupCounterPair << std::endl;
    fInterestingPDGs = pset.get<std::vector<int> >("InterestingPDGs");
    std::cout<<"NInteresting PDGs: " << fInterestingPDGs.size() << std::endl;
    for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
      std::cout<<"-- PDG: " << fInterestingPDGs[i] << std::endl; 
    }
    fParticleMinEnergy = pset.get<double>("ParticleMinEnergy",0);
    std::cout<<"Particle min energy: " << fParticleMinEnergy << std::endl;
    fParticleMaxEnergy = pset.get<double>("ParticleMaxEnergy",999999999);
    std::cout<<"Particle max energy: " << fParticleMaxEnergy << std::endl;
    fParticleMinEDep = pset.get<double>("ParticleMinEDep",1e-9);
    std::cout<<"Particle min E deposit: " << fParticleMinEDep << std::endl;
  }

  bool LArG4Filter::filter(art::Event & e){

    //Get the simulated counters
    std::vector<const sim::AuxDetSimChannel*> auxDetSimChannels;
    e.getView("largeant",auxDetSimChannels);

    //Get the simulated particles
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mclists;
    e.getManyByType(mclists);
    //Loop through the list of particles
    for (unsigned int i = 0; i < mclists.size() ; i++){
      for (unsigned int j = 0; j < mclists[i]->size(); j++){
        const art::Ptr<simb::MCTruth> mc_truth(mclists[i],j);
        for (int part = 0; part < mc_truth->NParticles(); part++){
          //Get an individual particle
          const simb::MCParticle particle = mc_truth->GetParticle(part);
          //Check if the particle is initially interested
          if (IsInterestingParticle(particle) == true){
            //Now see if it goes through a pairwise counter set
            for (unsigned int pair_i = 0; pair_i < fCounterSetPairs.size(); pair_i++){
              if (!fCounterSetPairs[pair_i].IsRequested) continue;
              //If the particle uses deposits energy in both of the opposing counter sets, save it
              if (UsesCounterSetPair(fCounterSetPairs[pair_i], auxDetSimChannels, particle)){
                //std::cout<<"DID USE COUNTER"<<std::endl;
                std::cout<<"SELECTED EVENT"<<std::endl;
                return true;
              }
              //else std::cout<<"DID NOT USE COUNTER"<<std::endl;
            }
          }
        }
      }
    }

    //Assume that the event is not worth saving
    return false;
  }

  void LArG4Filter::beginJob() {

    //Need to get the counter information.  By doing this at the start of the job, rather than per event, the code assumes the geomtry is not going to change between events
    art::ServiceHandle<geo::Geometry> geom;
    //Create the pairwise counter sets
    CounterSetPair EWCounterSetPair;
    CounterSetPair NupSdownCounterSetPair;
    CounterSetPair NdownSupCounterSetPair;
    //A stupid way of storing the IDs of the counters, this REALLY needs changing
    //The code loops through all of the counters in the geomtry, and if the number matches a particular counter number e.g. one of the east counters, store it in the correct, pairwise set
    for (unsigned int i = 0; i < geom->AuxDetGeoVec().size(); i++){
      //The WE counter pairs
      if (i >=6 && i <= 15) EWCounterSetPair.setA.push_back(i);
      else if (i >= 28 && i <=37) EWCounterSetPair.setB.push_back(i);
      //The N (up) S (down) counter pairs
      else if (i >= 22 && i <= 27) NupSdownCounterSetPair.setA.push_back(i);
      else if (i >= 0 && i <= 5) NupSdownCounterSetPair.setB.push_back(i);
      //The N (down) S (up) counter pairs
      else if (i >= 16 && i <= 21) NdownSupCounterSetPair.setA.push_back(i);
      else if (i >= 38 && i <= 43) NdownSupCounterSetPair.setB.push_back(i);
    }

    //Enable/disable the counter set pairs
    EWCounterSetPair.IsRequested = fUseEWCounterPair;
    NupSdownCounterSetPair.IsRequested = fUseNupSdownCounterPair;
    NdownSupCounterSetPair.IsRequested = fUseNdownSupCounterPair;

    //Store them onto the vector of counter set pairs
    fCounterSetPairs.push_back(EWCounterSetPair);
    fCounterSetPairs.push_back(NupSdownCounterSetPair);
    fCounterSetPairs.push_back(NdownSupCounterSetPair);


    for (unsigned int i = 0; i < fCounterSetPairs.size(); i++){
      std::cout<<"Counter set pair: "<<i<<std::endl;
      std::cout<<"--setA size: " << fCounterSetPairs[i].setA.size() << std::endl;
      std::cout<<"--setB size: " << fCounterSetPairs[i].setB.size() << std::endl;

    }
  }

  bool LArG4Filter::IsInterestingParticle(const simb::MCParticle &particle){
    //Loop over the list of requested PDGs.  See if that matches the particle under consideration
    for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
      //Check if the particle under consideration has a requested PDG
      if (particle.PdgCode() == fInterestingPDGs[i]){
        //Got a requested PDG,  now check that the energy matches the requested range
        TLorentzVector mom4 = particle.Momentum();
        if (mom4.T() > fParticleMinEnergy && mom4.T() < fParticleMaxEnergy){
          std::cout<<"FOUND INTERESTING PARTICLE"<<std::endl;
          return true;
        }
      }
    }
    return false;
  }

  bool LArG4Filter::UsesCounterSetPair(const CounterSetPair &CSP, const std::vector<const sim::AuxDetSimChannel*> &auxDetSimChannels, const simb::MCParticle &particle){
    //Check if particle uses one of the counters in the first pair
    if (UsesCounterSet(CSP.setA, auxDetSimChannels, particle)){
      //OK, it used the first set of the pair.  Now see if it uses the second set
      if (UsesCounterSet(CSP.setB, auxDetSimChannels, particle)) return true;
    }
    return false;
  }

  bool LArG4Filter::UsesCounterSet(const std::vector<unsigned int> &counterIDs, const std::vector<const sim::AuxDetSimChannel*> &auxDetSimChannels, const simb::MCParticle &particle){
    //The bulk of the filter code is here.
    //This code works by looking at which particle deposited energy in a particular counter.  If the ID of the energy deposit matches the particle passed as an argument, then the requested particle used one of the counters 
    //Loop through the counter IDs in the set
    for (unsigned int id_i = 0; id_i < counterIDs.size(); id_i++){
      //Get the ID of the counter
      int counterID = counterIDs[id_i];
      //Now we need to get the corresponding AuxDetSimChannel
      const sim::AuxDetSimChannel *channel = auxDetSimChannels[counterID];
      //Now need the vector of IDEs for this channel
      const std::vector<sim::AuxDetIDE>& auxDetIDEs = channel->AuxDetIDEs();
      //Loop over them
      for (unsigned int ide_i = 0; ide_i < auxDetIDEs.size(); ide_i++){
        //Check to see if the particle resposible for the deposited energy is the same as the one passed to this function
        if (auxDetIDEs[ide_i].trackID == particle.TrackId() && auxDetIDEs[ide_i].energyDeposited > fParticleMinEDep) {
          return true;
        }
      }
    }

    //Assume the particle did not deposit energy in any of the counters in the set
    return false;
  }

  DEFINE_ART_MODULE(LArG4Filter)
}
