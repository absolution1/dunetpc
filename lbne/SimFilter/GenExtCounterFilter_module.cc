#include <iostream>
#include <utility>
#include <complex>
#include <algorithm>


#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "Geometry/Geometry.h"
#include "SimulationBase/MCTruth.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"

namespace filt{

  class GenFilter : public art::EDFilter {
    public:
      explicit GenFilter(fhicl::ParameterSet const & pset) {};
      virtual ~GenFilter() {};
      virtual bool filter(art::Event& e);
      void reconfigure(fhicl::ParameterSet const& pset) {};
      void beginJob() ;

    private:

      struct CounterSetPair{
        bool isRequested;
        std::vector<geo::AuxDetGeo*> setA;
        std::vector<geo::AuxDetGeo*> setB;
        double normalVec[3];

      };
      std::vector<CounterSetPair> fCounterSetPairs;

      bool IsInterestingParticle(const simb::MCParticle &particle);
      bool ParticleHitsCounterSetPairs(const simb::MCParticle &particle, const CounterSetPair &CSP);
      bool ParticleHitsCounterSet(const simb::MCParticle &particle, const std::vector<geo::AuxDetGeo*> &counters, const TVector3 &counter_norm);
  
  };


  bool GenFilter::filter(art::Event & e){

    std::vector< art::Handle< std::vector<simb::MCTruth> > > mclists;
    e.getManyByType(mclists);
    for (unsigned int i = 0; i < mclists.size() ; i++){
      for (unsigned int j = 0; j < mclists[i]->size(); j++){
        //Should have the truth record for the event now
        const art::Ptr<simb::MCTruth> mc_truth(mclists[i],j);
        for (int part = 0; part < mc_truth->NParticles(); part++){
          const simb::MCParticle particle = mc_truth->GetParticle(part);
          if (!IsInterestingParticle(particle)) continue;
          TVector3 particle_pos = particle.Position().Vect();
          TVector3 particle_dir = particle.Momentum().Vect().Unit();
          for (unsigned CSP_i = 0; CSP_i < fCounterSetPairs.size(); CSP_i++){
            if (!fCounterSetPairs[CSP_i].isRequested) continue;
            if (ParticleHitsCounterSetPairs(particle,fCounterSetPairs[i])){
              std::cout<<"HIT COUNTER SET"<<std::endl;
              return true;
            }
            /*
            TVector3 counter_norm(fCounterSetPairs[i].normalVec[0],fCounterSetPairs[i].normalVec[1],fCounterSetPairs[i].normalVec[2]);
            double counter_pos_array[3];
            fCounterSetPairs[CSP_i].setA.front()->GetCenter(counter_pos_array);
            TVector3 counter_pos(counter_pos_array[0], counter_pos_array[1], counter_pos_array[2]);
            double scale_factor = counter_norm.Dot(counter_pos - particle_pos)/(counter_norm.Dot(particle_dir));
            TVector3 particle_pos_in_plane = particle_pos + scale_factor * particle_dir;
            std::cout<<"Particle pos: " << "("<<particle_pos_in_plane.X()<<","<<particle_pos_in_plane.Y()<<","<<particle_pos_in_plane.Z()<<")"<<std::endl;
            std::cout<<"Counter pos: " << "("<<counter_pos.X()<<","<<counter_pos.Y()<<","<<counter_pos.Z()<<")"<<std::endl;
            */

          }
          /*
          geo::AuxDetGeo *counter = fCounterSetPairs.front().setA.front();
          std::cout<<"Length: " << counter->Length() << std::endl;
          std::cout<<"HalfWidth1: " << counter->HalfWidth1() << std::endl;
          std::cout<<"HalfWidth2: " << counter->HalfWidth2() << std::endl;

          std::cout<<"HalfHeight: " << counter->HalfHeight() << std::endl;
          */
        }
      }
    }

    return false;
  }

  void GenFilter::beginJob() {
    art::ServiceHandle<geo::Geometry> geom;

    CounterSetPair fEWCounterSetPair;
    for (unsigned int i = 0; i < geom->AuxDetGeoVec().size(); i++){

      if (i >= 6 && i <= 15) fEWCounterSetPair.setA.push_back(geom->AuxDetGeoVec()[i]);
      else if (i >= 28 && i <= 37) fEWCounterSetPair.setB.push_back(geom->AuxDetGeoVec()[i]);
    }
    fEWCounterSetPair.normalVec[0] = 0.;
    fEWCounterSetPair.normalVec[1] = 0.;
    fEWCounterSetPair.normalVec[2] = 1.;
    fEWCounterSetPair.isRequested = true;
    fCounterSetPairs.push_back(fEWCounterSetPair);

  }

  bool GenFilter::IsInterestingParticle(const simb::MCParticle &particle){

    if (particle.PdgCode() == 13) return true;
    return false;
  }

  bool GenFilter::ParticleHitsCounterSetPairs(const simb::MCParticle &particle, const CounterSetPair &CSP){
    //Need the normal to the counters
    TVector3 counter_norm(CSP.normalVec[0],CSP.normalVec[1],CSP.normalVec[2]);

    //Loop through one of the counter sets
    if (ParticleHitsCounterSet(particle, CSP.setA, counter_norm)){
      if (ParticleHitsCounterSet(particle, CSP.setB, counter_norm)){
        return true;
      }
    }
    return false;
  }

  bool GenFilter::ParticleHitsCounterSet(const simb::MCParticle &particle, const std::vector<geo::AuxDetGeo*> &counters, const TVector3 &counter_norm){

    //Loop through the counters
    for (unsigned int i = 0; i < counters.size(); i++){
      //First step is to push the particle to the counter plane
      geo::AuxDetGeo *counter = counters[i];

      double counter_pos_array[3];
      //fCounterSetPairs[CSP_i].setA.front()->GetCenter(counter_pos_array);
      counter->GetCenter(counter_pos_array);
      TVector3 counter_pos(counter_pos_array[0], counter_pos_array[1], counter_pos_array[2]);
      TVector3 particle_pos = particle.Position().Vect();
      TVector3 particle_dir = particle.Momentum().Vect().Unit();

      /*
          Length: 62.992
          HalfWidth1: 16.2814
          HalfWidth2: 13.5255
          HalfHeight: 0.475
          */

      double scale_factor = counter_norm.Dot(counter_pos - particle_pos)/(counter_norm.Dot(particle_dir));

      TVector3 particle_pos_in_plane = particle_pos + scale_factor * particle_dir;

      //We now have the particle position in the plane of the counter.  We now just need to calculate whether it sits inside the box
      //Create two TVector3s, each will hold the coordinates of opposing corners of the counter in question
      //We need to use the normals associated with the counter.  We already have two, one is a member of the C++ object and the other was passed to this function
      //Create the two TVector3s
      TVector3 pos_corner, neg_corner;
      //Lets start with the one passed to this function
      //The thin dimension of the counter is the one associated with normal passed to this function
      pos_corner += counter->HalfHeight()*counter_norm;
      neg_corner += -1.*counter->HalfHeight()*counter_norm; 

      //Now lets to the same for the normal stored in the C++ object (this "normal" actually points out the top of a counter)
      double counter_top_norm_array[3];
      counter->GetNormalVector(counter_top_norm_array);
      //Now package this up a TVector
      TVector3 counter_top_norm;
      counter_top_norm.SetX(counter_top_norm_array[0]);
      counter_top_norm.SetY(counter_top_norm_array[1]);
      counter_top_norm.SetZ(counter_top_norm_array[2]);
      //OK now we can add the dimension to the corner vectors.  The relevant dimension this time is Length/2
      pos_corner += counter->Length()*counter_top_norm*0.5;
      neg_corner += -1.*counter->Length()*counter_top_norm*0.5;
      //Now we need to the same for the vector pointing along the counter.  
      //The relevant dimension in this case in HalfWidth1 (going to assume the counter is a square and take the bigger width)
      //Because we have the other two vectors already, we can very easily get the final one by taking the cross product of them
      TVector3 counter_side_norm = counter_norm.Cross(counter_top_norm); 
      //now add the dimenions onto the corner vectors
      pos_corner += counter->HalfWidth1()*counter_side_norm;
      neg_corner += -1.*counter->HalfWidth1()*counter_side_norm;

      //Almost ready
      //We have calculated the corners assuming the centre of the counter is the origin.  Two choices, either translate the corners OR translate the propagated particle position
      //The latter is less lines of code so lets to that
      particle_pos_in_plane += -1.*counter_pos;

      //We are now ready to check if the particle sits in the counter
      //We don't know by default which way round the corners are oriented, but that doesn't matter as we can take abs, max and min
      //We are going to take abs of the particle position, the pos corner and the neg corner first
      /*
      pos_corner.SetXYZ(std::abs(pos_corner.X()),std::abs(pos_corner.Y()),std::abs(pos_corner.Z()));
      neg_corner.SetXYZ(std::abs(neg_corner.X()),std::abs(neg_corner.Y()),std::abs(neg_corner.Z()));
      particle_pos_in_plane.SetXYZ(std::abs(particle_pos_in_plane.X()),std::abs(particle_pos_in_plane.Y()),std::abs(particle_pos_in_plane.Z()));
      */

      /*
      //Dump the corners and particle pos
      std::cout<<"pos_corner: " << pos_corner.X()<<","<<pos_corner.Y()<<","<<pos_corner.Z()<<std::endl;
      std::cout<<"neg_corner: " << neg_corner.X()<<","<<neg_corner.Y()<<","<<neg_corner.Z()<<std::endl;
      std::cout<<"particle_pos_in_plane: " << particle_pos_in_plane.X()<<","<<particle_pos_in_plane.Y()<<","<<particle_pos_in_plane.Z()<<std::endl;
      */



      //Now we can check
      //This is going to be one huge if statement...
      if (particle_pos_in_plane.X() > std::min(pos_corner.X(), neg_corner.X()) && particle_pos_in_plane.X() < std::max(pos_corner.X(),neg_corner.X()) && particle_pos_in_plane.Y() > std::min(pos_corner.Y(), neg_corner.Y()) && particle_pos_in_plane.Y() < std::max(pos_corner.Y(),neg_corner.Y()) && particle_pos_in_plane.Z() > std::min(pos_corner.Z(), neg_corner.Z()) && particle_pos_in_plane.Z() < std::max(pos_corner.Z(),neg_corner.Z())){
        std::cout<<"Particle uses counter in set"<<std::endl;
        return true;
      }

    }

    return false;
  }


  DEFINE_ART_MODULE(GenFilter)

}
