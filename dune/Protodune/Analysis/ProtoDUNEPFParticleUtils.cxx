#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

protoana::ProtoDUNEPFParticleUtils::ProtoDUNEPFParticleUtils(){

}

protoana::ProtoDUNEPFParticleUtils::~ProtoDUNEPFParticleUtils(){

}

// Return a map of particles grouped by their reconstructed slice. Useful for finding slices with multiple particles
std::map<unsigned int,std::vector<recob::PFParticle*>> protoana::ProtoDUNEPFParticleUtils::GetPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const{

  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  // And their meta data
  const art::FindManyP<larpandoraobj::PFParticleMetadata> findMetaData(pfParticles,evt,particleLabel);

  std::map<unsigned int, std::vector<recob::PFParticle*>> sliceMap;

  for(unsigned int p = 0; p < pfParticles->size(); ++p){
    recob::PFParticle* particle = const_cast<recob::PFParticle*>(&(pfParticles->at(p)));

    //  Only the primary particles have the slice association
    if(!particle->IsPrimary()) continue;

    // If there is no metadata then just move on
    if(findMetaData.at(p).size() == 0) continue;

    const larpandoraobj::PFParticleMetadata metaData = *((findMetaData.at(p)).at(0));

    // If there is no SliceIndex entry in the meta data then carry on
    if(metaData.GetPropertiesMap().find("SliceIndex") == metaData.GetPropertiesMap().end()) continue;

    unsigned int thisSlice = static_cast<unsigned int>(metaData.GetPropertiesMap().at("SliceIndex"));
  
    sliceMap[thisSlice].push_back(particle);

  }

  return sliceMap;

}

// Get the cosmic tag(s) from a given PFParticle
std::vector<anab::CosmicTag> protoana::ProtoDUNEPFParticleUtils::GetPFParticleCosmicTag(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const{

  try{
    auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle> >(particleLabel);
  
    unsigned int pIndex = particle.Self();
  
    std::vector<anab::CosmicTag> pTags;
  
    const art::FindManyP<anab::CosmicTag> findCosmicTags(recoParticles,evt,particleLabel);
    for(unsigned int p = 0; p < findCosmicTags.at(pIndex).size(); ++p){
      pTags.push_back((*(findCosmicTags.at(pIndex)[p])));
    }

    return pTags;
  }
  catch(...){
    return std::vector<anab::CosmicTag>();
  }
}

// Get the T0(s) from a given PFParticle
std::vector<anab::T0> protoana::ProtoDUNEPFParticleUtils::GetPFParticleT0(const recob::PFParticle &particle, art::Event const &evt, std::string particleLabel) const{

  try{
    auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle> >(particleLabel);
  
    unsigned int pIndex = particle.Self();
  
    std::vector<anab::T0> pT0s;
  
    const art::FindManyP<anab::T0> findParticleT0s(recoParticles,evt,particleLabel);
    for(unsigned int p = 0; p < findParticleT0s.at(pIndex).size(); ++p){
      pT0s.push_back((*(findParticleT0s.at(pIndex)[p])));
    }

    return pT0s;
  }
  catch(...){
    return std::vector<anab::T0>();
  }

}


