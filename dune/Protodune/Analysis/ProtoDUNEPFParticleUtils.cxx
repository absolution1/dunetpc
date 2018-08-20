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

