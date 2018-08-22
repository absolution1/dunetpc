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

  std::map<unsigned int, std::vector<recob::PFParticle*>> sliceMap;

  for(unsigned int p = 0; p < pfParticles->size(); ++p){
    recob::PFParticle* particle = const_cast<recob::PFParticle*>(&(pfParticles->at(p)));

    //  Only the primary particles have the slice association
    if(!particle->IsPrimary()) continue;

    unsigned int thisSlice = GetPFParticleSliceIndex(*particle,evt,particleLabel);

    if(thisSlice != 9999){
      sliceMap[thisSlice].push_back(particle);
    }
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

// Try to get the slice tagged as beam
unsigned short protoana::ProtoDUNEPFParticleUtils::GetBeamSlice(art::Event const &evt, const std::string particleLabel) const{

  std::map<unsigned int, std::vector<recob::PFParticle*>> sliceMap = GetPFParticleSliceMap(evt,particleLabel);

  for(auto slice : sliceMap){
    for(auto particle : slice.second){
      if(IsBeamParticle(*particle,evt,particleLabel)){
        return slice.first;
      }
    }
  }

  return 9999;

}

// Returns pointers for the PFParticles in a slice
std::vector<recob::PFParticle*> protoana::ProtoDUNEPFParticleUtils::GetPFParticlesFromSlice(const unsigned short slice, art::Event const &evt, const std::string particleLabel) const{

  std::map<unsigned int, std::vector<recob::PFParticle*>> sliceMap = GetPFParticleSliceMap(evt,particleLabel);
  
  if(sliceMap.find(slice) != sliceMap.end()){
    return sliceMap.at(slice);
  }
  else{
    return std::vector<recob::PFParticle*>();
  }

}

// Return the pointers for the PFParticles in the beam slice
std::vector<recob::PFParticle*> protoana::ProtoDUNEPFParticleUtils::GetPFParticlesFromBeamSlice(art::Event const &evt, const std::string particleLabel) const{

  unsigned short beamSlice = GetBeamSlice(evt,particleLabel);

  return GetPFParticlesFromSlice(beamSlice,evt,particleLabel);
}

// Access the BDT output used to decide if a slice is beam-like or cosmic-like
float protoana::ProtoDUNEPFParticleUtils::GetBeamCosmicScore(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);

  std::string search = "BeamScore";
  if(mdMap.find(search) != mdMap.end()){
    return mdMap.at(search);
  }
  else{
//    std::cerr << "Object has no beam score... returning -999." << std::endl;
    return -999.;
  }
}

// Use the pandora metadata to tell us if this is a beam particle or not
bool protoana::ProtoDUNEPFParticleUtils::IsBeamParticle(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);

  // IsTestBeam only appears for beam particles
  std::string search = "IsTestBeam";
  if(mdMap.find(search) != mdMap.end()){
    return true;
  }
  else{
    return false;
  }

}

// Get the reconstructed slice associated with a particle
unsigned short protoana::ProtoDUNEPFParticleUtils::GetPFParticleSliceIndex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);

  std::string search = "SliceIndex";
  if(mdMap.find(search) != mdMap.end()){
    return static_cast<unsigned short>(mdMap.at(search));
  }
  else{
//    std::cerr << "Object has no slice index... returning 9999" << std::endl;
    return 9999;
  }
}

const std::map<std::string,float> protoana::ProtoDUNEPFParticleUtils::GetPFParticleMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const {

  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  // And their meta data
  const art::FindManyP<larpandoraobj::PFParticleMetadata> findMetaData(pfParticles,evt,particleLabel);

  const larpandoraobj::PFParticleMetadata metaData = *((findMetaData.at(particle.Self())).at(0));

  return metaData.GetPropertiesMap();
}

// Pandora tags and removes clear cosmics before slicing, so check if this particle is a clear cosmic
bool protoana::ProtoDUNEPFParticleUtils::IsClearCosmic(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);

  // IsTestBeam only appears for beam particles
  std::string search = "IsClearCosmic";
  if(mdMap.find(search) != mdMap.end()){
    return true;
  }
  else{
    return false;
  }

}

// Get all of the clear cosmic ray particles
std::vector<recob::PFParticle*> protoana::ProtoDUNEPFParticleUtils::GetClearCosmicPFParticles(art::Event const &evt, const std::string particleLabel) const{

  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);

  std::vector<recob::PFParticle*> cosmicParticles;

  for(unsigned int p = 0; p < pfParticles->size(); ++p){
    recob::PFParticle* particle = const_cast<recob::PFParticle*>(&(pfParticles->at(p)));

    //  Only consider primary particles
    if(!particle->IsPrimary()) continue;

    if(IsClearCosmic(*particle,evt,particleLabel)){
      cosmicParticles.push_back(particle);
    }

  }

  return cosmicParticles;

}


