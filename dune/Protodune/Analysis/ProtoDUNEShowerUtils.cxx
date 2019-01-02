#include "dune/Protodune/Analysis/ProtoDUNEShowerUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

protoana::ProtoDUNEShowerUtils::ProtoDUNEShowerUtils(){

}

protoana::ProtoDUNEShowerUtils::~ProtoDUNEShowerUtils(){

}

// Get the hits from a given reco shower
const std::vector<const recob::Hit*> protoana::ProtoDUNEShowerUtils::GetRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);
  art::FindManyP<recob::Hit> findHits(recoShowers,evt,showerModule);

  // Shower.ID is sometimes at a default value - make sure we get the correct one
  int actualIndex = GetShowerIndex(shower,evt,showerModule);

  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(actualIndex);
  std::vector<const recob::Hit*> showerHits;

  for(const art::Ptr<recob::Hit> hit : inputHits){

    showerHits.push_back(hit.get());

  }

  return showerHits;  

}

// Get the hits from a given reco shower
unsigned int protoana::ProtoDUNEShowerUtils::GetNumberRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  return GetRecoShowerHits(shower,evt,showerModule).size();

}

// Get the PCAxis object from the reco shower
std::vector<const recob::PCAxis*> protoana::ProtoDUNEShowerUtils::GetRecoShowerPCAxis(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);
  art::FindManyP<recob::PCAxis> findPCA(recoShowers,evt,showerModule);

  // Shower.ID is sometimes at a default value - make sure we get the correct one
  int actualIndex = GetShowerIndex(shower,evt,showerModule);

  std::vector<const recob::PCAxis*> pcaVec;
  for(auto const pca : findPCA.at(actualIndex)){
    pcaVec.push_back(pca.get());
  }

  return pcaVec;
}

// If the shower.ID() isn't filled we must find the actual shower index ourselves
int protoana::ProtoDUNEShowerUtils::GetShowerIndex(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) const{

  if(shower.ID() != -999) return shower.ID();

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);

  // Iterate through all showers to find the matching one to our shower
  int actualIndex = shower.ID();
  if(shower.ID() < 0){
    for(unsigned int s = 0; s < recoShowers->size(); ++s){
      const recob::Shower thisShower = (*recoShowers)[s];
      // Can't compare actual objects so look at a property
      if(fabs(thisShower.Length() - shower.Length()) < 1e-5){
        actualIndex = s;
        continue;
      }
    }
  }

  return actualIndex;

}




