////////////////////////////////////////////////////////////////////////
//
//  Very simple module to filter events without all of the expected
//  products in the event.
//
//  Leigh Whitehead - leigh.howard.whitehead@cern.ch
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace mpf{
  class MissingProductFilter;
}

class mpf::MissingProductFilter : public art::EDFilter{
public:
  
  explicit MissingProductFilter(fhicl::ParameterSet const& pset);
  virtual ~MissingProductFilter();
  
  void beginJob() override;
  bool filter(art::Event& evt) override;
  bool beginRun(art::Run& r) override;
  void endJob() override;

private:

  std::vector<std::string> fModules;
};
  
//-----------------------------------------------------------------------
mpf::MissingProductFilter::MissingProductFilter(fhicl::ParameterSet const& pset):
  EDFilter(pset)
{
  fModules = pset.get<std::vector<std::string>>("HitModule");
}

//-----------------------------------------------------------------------
mpf::MissingProductFilter::~MissingProductFilter(){}

//-----------------------------------------------------------------------
void mpf::MissingProductFilter::beginJob() {}

//-----------------------------------------------------------------------
bool mpf::MissingProductFilter::beginRun(art::Run& r){
  if (fScaleThresholdForReadoutWindow){
    unsigned int fSize = art::ServiceHandle<detinfo::DetectorPropertiesService const>{}->provider()->ReadOutWindowSize();
    fHitLimit = (unsigned int)(fHitLimit*fSize/6000.);
    std::cout<<"Scale HitLimit based on readout window size "<<fSize<<std::endl;
    std::cout<<"HitLimit = "<<fHitLimit<<std::endl;
  }
  return true;
}

//-----------------------------------------------------------------------
bool mpf::MissingProductFilter::filter(art::Event& evt){

  for(const std::string &s : fModules){

    try{
      auto product = evt.get
    }

  }


  // Get the mpf collection from the event 
  auto allHits = evt.getValidHandle<std::vector<recob::Hit> >(fHitModule);  
 
  bool result = true;

  if(fLimitPerTPC){
    // Find the number of mpfs per TPC and then filter based on a large value
    std::map<unsigned int,unsigned int> mpfsPerTPC;

    for(auto const mpf : *allHits){
      mpfsPerTPC[mpf.WireID().TPC]++;
    }
    
    for(auto const m:  mpfsPerTPC){
      if(m.second > fHitLimit){
        result = false;
        break;
      }
    }
  }
  else{
    // This is the simplest thing we can do, just cut on the total number of mpfs
    if (allHits->size() > fHitLimit){
      result = false;
    }
  }


  return result;
}

void mpf::MissingProductFilter::endJob() {}
 
DEFINE_ART_MODULE(mpf::MissingProductFilter)
