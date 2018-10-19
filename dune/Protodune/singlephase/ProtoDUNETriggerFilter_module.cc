// A very simple module to filter out the events with beam triggers
// leigh.howard.whitehead@cern.ch

#include <iostream>
#include <utility>
#include <set>

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 

#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"

namespace filt{

  class ProtoDUNETriggerFilter : public art::EDFilter {
    public:
      explicit ProtoDUNETriggerFilter(fhicl::ParameterSet const & pset);
      virtual ~ProtoDUNETriggerFilter() {};
      virtual bool filter(art::Event& e);

    private:

  };

  ProtoDUNETriggerFilter::ProtoDUNETriggerFilter::ProtoDUNETriggerFilter(fhicl::ParameterSet const & pset){}

  bool ProtoDUNETriggerFilter::filter(art::Event & evt){

    // The ProtoDUNE data utility tells us if we have a beam trigger
    protoana::ProtoDUNEDataUtils dataUtil;
    return dataUtil.IsBeamTrigger(evt);
  }

  DEFINE_ART_MODULE(ProtoDUNETriggerFilter)

}
