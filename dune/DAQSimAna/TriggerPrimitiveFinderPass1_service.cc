////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveFinderPass1
// Plugin Type: service (art v2_10_03)
// File:        TriggerPrimitiveFinderPass1_service.cc
//
// Generated at Tue Jun  5 07:51:38 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "dune/DAQSimAna/TriggerPrimitiveFinderService.h"

class TriggerPrimitiveFinderPass1;


class TriggerPrimitiveFinderPass1 : public TriggerPrimitiveFinderService {
public:
  explicit TriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p, art::ActivityRegistry & areg);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

    virtual std::vector<TriggerPrimitiveFinderService::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples);
    

private:

  // Declare member data here.

};


TriggerPrimitiveFinderPass1::TriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p, art::ActivityRegistry & areg)
// :
// Initialize member data here.
{
}

std::vector<TriggerPrimitiveFinderService::Hit>
TriggerPrimitiveFinderPass1::findHits(const std::vector<unsigned int>& channel_numbers, 
                                      const std::vector<std::vector<short>>& collection_samples)
{
    return std::vector<TriggerPrimitiveFinderService::Hit>();
}

DECLARE_ART_SERVICE(TriggerPrimitiveFinderPass1, LEGACY)
DEFINE_ART_SERVICE(TriggerPrimitiveFinderPass1)
