// AdcWireBuildingService.h

#ifndef AdcWireBuildingService_H
#define AdcWireBuildingService_H

// David Adams
// June 2016
//
// Interface for a service that builds a recob::Wire in ADC channel data.

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "lardataobj/RecoBase/Wire.h"

class AdcWireBuildingService {

public:

  typedef std::vector<recob::Wire> WireVector;

  // Build wire info.
  // If wires is not null, the Wire object is stored and managed there.
  // Otherwise the caller must manage the wire whose pointer is data.wire.
  virtual int build(AdcChannelData& data, WireVector* wires) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

DECLARE_ART_SERVICE_INTERFACE(AdcWireBuildingService, LEGACY)

#endif
