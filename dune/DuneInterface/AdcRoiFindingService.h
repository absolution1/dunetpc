// AdcRoiFindingService.h

#ifndef AdcRoiFindingService_H
#define AdcRoiFindingService_H

// David Adams
// June 2016
//
// Interface for a service that finds signals in ADC channel data.

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

class AdcRoiFindingService {

public:

  // Find ROIs and set data.rois correspondingly.
  virtual int find(AdcChannelData& data) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

DECLARE_ART_SERVICE_INTERFACE(AdcRoiFindingService, LEGACY)

#endif
