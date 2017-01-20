// AdcRoiBuildingService.h

#ifndef AdcRoiBuildingService_H
#define AdcRoiBuildingService_H

// David Adams
// June 2016
//
// Interface for a service that builds ROIs in ADC channel data.

#include <iostream>
#include "dune/DuneInterface/AdcChannelData.h"

class AdcRoiBuildingService {

public:

  // Build ROIs and set data.rois correspondingly.
  virtual int build(AdcChannelData& data) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcRoiBuildingService, LEGACY)
#endif

#endif
