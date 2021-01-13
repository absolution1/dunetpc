// AdcSignalFindingService.h

#ifndef AdcSignalFindingService_H
#define AdcSignalFindingService_H

// David Adams
// June 2016
//
// Interface for a service that finds signals in ADC channel data.

#include <iostream>
#include "dune/DuneInterface/Data/AdcChannelData.h"

class AdcSignalFindingService {

public:

  virtual ~AdcSignalFindingService() =default;

  // Find signals and set data.signal correspondingly.
  virtual int find(AdcChannelData& data) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcSignalFindingService, LEGACY)
#endif

#endif
