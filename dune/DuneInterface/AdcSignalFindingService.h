// AdcSignalFindingService.h

#ifndef AdcSignalFindingService_H
#define AdcSignalFindingService_H

// David Adams
// June 2016
//
// Interface for a service that finds signals in ADC channel data.
// of ADC channel signals.

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

class AdcSignalFindingService {

public:

  // Extract the channel, signals and flags. All are pointer types so the caller can use null to
  // indicate any of the fields are not of interest.
  virtual int find(AdcChannelData& data) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

DECLARE_ART_SERVICE_INTERFACE(AdcSignalFindingService, LEGACY)

#endif
