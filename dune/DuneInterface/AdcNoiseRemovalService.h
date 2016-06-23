// AdcNoiseRemovalService.h

#ifndef AdcNoiseRemovalService_H
#define AdcNoiseRemovalService_H

// David Adams
// May 2016
//
// Interface for a service that removes coherent noise from a collection
// of ADC channel signals.

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

class AdcNoiseRemovalService {

public:

  // Extract the channel, signals and flags. All are pointer types so the caller can use null to
  // indicate any of the fields are not of interest.
  virtual int update(AdcChannelDataMap& datamap) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

DECLARE_ART_SERVICE_INTERFACE(AdcNoiseRemovalService, LEGACY)

#endif
