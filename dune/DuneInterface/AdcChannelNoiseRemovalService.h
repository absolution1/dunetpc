// AdcChannelNoiseRemovalService.h

#ifndef AdcChannelNoiseRemovalService_H
#define AdcChannelNoiseRemovalService_H

// David Adams
// August 2016
//
// Interface for a service that removes noise from an ADC channel.
// Note that AdcNoiseRemovalService removes noise for multiple channels.

#include <iostream>
#include "dune/DuneInterface/AdcChannelData.h"

class AdcChannelNoiseRemovalService {

public:

  // Remove the noise for selected channels from data.samples and store results
  // in that same vector.
  // Return 0 for success.
  virtual int update(AdcChannelData& data) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcChannelNoiseRemovalService, LEGACY)
#endif

#endif
