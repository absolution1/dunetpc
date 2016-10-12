// AdcDeconvolutionService.h

#ifndef AdcDeconvolutionService_H
#define AdcDeconvolutionService_H

// David Adams
// May 2016
//
// Interface for a service that modifies the signals for an ADC channel.

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

class AdcDeconvolutionService {

public:

  // Extract the channel, signals and flags. All are pointer types so the caller can use null to
  // indicate any of the fields are not of interest.
  virtual int update(AdcChannelData& data) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

DECLARE_ART_SERVICE_INTERFACE(AdcDeconvolutionService, LEGACY)

#endif
