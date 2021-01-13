// AdcDeconvolutionService.h

#ifndef AdcDeconvolutionService_H
#define AdcDeconvolutionService_H

// David Adams
// May 2016
//
// Interface for a service that modifies the signals for an ADC channel.

#include <iostream>
#include "dune/DuneInterface/Data/AdcChannelData.h"
namespace detinfo {
  class DetectorClocksData;
}

class AdcDeconvolutionService {

public:

  virtual ~AdcDeconvolutionService() = default;

  // Extract the channel, signals and flags. All are pointer types so the caller can use null to
  // indicate any of the fields are not of interest.
  virtual int update(detinfo::DetectorClocksData const& clockData,
                     AdcChannelData& data) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcDeconvolutionService, LEGACY)
#endif

#endif
