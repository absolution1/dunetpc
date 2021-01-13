// AdcDistortionService.h

// David Adams
// December 2015
//
// Interface for a service that distorts the counts in an ADC vector.

#ifndef AdcDistortionService_H
#define AdcDistortionService_H

#include <vector>
#include <iostream>
#include "dune/DuneInterface/AdcTypes.h"

class AdcDistortionService {

public:

  typedef unsigned int Channel;

  virtual ~AdcDistortionService() =default;

  // Modify an input ADC vector.
  virtual int modify(Channel chan, AdcCountVector& adcvec) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcDistortionService, LEGACY)
#endif

#endif

