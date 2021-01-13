// MultiAdcSuppressService.h
//
// David Adams
// December 2015
//
// Service interface for identifying zero-suppressed ticks for multiple
// ADC channels, e.g. for a ROP.

#ifndef MultiAdcSuppressService_H
#define MultiAdcSuppressService_H

#include <string>
#include <iostream>
#include <memory>
#include <string>
#include <iostream>
#include "dune/DuneInterface/AdcTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

namespace fhicl {
class ParameterSet;
}
namespace art {
class ActivityRegistry;
}

class MultiAdcSuppressService {

public:

  // Dtor.
  virtual ~MultiAdcSuppressService() =default;

  // Update the filters for a collection of ADC channels with selection.
  virtual int filter(const AdcCountSelectionVector& acss) =0;

  // Print the configuration.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="  ") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(MultiAdcSuppressService, LEGACY)
#endif

#endif
