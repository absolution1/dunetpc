// AdcSuppressService.h
//
// David Adams
// December 2015
//
// Service interface for identifying zero-suppressed ticks in an ADC vector.

#ifndef AdcSuppressService_H
#define AdcSuppressService_H

#include <string>
#include <iostream>
#include <memory>
#include <string>
#include <iostream>
#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcCountSelection.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

namespace fhicl {
class ParameterSet;
}
namespace art {
class ActivityRegistry;
}

class AdcSuppressService {

public:

  typedef unsigned int Channel;

  // Dtor.
  virtual ~AdcSuppressService() =default;

  // Suppress a vector of signals.
  //     sigs: Input ADC vector.
  //      ped: Input pedestal.
  //     keep: I/O boolean vector indicating which signals are retained.
  virtual int filter(const AdcCountVector& sigs,
                     Channel chan,
                     AdcPedestal ped,
                     AdcFilterVector& keep) const =0;

  // Print the configuration.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="  ") const =0;

  // Alternate interface.
  int filter(AdcCountSelection& acs) {
    return filter(acs.counts, acs.channel, acs.pedestal, acs.filter);
  }

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcSuppressService, LEGACY)
#endif

#endif
