// AdcCompressService.h
//
// David Adams
// December 2015
//
// Service interface for applying suppression and compressing an ADC vector.

#ifndef AdcCompressService_H
#define AdcCompressService_H

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

class AdcCompressService {

public:

  // Dtor.
  virtual ~AdcCompressService() =default;

  // Compress a vector of signals.
  //     sigs: I/O ADC vector to be compressed
  //     keep: Input boolean vector indicating which signals are retained.
  //   offset: Input pedestal. If used, the same value must be supplied for uncompression.
  //     comp: Output value indicating the compression strategy. Needed for uncompression.
  virtual int compress(AdcCountVector& sigs,
                       const AdcFilterVector& keep,
                       AdcCount offset,
                       raw::Compress_t& comp) const =0;

  // Print the configuration.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="  ") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcCompressService, LEGACY)
#endif

#endif
