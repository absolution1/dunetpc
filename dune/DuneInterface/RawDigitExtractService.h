// RawDigitExtractService.h

#ifndef RawDigitExtractService_H
#define RawDigitExtractService_H

// David Adams
// May 2016
//
// Interface for a service that extracts the ADC signal vector from
// a larsoft raw digit. The latter holds the TPC samples for one channel.

#include "dune/DuneInterface/AdcTypes.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace raw {
class RawDigit;
}

class RawDigitExtractService {

public:

  // Extract the signals and flags. Both are pointer types so the caller can use null to
  // indicate one of the fields is not of interest.
  virtual int extract(const raw::RawDigit& rawin,
                      AdcSignalVector* psigs, AdcFlagVector* pflags =nullptr) =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

DECLARE_ART_SERVICE_INTERFACE(RawDigitExtractService, LEGACY)

#endif
