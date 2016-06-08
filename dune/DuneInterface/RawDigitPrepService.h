// RawDigitPrepService.h

#ifndef RawDigitPrepService_H
#define RawDigitPrepService_H

// David Adams
// May 2016
//
// Interface for a service that extracts the ADC signal vector from
// a larsoft raw digit. The latter holds the TPC samples for one channel.

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace raw {
class RawDigit;
}

class RawDigitPrepService {

public:

  // Prepare the data for a vector of TPC digits. Each digit holds the time samples for one channel.
  // The preparation includes extraction (float conversion and pedstal subtraction) and optionally
  // actions such as stuck bit mitigation and coherent noise removal.
  virtual int prepare(const std::vector<raw::RawDigit>& digs, AdcChannelDataMap& prepdigs) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

DECLARE_ART_SERVICE_INTERFACE(RawDigitPrepService, LEGACY)

#endif
