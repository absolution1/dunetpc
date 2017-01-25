// RawDigitExtractService.h

#ifndef RawDigitExtractService_H
#define RawDigitExtractService_H

// David Adams
// May 2016
//
// Interface for a service that extracts the ADC signal vector from
// a larsoft raw digit. The latter holds the TPC samples for one channel.

#include <iostream>
#include "dune/DuneInterface/AdcChannelData.h"

namespace raw {
class RawDigit;
}

class RawDigitExtractService {

public:

  // Extract the data for a channel.
  // The digit is assigned to the chanel data before the call.
  virtual int extract(AdcChannelData& acd) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(RawDigitExtractService, LEGACY)
#endif

#endif
