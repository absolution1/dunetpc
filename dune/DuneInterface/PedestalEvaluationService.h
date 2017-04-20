// PedestalEvaluationService.h

#ifndef PedestalEvaluationService_H
#define PedestalEvaluationService_H

// David Adams
// May 2016
//
// Interface for a service that uses a vector of ADC signals to estimate
// the pedestal for a channel.

#include <iostream>
#include "dune/DuneInterface/AdcChannelData.h"

class PedestalEvaluationService {

public:

  // Extract the channel, signals and flags. All are pointer types so the caller can use null to
  // indicate any of the fields are not of interest.
  //      sigs - Input signal vector
  //     *pped - Pedestal
  //     *prem - Pedestal RMS
  //  *ppederr - Uncertainty in the pedestal
  //  *premerr - Uncertainty in the pedestal RMS
  // The last four arguments are output values. Null values may be used to indicate the value
  // is not used and need not be calculated.
  virtual int evaluate(const AdcChannelData& data,
                       AdcSignal* pped =nullptr, AdcSignal* prms =nullptr,
                       AdcSignal* ppederr =nullptr, AdcSignal* prmserr =nullptr) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(PedestalEvaluationService, LEGACY)
#endif

#endif
