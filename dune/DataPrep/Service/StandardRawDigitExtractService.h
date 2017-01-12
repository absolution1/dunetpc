// StandardRawDigitExtractService
//
// David Adams
// May 2016
//
// Implementation of service that is first step in processing TPC data.
// It uncompressses the data, converts to float, subtracts pedestals and
// flags underflows, overflows and stuck bits.
//
// Flags are set assuming a 12-bit ADC where the low six bits can be stuck.
//
// Configuration:
//   LogLevel - message logging level: 0=none, 1=initialization, 2+=every event
//   PedestalOption
//     1: take from digit
//     2: take from service
//     3: evaluate with PedestalEvaluationService
//     otherwise no subtraction
//   FlagStuckOff - flag samples for which the low bits are all 0
//   FlagStuckOn  - flag samples for which the low bits are all 1


#ifndef StandardRawDigitExtractService_H
#define StandardRawDigitExtractService_H

#include "dune/DuneInterface/RawDigitExtractService.h"
#include "dune/DuneInterface/PedestalEvaluationService.h"

namespace lariov {
  class DetPedestalProvider;
}

class StandardRawDigitExtractService : public RawDigitExtractService {

public:

  StandardRawDigitExtractService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int extract(AdcChannelData& acd) const;

  int extract(const raw::RawDigit& rawin,
              AdcChannel* pchan =nullptr,
              AdcSignal* pped =nullptr,
              AdcCountVector* praw =nullptr,
              AdcSignalVector* psigs =nullptr,
              AdcFlagVector* pflags =nullptr) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int m_LogLevel;
  int m_PedestalOption;
  bool m_FlagStuckOff;
  bool m_FlagStuckOn;

  const lariov::DetPedestalProvider* m_pPedProv;
  PedestalEvaluationService* m_PedestalEvaluationService;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitExtractService, RawDigitExtractService, LEGACY)

#endif
