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
//   DigitReadTool name for the tool that reads digits into ADC channel data
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

class AdcChannelTool;
namespace lariov {
  class DetPedestalProvider;
}

class StandardRawDigitExtractService : public RawDigitExtractService {

public:

  StandardRawDigitExtractService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int extract(AdcChannelData& acd) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  using AcdModifierPtr = std::unique_ptr<const AdcChannelTool>;

  // Configuration parameters.
  int         m_LogLevel;
  std::string m_ROIBuilderTool;
  std::string m_DigitReadTool;
  int         m_PedestalOption;
  bool        m_FlagStuckOff;
  bool        m_FlagStuckOn;

  AcdModifierPtr m_pDigitReadTool;
  AcdModifierPtr m_pROIBuilderTool;

  const lariov::DetPedestalProvider* m_pPedProv;
  PedestalEvaluationService* m_PedestalEvaluationService;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitExtractService, RawDigitExtractService, LEGACY)

#endif
