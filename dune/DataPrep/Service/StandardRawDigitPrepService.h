// StandardRawDigitPrepService.h
//
// David Adams
// May 2016
//
// Implementation of service that prepares TPC raw digits for reconstruction.
// The data is converted to float, pedestals are subtracted and optionally
// stuck bits are mitigated and coherent noise is removed.
//
// Configuration:
//   LogLevel - message logging level: 0=none, 1=initialization, 2+=every event
//   DoMitigation - Run mitigation (e.g. stuck bit removal) on extracted data.


#ifndef StandardRawDigitPrepService_H
#define StandardRawDigitPrepService_H

#include "dune/DuneInterface/RawDigitPrepService.h"

class RawDigitExtractService;
class AdcMitigationService;

class StandardRawDigitPrepService : public RawDigitPrepService {

public:

  StandardRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int prepare(const std::vector<raw::RawDigit>& digs, AdcChannelDataMap& prepdigs) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int m_LogLevel;
  bool m_DoMitigation;

  const RawDigitExtractService* m_pExtractSvc;
  const AdcMitigationService* m_pmitigateSvc;


};

DECLARE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitPrepService, RawDigitPrepService, LEGACY)

#endif
