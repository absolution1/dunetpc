// StandardRawDigitPrepService.h
//
// David Adams
// July 2016
//
// Implementation of service that prepares TPC raw digits for reconstruction.
// The data is converted to float, pedestals are subtracted and optionally
// stuck bits are mitigated and coherent noise is removed.
//
// Configuration:
//   LogLevel - message logging level: 0=none, 1=initialization, 2+=every event
//   DoMitigation - Run mitigation (e.g. stuck bit removal) on extracted data.
//   DoEarlySignalFinding - Run signal building before noise removal.
//   DoNoiseRemoval - Run coherent noise suppression.
//   DoPedestalAdjustment - Do dynamic pedestal adjustment.
//   DoDeconvolution - Deconvolute the signal.
//   DoROI - Build ROIs.

#ifndef StandardRawDigitPrepService_H
#define StandardRawDigitPrepService_H

#include "dune/DuneInterface/RawDigitPrepService.h"

class RawDigitExtractService;
class AdcMitigationService;
class AdcSignalFindingService;
class AdcNoiseRemovalService;
class PedestalEvaluationService;
class AdcDeconvolutionService;
class AdcRoiBuildingService;

class StandardRawDigitPrepService : public RawDigitPrepService {

public:

  StandardRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int prepare(const std::vector<raw::RawDigit>& digs, AdcChannelDataMap& prepdigs) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int m_LogLevel;
  bool m_DoMitigation;
  bool m_DoEarlySignalFinding;
  bool m_DoNoiseRemoval;
  bool m_DoPedestalAdjustment;
  bool m_DoDeconvolution;
  bool m_DoROI;

  const RawDigitExtractService* m_pExtractSvc;
  const AdcMitigationService* m_pmitigateSvc;
  const AdcSignalFindingService* m_pAdcSignalFindingService;
  const AdcNoiseRemovalService* m_pNoiseRemoval;
  const PedestalEvaluationService* m_pPedestalEvaluation;
  const AdcDeconvolutionService* m_pDeconvolutionService;
  const AdcRoiBuildingService* m_pRoiBuildingService;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitPrepService, RawDigitPrepService, LEGACY)

#endif
