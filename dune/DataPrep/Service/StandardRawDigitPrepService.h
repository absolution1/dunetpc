// StandardRawDigitPrepService.h
//
// David Adams
// July 2016
//
// Implementation of service that prepares TPC raw digits for reconstruction.
// The data is converted to float, pedestals are subtracted and optionally
// stuck bits are mitigated and coherent noise is removed.
//
// There is an option to construct intermediate states with the following names:
//   extracted - After pedestal subtraction
//   mitigated - After mitigation (e.g. stuck bit interpolation)
//   noiseRemoved - After noise removal
//
// Configuration:
//   LogLevel - message logging level: 0=none, 1=initialization, 2+=every event
//   SkipBad - Skip bad channels as reported by ChannelStatusService.
//   SkipNoisy - Skip noisy channels as reported by ChannelStatusService.
//   ChannelStatusOnline [false] - If true, status is retrieved with online channel number.
//                       The current (Aug 2016) DUNE convention is to use offline numbering.
//   DoMitigation - Run mitigation (e.g. stuck bit removal) on extracted data.
//   DoEarlySignalFinding - Run signal building before noise removal.
//   DoNoiseRemoval - Run coherent noise suppression.
//   DoPedestalAdjustment - Do dynamic pedestal adjustment.
//   DoDeconvolution - Deconvolute the signal.
//   DoROI - Build ROIs.
//   DoWires - Build wires.
//   DoIntermediateStates - Build intermediate states requested by caller.
//   DoDump [false] - If true, the info for one tick is displayed in the log.
//   WiresWithoutROIFlag [2] - How to handle a request for wire building without ROI building.
//                               0 - Silently ignore
//                               1 - Warn in ctor only.
//                               2 - Warn in ctor and every event.
//                               3 - Abort after message.
//   DumpChannel [0] - The channel that is dumped.
//   DumpTick [0] - The tick that is dumped.

#ifndef StandardRawDigitPrepService_H
#define StandardRawDigitPrepService_H

#include "dune/DuneInterface/RawDigitPrepService.h"

class ChannelMappingService;
namespace lariov {
class ChannelStatusProvider;
}
class RawDigitExtractService;
class AdcMitigationService;
class AdcSignalFindingService;
class AdcNoiseRemovalService;
class PedestalEvaluationService;
class AdcDeconvolutionService;
class AdcRoiBuildingService;
class AdcWireBuildingService;
class AdcChannelDataCopyService;

class AdcChannelTool;

class StandardRawDigitPrepService : public RawDigitPrepService {

public:

  using AdcChannelToolPtr = std::unique_ptr<AdcChannelTool>;

  StandardRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int prepare(AdcChannelDataMap& prepdigs,
              std::vector<recob::Wire>* pwires,
              WiredAdcChannelDataMap* pintStates) const override;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const override;

private:

  // Configuration parameters.
  int m_LogLevel;
  bool m_SkipBad;
  bool m_SkipNoisy;
  bool m_ChannelStatusOnline;
  bool m_DoMitigation;
  bool m_DoEarlySignalFinding;
  bool m_DoNoiseRemoval;
  bool m_DoPedestalAdjustment;
  bool m_DoDeconvolution;
  bool m_DoROI;
  bool m_DoWires;
  unsigned int m_WiresWithoutROIFlag;
  bool m_DoDump;
  bool m_DoIntermediateStates;
  unsigned int m_DumpChannel;
  unsigned int m_DumpTick;
  std::vector<std::string> m_DisplayTools;

  const ChannelMappingService* m_pChannelMappingService;
  const lariov::ChannelStatusProvider* m_pChannelStatusProvider;
  const RawDigitExtractService* m_pExtractSvc;
  const AdcMitigationService* m_pmitigateSvc;
  const AdcSignalFindingService* m_pAdcSignalFindingService;
  const AdcNoiseRemovalService* m_pNoiseRemoval;
  const PedestalEvaluationService* m_pPedestalEvaluation;
  const AdcDeconvolutionService* m_pDeconvolutionService;
  const AdcRoiBuildingService* m_pRoiBuildingService;
  const AdcWireBuildingService* m_pWireBuildingService;
  const AdcChannelDataCopyService* m_pAdcChannelDataCopyService;

  std::vector<AdcChannelToolPtr> m_DisplayToolPtrs;
  std::vector<std::string> m_FoundDisplayToolNames;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitPrepService, RawDigitPrepService, LEGACY)

#endif
