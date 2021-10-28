// AdcNoiseSignalFinder.h

// David Adams
// July 2019
//
// Tool to find and flag signal regions in ADC data using a threshold
// based on the noise level.
//
// For each tick with sample above threshold and/or below -threshold, the region
// [bin-BinsBefore, bin+binsAfter] is flagged as signal.
// An ROI is created for each range of contiguous signal samples.
//
// For
//   sigfrac = fraction of samples flagged as signal
//   noise = sample RMS for non-signal samples
//   tr = threshold/noise
//   trmin = ThresholdRatio - ThresholdRatioTol and
//   trmax = ThresholdRatio + ThresholdRatioTol:
// If sigfrac > SigFracMax or tr < trmin, the threshold is increased
// Else if tr > trmax and threshold > thresholdMin, the threshold is decreased
//
// The signal finding is repeated until there is no change in threshold or
// the loop count reaches MaxLoop.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 3=each event, >4=more
//              1 shows warnings for good channels, 2 for all
//   SigFracMax - maximum allowed value for sigfrac
//   ThresholdMin  - minimum and starting threshold
//   ThresholdRatio - Target threshold/noise
//   ThresholdRatioTol - Tolerance for threshold/noise
//   MaxLoop - Maximum # threshold iterations
//   BinsBefore - lower limit for signal range
//   BinsAfter  - upper limit for signal range
//   FlagPositive - Flag signals above Threshold
//   FlagNegative - Flag signals below Threshold
//
// The following are returnd and added to the ADC channel metadata:
//   float nsfSigFrac - fraction of samples above threshold
//   float nsfNoise - sample RMS outside of signal range
//   float nsfThreshold - final threshold
//   int nsfLoopCount - # threshold iterations
//   int nsfRoiCount - # ROIs

#ifndef AdcNoiseSignalFinder_H
#define AdcNoiseSignalFinder_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <string>
#include <vector>

class AdcNoiseSignalFinder
: public TpcDataTool {

public:

  AdcNoiseSignalFinder(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Name = std::string;

  // Configuration data.
  int m_LogLevel;
  float m_SigFracMax;
  float m_ThresholdMin;
  float m_ThresholdRatio;
  float m_ThresholdRatioTol;
  unsigned int m_MaxLoop;
  unsigned int m_BinsBefore;
  unsigned int m_BinsAfter;
  bool m_FlagPositive;
  bool m_FlagNegative;

};


#endif
