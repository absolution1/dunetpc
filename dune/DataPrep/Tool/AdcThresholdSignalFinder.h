// AdcThresholdSignalFinder.h

// David Adams
// October 2017
//
// Tool to find and flag signal regions in ADC data.
//
// For each tick with sample above Threshold and/or below -Threshold, the region
// [bin-BinsBefore, bin+binsAfter] is flagged as signal.
// An ROI is created for each range of contiguous signals.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   Threshold  - threshold for signal finding
//   BinsBefore - lower limit for signal range
//   BinsAfter  - upper limit for signal range
//   FlagPositive - Flag signals above Threshold
//   FlagNegative - Flag signals below Threshold
//
// The output results holds:
//   int nThresholdBins - # bins above threshold

#ifndef AdcThresholdSignalFinder_H
#define AdcThresholdSignalFinder_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <string>
#include <vector>

class AdcThresholdSignalFinder
: public AdcChannelTool {

public:

  AdcThresholdSignalFinder(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Name = std::string;

  // Configuration data.
  int m_LogLevel;
  float m_Threshold;
  unsigned int m_BinsBefore;
  unsigned int m_BinsAfter;
  bool m_FlagPositive;
  bool m_FlagNegative;

};

DEFINE_ART_CLASS_TOOL(AdcThresholdSignalFinder)

#endif
