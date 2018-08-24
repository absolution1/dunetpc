// AdcKeepAllSignalFinder.h

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
//
// The output results holds the status and
//  nroi - # rois (always 1)

#ifndef AdcKeepAllSignalFinder_H
#define AdcKeepAllSignalFinder_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <string>
#include <vector>

class AdcKeepAllSignalFinder
: public AdcChannelTool {

public:

  AdcKeepAllSignalFinder(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Name = std::string;

  // Configuration data.
  int m_LogLevel;

};

DEFINE_ART_CLASS_TOOL(AdcKeepAllSignalFinder)

#endif
