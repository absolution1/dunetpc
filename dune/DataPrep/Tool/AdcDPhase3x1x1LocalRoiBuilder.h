// AdcDPhase3x1x1LocalRoiBuilder.h

// Christoph Alt
// October 2017
//
// Tool to find signals and build ROI's with respect to a local baseline (instead of a global one which is typically 0)
// Introduced because of slow baseline fluctuations in 3x1x1 dual phase data
//


#ifndef AdcDPhase3x1x1LocalRoiBuilder_H
#define AdcDPhase3x1x1LocalRoiBuilder_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneInterface/AdcTypes.h"

class AdcDPhase3x1x1LocalRoiBuilder : AdcChannelTool {

public:

  AdcDPhase3x1x1LocalRoiBuilder(fhicl::ParameterSet const& ps);

  ~AdcDPhase3x1x1LocalRoiBuilder() override;

  // Build ROIs
  DataMap update(AdcChannelData& acd) const override;

private:

  // Parameters.
  int m_LogLevel;
  AdcIndex m_BinsToAverageForPedestal;
  AdcIndex m_BinsToSkip;
  bool m_UseStandardDeviation;
  AdcIndex m_NConsecBinsAboveThreshold1;
  AdcSignal m_NSigmaStart1;
  AdcSignal m_NSigmaEnd1;
  AdcIndex m_NConsecBinsAboveThreshold2;
  AdcSignal m_NSigmaStart2;
  AdcSignal m_NSigmaMax;
  AdcIndex m_PadLow;
  AdcIndex m_PadHigh;
};

DEFINE_ART_CLASS_TOOL(AdcDPhase3x1x1LocalRoiBuilder)

#endif
