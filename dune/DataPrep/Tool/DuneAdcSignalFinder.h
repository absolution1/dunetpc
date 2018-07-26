// DuneAdcSignalFinder.h

// David Adams
// May 2018
//
// Tool to find and flag signal regions in ADC data using the
// algorithm developed for the 35t prototype TPC. An unpadded
// ROI starts for any signal above NSigmaStart sigma above the noise
// and continues until the level falls below NSigmaEnd sigma.
// NSigmaEnd was fixed to 1.0 in the original code.
//
// Configuration:
//   LogLevel    - usual log level
//   NoiseSigma  - Assumed noise level (sample units).
//                 If <=0, AdcChannelData::sampleNoise is used.
//   NSigmaStart - Level in sigma at which an unpadded signal starts.
//   NSigmaEnd   - Level in sigma at which an unpadded signal ends.
//   TicksBefore - Number of ticks to retain before unpadded region.
//   TicksAfter  - Number of ticks to retain after unpadded region.
//
// The output results holds:
//   nroi - # ROIs found

#ifndef DuneAdcSignalFinder_H
#define DuneAdcSignalFinder_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <string>
#include <vector>

class DuneAdcSignalFinder
: public AdcChannelTool {

public:

  DuneAdcSignalFinder(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Name = std::string;

  // Configuration data.
  int m_LogLevel;
  AdcSignal m_NoiseSigma;
  AdcSignal m_NSigmaStart;
  AdcSignal m_NSigmaEnd;
  AdcIndex m_TicksBefore;
  AdcIndex m_TicksAfter;

};

DEFINE_ART_CLASS_TOOL(DuneAdcSignalFinder)

#endif
