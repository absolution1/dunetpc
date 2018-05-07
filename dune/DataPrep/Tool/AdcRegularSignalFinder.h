// AdcRegularSignalFinder.h

// David Adams
// December 2017
//
// Tool to find and flag regularly-spaced signal regions in ADC data.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   Period - Period for ROIs, i.e. one ROI is found every Period ticks
//   Length - # ticks in each ROI (0 means same as Period)
//
// If Period is zero, then the period and length are obtained from the
// data passed in each call to update (or view). The period is the
// length of the first ROI (roi.second + 1 - roi.first).
// The configured length is the length of the second ROI, if present,
// or otherwise is the same as the period.
//
// The first ROI starts at tick 0. We may later add a parameter "Offset" to
// start at a later tick. The last ROI ends at the last tick and may be
// shorter than the others.
//
// The output results holds:
//   int roiPeriod - period (actual)
//   int roiLength - length (actual)
//   int roiCount  - # ROIs found

#ifndef AdcRegularSignalFinder_H
#define AdcRegularSignalFinder_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <string>
#include <vector>

class HistogramManager;
class TH1;

class AdcRegularSignalFinder
: public AdcChannelTool {

public:

  AdcRegularSignalFinder(AdcIndex per, AdcIndex len, int lev);
  AdcRegularSignalFinder(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Name = std::string;

  // Configuration data.
  int m_LogLevel;
  AdcIndex m_Period;
  AdcIndex m_Length;

};

DEFINE_ART_CLASS_TOOL(AdcRegularSignalFinder)

#endif
