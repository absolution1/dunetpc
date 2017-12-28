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
// The output results holds:
//   roiPeriod - Period
//   roiLength - period
//   roiCount - # ROIs found

#ifndef AdcRegularSignalFinder_H
#define AdcRegularSignalFinder_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelDataModifier.h"
#include <string>
#include <vector>

class HistogramManager;
class TH1;

class AdcRegularSignalFinder
: public AdcChannelDataModifier {

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
