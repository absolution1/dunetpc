// AdcChannelDftPlotter.h

// David Adams
// April 2018
//
// Tool to plot DFT info for ADC channel data.
//
// Four plots are supported: graphs of DFT magnitude and phase and
// histograms of power and power/tick. All are shown as a function of frequency.
//
// The DFT magnitude and phase are drawn with graphs with one point for
// each relevant frequency, i.e. (Nsam + 2)/2 values.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   Variable - Variable to plot:
//                magnitude  - graph of DFT magnitude
//                phase      - graph of DFT phase
//                power      - histogram of power
//                power/tick - histogram of power/tick
//   SampleFreq - Sampling frequency in kHz used to define the x-axis.
//                If zero, frequency index is used instead.
//   YMax - Specifies the maximum value for the Y-axis in the plot:
//              0: Use automatic value
//            > 0: Use YMax
//            < 0: Use -YMax if there are no higher values, otherwise automatic
//   NbinX - # bins in the power histogram
//   HistName - Histogram name
//   HistTitle - Histogram title
//   PlotName - Plot name. If not blank, plot file is created.

#ifndef AdcChannelDftPlotter_H
#define AdcChannelDftPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AdcChannelStringTool;

class AdcChannelDftPlotter : AdcChannelTool {

public:

  using Index = unsigned int;
  using Name = std::string;

  AdcChannelDftPlotter(fhicl::ParameterSet const& ps);

  ~AdcChannelDftPlotter() override =default;

  DataMap view(const AdcChannelData& acd) const override;
  bool updateWithView() const override { return true; }

private:

  // Configuration data.
  int    m_LogLevel;
  Name   m_Variable;
  float  m_SampleFreq;
  float  m_YMax;
  Index  m_NBinX;
  Name   m_HistName;
  Name   m_HistTitle;
  Name   m_PlotName;

  // ADC string tools.
  const AdcChannelStringTool* m_adcNameBuilder;
  const AdcChannelStringTool* m_adcTitleBuilder;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, bool isTitle) const;

};

DEFINE_ART_CLASS_TOOL(AdcChannelDftPlotter)

#endif
