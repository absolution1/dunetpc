// AdcMultiChannelPlotter.h
//
// AdcChannelTool intermediate base that constructs plots for
// multiple channels. It implements the viewMap methods as a loop
// over calls to view and uses the result from each of those to
// add to a pad with combined results.
//
// Each plot file shows up to NX*NY channels in a NY x NX array.
// The overall plot size is YSize pixels by XSize pixels or Root
// defualt (500 x 700?) if either is zero.
//
// Configuration (XXX is prefix supplied in ctor):
//   LogLevel  - Logginf opiot (0=none, 1=init only, 2=every call, ...)
//   XXXName   - Name for the multi-plot file
//   XXXSizeX  - XSize in pixels of the multi-plot file.
//   XXXSizeY  - YSize in pixels of the multi-plot file.
//   XXXSplitX - NX
//   XXXSplitY - NY

#ifndef AdcMultiChannelPlotter_H
#define AdcMultiChannelPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AdcChannelStringTool;
class TPadManipulator;

class AdcMultiChannelPlotter : public AdcChannelTool {

public:

  using Name = std::string;
  using Index = unsigned int;

  AdcMultiChannelPlotter(const fhicl::ParameterSet& ps, Name prefix ="Plot");

  ~AdcMultiChannelPlotter() override =default;

  // Loop over channels, create plots and call viewMapChannel for each.
  // Return includes data from those calls plus 
  //   multiChannelNChannel - # channels processed
  //   multiChannelNPlot    - # plots produced
  DataMap viewMap(const AdcChannelDataMap& acds) const override;

  // Subclass provides this method to view one channel in a map.
  //   acd - Data for the current channel
  //   ret - Result of processing previous channels. Append channel results here.
  //   man - Empty pad to be filled with the plot for this channel.
  virtual int viewMapChannel(const AdcChannelData& acd, DataMap& ret, TPadManipulator& man) const =0;

  // Provide read access to configuration.
  Index getLogLevel() const { return m_LogLevel; }
  Name  getPlotName() const { return m_PlotName; }
  Index getPlotSizeX() const { return m_PlotSizeX; }
  Index getPlotSizeY() const { return m_PlotSizeY; }
  Index getPlotSplitX() const { return m_PlotSplitX; }
  Index getPlotSplitY() const { return m_PlotSplitY; }

private:

  // Configuration data.
  Index m_LogLevel;
  Name  m_PlotName;
  Index m_PlotSizeX;
  Index m_PlotSizeY;
  Index m_PlotSplitX;
  Index m_PlotSplitY;

  // ADC string tools.
  const AdcChannelStringTool* m_adcStringBuilder;

};

#endif
