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
//   LogLevel  - Logging opt (0=none, 1=init only, 2=every call, ...)
//   XXXChannelRanges - If not empty, a plot is created for each named range.
//                      If empty, a plot is created for each channel.
//   XXXName        - Name for the multi-plot event file
//   XXXSummaryName - Name for the multi-plot summary file
//   XXXSizeX       - XSize in pixels of the multi-plot file.
//   XXXSizeY       - YSize in pixels of the multi-plot file.
//   XXXSplitX      - NX
//   XXXSplitY      - NY

#ifndef AdcMultiChannelPlotter_H
#define AdcMultiChannelPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AdcChannelStringTool;
class TPadManipulator;

class AdcMultiChannelPlotter : public AdcChannelTool {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using AcdVector = std::vector<const AdcChannelData*>;
  using ChannelRangeMap = std::map<Name, IndexRange>;

  AdcMultiChannelPlotter(const fhicl::ParameterSet& ps, Name prefix ="Plot");

  ~AdcMultiChannelPlotter() override;

  // Loop over channels, create plots and call viewMapChannel to fill them.
  // Return could includes data from those calls (not implemented) and does include
  //   multiChannelNChannel - # channels processed
  //   multiChannelNPlot    - # plots produced
  DataMap viewMap(const AdcChannelDataMap& acds) const override;

  // Loop over channels, create summary plots and call viewMapSummary to fill them.
  // Subclass must call this from its dtor.
  void viewSummary() const;

  // Subclass provides this method to view one channel in a map.
  //   acd - Data for the current channel
  //   ret - Result of processing previous channels. Append channel results here.
  //   man - Empty pad to be filled with the plot for this channel.
  //virtual int viewMapChannel(const AdcChannelData& acd, DataMap& ret, TPadManipulator& man) const =0;

  // Subclass provides this method to process the channels in one range from one call to view.
  //    crn - Name for this set of channels
  //   acds - Data for the channels
  //    man - Empty pad to be filled with the plot for this channel.
  virtual int viewMapChannels(Name crn, const AcdVector& acds, TPadManipulator& man) const =0;

  // Subclass provides this method to process the channels in one range for a summary plot,
  // i.e. combining data from preceding view calls.
  //    crn - Name for this set of channels
  //   acds - Data for the channels
  //    man - Empty pad to be filled with the plot for this channel.
  virtual int viewMapSummary(Name crn, TPadManipulator& man) const =0;

  // Provide read access to configuration.
  Index getLogLevel() const { return m_LogLevel; }
  const NameVector& getChannelRangeNames() const { return m_PlotChannelRanges; }
  bool haveChannelRanges() const { return getChannelRangeNames().size(); }
  Name  getPlotName() const { return m_PlotName; }
  Name  getPlotSummaryName() const { return m_PlotSummaryName; }
  Index getPlotSizeX() const { return m_PlotSizeX; }
  Index getPlotSizeY() const { return m_PlotSizeY; }
  Index getPlotSplitX() const { return m_PlotSplitX; }
  Index getPlotSplitY() const { return m_PlotSplitY; }

private:

  // Configuration data.
  Index m_LogLevel;
  NameVector m_PlotChannelRanges;
  Name  m_PlotName;
  Name  m_PlotSummaryName;
  Index m_PlotSizeX;
  Index m_PlotSizeY;
  Index m_PlotSplitX;
  Index m_PlotSplitY;

  // Derived from configuration.
  NameVector m_crns;
  ChannelRangeMap m_crmap;

  // ADC string tools.
  const AdcChannelStringTool* m_adcStringBuilder;

  // State data.
  //   run = run number from first call to view.
  class BaseState {
    IndexVector m_runs;
  public:
    void setRun(Index irun) { m_runs.clear(); m_runs.resize(1, irun); }
    bool hasRun() const { return m_runs.size(); }
    Index run() const { return hasRun() ? m_runs[0] : 0; }
  };
  BaseState m_baseState;

protected:

  BaseState& getBaseState() const { return const_cast<BaseState&>(m_baseState); }

};

#endif
