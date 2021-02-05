// AdcMultiChannelPlotter.h
//
// AdcMultiChannelPlotter is an intermediate base for 
// AdcChannelTool that implements viewMap by looping over
// channel ranges and calling the subclass method
// viewMapChannels with the data for the channels in that group.
//
// Each call includes a pad that the receiver is expected to fill.
// Plots with names derived from the template XXXName are
// produced from these pads. There is one pad for each channel group.
// 
// In addition, this class provides the method viewSummary(ilev) which
// produces plots with names derived from XXXSummaryNames[ilev].
// The pads are filled by calls to the subclass method viewMapSummary
// for each channel range with again one pad for each channel group.
// Calls to viewSummary are made by the subclass. E.g it might create
// event-level plots from endEvent or job level plots from its dtor.
//
// Both channel ranges and range groups may be specified and the subclass
// calls are made for the ranges in each.
// The subclass provides the method viewMapChannels that is called with the
// channel data for
//   each range in the range list and
//   each range in each of the groups in the group list.
// Each range in the range list is drawn on a separate pad.
// If XXXOverlayGroups is true, there is one drawing pad per channel group.
// Otherwise, there is one for each range in each group.
// If no ranges or groups are specified, a separate plot is made for each channel.
//
// Each plot file shows up to NX*NY channels in a NY x NX array.
// The overall plot size is YSize pixels by XSize pixels or Root
// default (500 x 700?) if either is zero.
//
// Configuration (XXX is prefix supplied in ctor):
//   LogLevel  - Logging opt (0=none, 1=init only, 2=every call, ...)
//   XXXChannelRanges - Channel ranges.
//   XXXChannelGroups - Channel range groups.
//   XXXOverlayGroups - Flag indicating if ranges in a group appear on the same pad.
//   XXXDataView      - Data view to plot ("" = top)
//   XXXName          - Name for the multi-plot file created for each call to view.
//   XXXSummaryNames  - Name for the multi-plot summary file
//   XXXSizeX         - XSize in pixels of the multi-plot file.
//   XXXSizeY         - YSize in pixels of the multi-plot file.
//   XXXSplitX        - NX
//   XXXSplitY        - NY
//
// To obtain summary plots, viewMap must called (view is no sufficient).

#ifndef AdcMultiChannelPlotter_H
#define AdcMultiChannelPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Data/IndexRangeGroup.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <set>

class AdcChannelStringTool;
class TPadManipulator;

class AdcMultiChannelPlotter : public TpcDataTool {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;
  using AcdVector = std::vector<const AdcChannelData*>;
  using ChannelRangeMap = std::map<Name, IndexRange>;
  using ChannelGroupMap = std::map<Name, IndexRangeGroup>;

  // Base class for the description of one drawing pad.
  class Pad {
  public:
    Name cgname;            // Channel group name
    NameVector crnames;     // Ordered vector of channel range names
    ChannelRangeMap crmap;  // Channel range associated with each name
    // Create a pad with no channel ranges.
    Pad(Name a_cgname) : cgname(a_cgname) { }
  };

  using PadVector = std::vector<Pad>;

  AdcMultiChannelPlotter(const fhicl::ParameterSet& ps, Name prefix ="Plot");

  ~AdcMultiChannelPlotter() override;

  // Loop over channels, create plots and call viewMapChannel to fill them.
  // Return could includes data from those calls (not implemented) and does include
  //   multiChannelNChannel - # channels processed
  //   multiChannelNPlot    - # plots produced
  DataMap viewMap(const AdcChannelDataMap& acds) const override;

  // Loop over channels, create summary plots at level ilev and call viewMapSummary to fill them.
  // Subclass calls this at appropriate levels.
  void viewSummary(Index ilev) const;

  // Subclass provides this method to process the channels in one range from one call to view.
  // If a data view has been specified, then there be zero, one or one data object for each channel.
  //    crn - Name for this set of channels
  //   acds - Data for the channels
  //    man - Pad to be filled with the plot for this channel.
  //          Subclass can use man.haveHistOrGraph() to see if it has previously filled the pad.
  //    ncr - # channel ranges in this plot
  //    icr - index of the channel ranges in this plot
  virtual int viewMapChannels(Name crn, const AcdVector& acds, TPadManipulator& man, Index ncr, Index icr) const =0;

  // Subclass provides this method to process the channels in one range for a summary plot,
  // i.e. combining data from preceding view calls.
  //    cgn - Name for the group holding the channel range
  //    crn - Name for the channel range
  //   acds - Data for the channels
  //    man - Pad to be filled with the plot for this channel.
  //          Subclass can use man.haveHistOrGraph() or icrn to see if it has previously filled the pad.
  //    ncr - # channel ranges in this plot
  //    icr - index of the channel ranges in this plot
  virtual int viewMapSummary(Index ilev, Name cgn, Name crn, TPadManipulator& man, Index ncr, Index icr) const =0;

  // Provide read access to configuration.
  Index getLogLevel() const { return m_LogLevel; }
  Name getDataView() const { return m_PlotDataView; }
  const NameVector& getChannelRangeNames() const { return m_PlotChannelRanges; }
  const NameVector& getChannelGroupNames() const { return m_PlotChannelGroups; }
  bool haveChannelRanges() const { return getChannelRangeNames().size(); }
  bool haveChannelGroups() const { return getChannelGroupNames().size(); }
  bool overlayGroups() const { return m_PlotOverlayGroups; }
  Name  getPlotName() const { return m_PlotName; }
  Name  getPlotSummaryName(Index ilev) const {
    if ( ilev >= m_PlotSummaryNames.size() ) return "";
    return m_PlotSummaryNames[ilev];
  }
  Index getPlotSizeX() const { return m_PlotSizeX; }
  Index getPlotSizeY() const { return m_PlotSizeY; }
  Index getPlotSplitX() const { return m_PlotSplitX; }
  Index getPlotSplitY() const { return m_PlotSplitY; }

  // Return the channel groups.
  const IndexRangeGroup& getChannelGroup(Name cgn) const;

private:

  // Configuration data.
  Index m_LogLevel;
  NameVector m_PlotChannelRanges;
  NameVector m_PlotChannelGroups;
  Index m_PlotOverlayGroups;
  Name  m_PlotDataView;
  Name  m_PlotName;
  NameVector  m_PlotSummaryNames;
  Index m_PlotSizeX;
  Index m_PlotSizeY;
  Index m_PlotSplitX;
  Index m_PlotSplitY;

  // Derived from configuration.
  NameVector m_crns;
  ChannelRangeMap m_crmap;
  NameVector m_cgns;
  ChannelGroupMap m_cgmap;
  PadVector m_pads;

  // ADC string tools.
  const AdcChannelStringTool* m_adcStringBuilder;

  // State data.
  //   run = run number from first call to view.
  class BaseState {
    IndexVector m_runs;
  public:
    Index event =0;
    void setRun(Index irun) { m_runs.clear(); m_runs.resize(1, irun); }
    bool hasRun() const { return m_runs.size(); }
    Index run() const { return hasRun() ? m_runs[0] : 0; }
    IndexSet channels;    // Channels for no ranges or groups
  };
  BaseState m_baseState;

protected:

  BaseState& getBaseState() const { return const_cast<BaseState&>(m_baseState); }

};

#endif
