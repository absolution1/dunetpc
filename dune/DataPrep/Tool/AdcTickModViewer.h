// AdcTickModViewer.h

// David Adams
// June 2018
//
// Tool to view  and fit ADC tickmods.
//
// A tickmod is the tick number modulus a specified period, e.g. that
// from a pulser run.
//
// If FitRmsMin < FitRmsMax, the the RMS is constrained to the range
// (FitRmsMin, FitRmsMax) in the fit.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   TickModPeriod - Period for tickmods [ticks], i.e. # ticks when wf repeats
//   TimeOffsetTool - Name of the tool used to get the tick offset.
//                    If blank, no offset is applied.
//                    Otherwise integral offset only is used. Unit must be tick.
//   FitRmsMin: Lower limit for RMS fit range.
//   FitRmsMax: Upper limit for RMS fit range.
//   HistName:  Name for the histogram.
//   HistTitle: Title for the histogram.
//   HistChannelCount: # channels shown in histogram
//   PlotFileName: If nonblank, histograms are displayed in this file.
//   RootFileName: If nonblank, histogram is copied to this file.
//   TreeFileName: If nonblank, a tickmod tree is created in this file.
//   PlotChannels: If not empty, only the listed channels are plotted.
//   PlotSizeX, PlotSizeY: Size in pixels of the plot file.
//                         Root default (700x500?) is used if either is zero.
//   PlotShowFit: Flag indicating how fit should be displayed.
//                  >= 1: Show final fit
//                  >= 2: Show starting fit function
//   PlotSplitX: If this is nonzero, plots are created in updateMap (not update)
//               and the drawing canvas is split NY x NX where NX = PlotSplitX.
//               If PlotSplitX == 0, one canvas/plot is created in update.
//   PlotSplitY: If PlotSplitY > 0, then the above NY = PlotSplitY. Otherwise
//               NY = PlotSplitX. No effect if PlotSplitX == 0.
//               and up to that many plots are shown on the screen.
//   PlotWhich: Bit pattern indicating which plots to draw
//                1 - All tickmods starting with 0
//                2 - NX*NY tickmods around the minimum.
//                4 - NX*NY tickmods around the maximum.
//   PlotFrequency: 0 - Only make plots at end of job (in dtor).
//                  1 - Make plots every event
//
// Tools:
//   adcStringBuilder is used to make the following
// substitutions in the names and title:
//      %RUN%     --> run number
//      %SUBRUN%  --> subrun number
//      %EVENT%   --> event number
//      %CHAN%    --> channel number
//      %TICKMOD% --> tickmod (or Min or Max for those plots)
//
// The single-channel methods return a data map with the following:
//   tmCount:        - The tickmod period, e.g. 497
//   tmPlotCount:    - The number of tickmod plots
//   tmHists:        - The processed and narrowed tickmod histograms for this channel
//   tmWideHists:    - The raw tickmod histograms for this channel

#ifndef AdcTickModViewer_H
#define AdcTickModViewer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DataPrep/Utility/TickModTreeData.h"
#include <string>
#include <vector>
#include <memory>

class AdcChannelStringTool;
class TimeOffsetTool;
class TH1;
class TFile;
class TTree;
class TPadManipulator;

class AdcTickModViewer
: public AdcChannelTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using Name = std::string;
  using HistPtr = std::shared_ptr<TH1>;
  using HistVector = std::vector<HistPtr>;
  using HistVectorMap = std::map<Index, HistVector>;

  AdcTickModViewer(fhicl::ParameterSet const& ps);

  ~AdcTickModViewer();

  DataMap view(const AdcChannelData& acd) const override;
  bool updateWithView() const override { return true; }

private:

  // Configuration data.
  int m_LogLevel;
  int m_TickModPeriod;
  Name m_TimeOffsetTool;
  float m_FitRmsMin;
  float m_FitRmsMax;
  Name m_HistName;
  Name m_HistTitle;
  Index m_HistChannelCount;
  Name m_PlotFileName;
  Name m_RootFileName;
  Name m_TreeFileName;
  IndexVector m_PlotChannels;
  Index m_PlotSizeX;
  Index m_PlotSizeY;
  Index m_PlotShowFit;
  Index m_PlotSplitX;
  Index m_PlotSplitY;
  Index m_PlotWhich;
  Index m_PlotFrequency;

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

  // Tick offset tool.
  const TimeOffsetTool* m_tickOffsetTool;

  // This subclass carries the state for this tool, i.e. data that can change
  // after initialization.
  class State {
  public:
    // Tickmod Histograms mapped by channel.
    //    Full hists have 4096 channels.
    //    Proc hists have the region from the sticky code anlysis.
    HistVectorMap ChannelTickModFullHists;
    HistVectorMap ChannelTickModProcHists;
    // Run number for plot names.
    int run = -1;
    // Tree.
    TickModTreeData treedata;
    TFile* pfile =nullptr;
    TTree* tickmodTree =nullptr;
  };

  // Return the state.
  using StatePtr = std::shared_ptr<State>;
  StatePtr m_state;
  State& state() const { return *m_state; }

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, Index itkm) const;

  // Fill the histogram for a tickmod.
  int fillChannelTickMod(const AdcChannelData& acd, Index itkm0, Index itkm) const;

  // Process the accumulated histogram for a channel.
  // StickyCodeMetrics is used to create limited-range histogram and evaluate
  // sticky code metrics for that region.
  int processAccumulatedChannel(Index icha, Index& nplot) const;

  // Process the accumulated histogram for all channels.
  int processAccumulation(Index& nplot) const;

  // Create plot files for the tickmod histograms for a channel.
  //   icha - channel number
  //   nplot - # plot files created
  int makeTickModPlots(Index icha, Index& nplot) const;

};

DEFINE_ART_CLASS_TOOL(AdcTickModViewer)

#endif
