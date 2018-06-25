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
//
// Tools:
//   adcStringBuilder is used to make the following
// substitutions in the names and title:
//      %RUN%     --> run number
//      %SUBRUN%  --> subrun number
//      %EVENT%   --> event number
//      %CHAN%    --> channel number
//
// The single-channel methods return a data map with the following:
//   tmHists:           - The tickmod histograms for this channel

#ifndef AdcTickModViewer_H
#define AdcTickModViewer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <string>
#include <vector>
#include <memory>

class AdcChannelStringTool;
class TimeOffsetTool;
class TH1;
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
  IndexVector m_PlotChannels;
  Index m_PlotSizeX;
  Index m_PlotSizeY;
  Index m_PlotShowFit;
  Index m_PlotSplitX;
  Index m_PlotSplitY;

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

  // Tick offset tool.
  const TimeOffsetTool* m_tickOffsetTool;

  // This subclass carries the state for this tool, i.e. data that can change
  // after initialization.
  class State {
  public:
    HistVectorMap ChannelTickModHists;
  };

  // Return the state.
  using StatePtr = std::shared_ptr<State>;
  StatePtr m_state;
  State& state() const { return *m_state; }

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, Index itkm) const;

  // Fill the histogram for a tickmod.
  int fillChannelTickMod(const AdcChannelData& acd, Index itkm0, Index itkm) const;

};

DEFINE_ART_CLASS_TOOL(AdcTickModViewer)

#endif
