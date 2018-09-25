// AdcEventViewer.h
//
// David Adams
// August 2018
//
// ADC tool to extract and display event-level information.
//
// Configuration:
//     LogLevel - Logging level: 0=none, 1=init, 2=call, ...
//   EventHists - Array of string histogram specifiers with format "name:nbin:xmin:xmax"
//                The field name must contain the histogramed variable. Allowed values:
//                  nfemb - # FEMBs with data 
//                  rmPedPower - sqrt(<(ped noise)^2>)
//                  meanPedPower - <ped>
//  EventGraphs - Array of graph specifiers with format xname:xmin:xmax:yname:ymin:ymax
//                The name fields may contain any of the above plus
//                  event - Event number
//                  clock - Timing clock count
//                For both histograms and graphs min >= max gives auto scaling
//  ChannelRanges - If this is not empty, then a separate histogram an graph plots
//                  are made for each channel range. Otherwise all channels are included.
//                  The tool channelRanges is used to map the names in this list to ranges.
//                  Special name "all" or "" plots all channels with label "All".
//  ChannelRangeLabel - Label for channel ranges. May include %XXX% for
//                      XXX = CRLABEL to replace with cr.label()
//                      XXX = CRLABEL1 to replace with cr.label(1)
//                      XXX = CRLABEL2 to replace with cr.label(2)
//      ClockUnit - Unit for plotting clock counts: tick, ktick or Mtick
//      ClockRate - # clock ticks per second. Used to conver to time.
#ifndef AdcEventViewer_H
#define AdcEventViewer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include <iostream>

class AdcEventViewer : AdcChannelTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;
  using LongIndex = unsigned long;
  using LongIndexVector = std::vector<LongIndex>;
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using HistVector = std::vector<TH1*>;
  using HistMap = std::map<Name, TH1*>;
  using NameMap = std::map<Name, Name>;
  using FloatVector = std::vector<float>;
  using GraphVector = std::vector<TGraph*>;
  using IndexRangeVector = std::vector<IndexRange>;

  // Subclass that associates a variable name with a histogram.
  //  vary != "" ==> 2D histo
  class HistInfo {
  public:
    TH1* ph = nullptr;
    Name varx;
    Name vary;
  };

  using HistInfoMap = std::map<Name, HistInfo>;

  // Subclass that describes a graph.
  struct GraphInfo {
    Name varx;
    Name xlab;
    Name xunit;
    FloatVector xvals;
    float xmin;
    float xmax;
    Name vary;
    FloatVector yvals;
    Name ylab;
    Name yunit;
    float ymin;
    float ymax;
    GraphInfo() { };
    GraphInfo(Name avarx, Name axlab, Name axunit, float axmin, float axmax,
              Name avary, Name aylab, Name ayunit, float aymin, float aymax)
    : varx(avarx), xlab(axlab), xunit(axunit), xmin(axmin), xmax(axmax),
      vary(avary), ylab(aylab), yunit(ayunit), ymin(aymin), ymax(aymax) { }
    // Add value to vector if name matches.
    void add(Name var, float val) {
      if ( var == varx ) xvals.push_back(val);
      if ( var == vary ) yvals.push_back(val);
    }
    // Return axis labels.
    Name xAxisLabel() const {
      return xlab + (xunit.size() ? " [" + xunit + "]" : "");
    }
    Name yAxisLabel() const {
      return ylab + (yunit.size() ? " [" + yunit + "]" : "");
    }
  };
  using GraphInfoVector = std::vector<GraphInfo>;

  // This subclass carries state data for one channel range.
  class ChannelRangeState {
  public:
    IndexSet fembIDSet;      // FEMBs for ths event.
    Index nchan =0;          // # channels processed for this event
    float pedSum =0.0;       // Sum over pedestals
    float pedPower =0.0;     // Sum over (pedestal noise)^2
    HistVector hists;        // Monitoring histograms.
    GraphInfoVector graphs;  // Monitoring graphs.
  };

  using ChannelRangeStates = std::map<Name,ChannelRangeState>;

  // This subclass carries the state for this tool, i.e. data that can change
  // after initialization.
  class State {
  public:
    Index run;
    Index event;             // Current event.
    LongIndex clock;         // Current timing clock.
    IndexVector events;      // Events in processed order.
    LongIndexVector clocks;  // Timing clocks in processed order.
    LongIndex firstClock;    // First timing clock.
    LongIndex minClock;      // Minimum timing clock.
    IndexSet eventSet;       // Events ordered.
    Index ngroup;            // # groups processed for this event
    ChannelRangeStates crstates;
  };

  using StatePtr = std::shared_ptr<State>;

  AdcEventViewer(fhicl::ParameterSet const& ps);

  ~AdcEventViewer() override;

  // AdcChannelTool methods.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

  // Return the state.
  State& state() const { return *m_state; }

  // Return the state for a channel range.
  ChannelRangeState& crstate(Name crn) const;

  // Initialize the state for a new event.
  void startEvent(const AdcChannelData& acd) const;

  // End current event in state.
  void endEvent() const;

  // Print a report.
  void printReport() const;

  // Display the histograms, i.e create pngs.
  void displayHists() const;

  // Display the graphs, i.e create pngs.
  void displayGraphs() const;

private:

  // Configuration data.
  int m_LogLevel;
  NameVector m_EventHists;
  NameVector m_EventGraphs;
  NameVector m_ChannelRanges;
  Name m_ChannelRangeLabel;
  Name m_ClockUnit;
  double m_ClockRate;

  // Channel ranges.
  IndexRangeVector m_crs;
  
  // Shared pointer so we can make sure only one reference is out at a time.
  StatePtr m_state;

};

DEFINE_ART_CLASS_TOOL(AdcEventViewer)

#endif
