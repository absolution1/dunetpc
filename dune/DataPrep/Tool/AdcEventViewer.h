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
//  EventGraphs - Array of graph specifiers with format xname:xmin:xmax:yname:ymin:ymax
//                The name fields may contain any of the above plus
//                  event - Event number
//                For both histograms and graphs min >= max gives auto scaling

#ifndef AdcEventViewer_H
#define AdcEventViewer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <iostream>

class AdcEventViewer : AdcChannelTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using HistVector = std::vector<TH1*>;
  using HistMap = std::map<Name, TH1*>;
  using NameMap = std::map<Name, Name>;
  using FloatVector = std::vector<float>;
  using GraphVector = std::vector<TGraph*>;

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
    FloatVector xvals;
    float xmin;
    float xmax;
    Name vary;
    FloatVector yvals;
    Name ylab;
    float ymin;
    float ymax;
    GraphInfo() { };
    GraphInfo(Name avarx, Name axlab, float axmin, float axmax,
              Name avary, Name aylab, float aymin, float aymax)
    : varx(avarx), xlab(axlab), xmin(axmin), xmax(axmax),
      vary(avary), ylab(aylab), ymin(aymin), ymax(aymax) { }
    // Add value to vector if name matches.
    void add(Name var, float val) {
      if ( var == varx ) xvals.push_back(val);
      if ( var == vary ) yvals.push_back(val);
    }
  };
  using GraphInfoVector = std::vector<GraphInfo>;

  // This subclass carries the state for this tool, i.e. data that can change
  // after initialization.
  class State {
  public:
    Index run;
    Index event;             // Current event
    IndexVector events;      // Events in processed order.
    IndexSet eventSet;       // Events ordered.
    Index ngroup;            // # groups processed for this event
    IndexSet fembIDSet;      // FEMBs for ths event.
    Index nchan;             // # channels processed for this event
    HistVector hists;        // Monitoring histograms.
    GraphInfoVector graphs;  // Monitoring graphs.
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

  // Shared pointer so we can make sure only one reference is out at a time.
  StatePtr m_state;

};

DEFINE_ART_CLASS_TOOL(AdcEventViewer)

#endif
