// AdcDetectorPlotter.h

// David Adams
// March 2018
//
// Tool to make detector-wide event displays of data from an ADC channel data maps.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   WireAngle - Include wires with this angle (for protoDUNE, 0, +/-0.623)
//   DataType - Which data to plot: 0=prepared, 1=raw-pedestal
//   Tick0 - Tick used ast t = 0 for drift calculations.
//   DriftSpeed - Drift speed in cm/tick.
//   XMin, XMax - Limits for the drift coordinate
//   ZMin, ZMax - Limits for the wire coordinate
//   SignalThreshold - Signals in a channel-tick bin above this value are plotted
//   FirstTick - First tick number to display
//   LastTick - Last+1 tick number to display
//   ShowWires - Also show anode wires on the plot.
//   ShowCathode - Also show cathode planes (one point for each wire) on the plot.
//   ShowTpcSets - If not empty, only show wires and cathodes for these TPC sets.
//   ShowGrid - Also show (Root default) grid.
//   Title - Title for the plot.
//   FileName - Name for output plot file.
//              If blank, no file is written.
//              Existing file with the same name is replaced.
// For the title and file names, the following sustitutions are made:
//     %RUN%    --> run number
//     %SUBRUN% --> subrun number
//     %EVENT%  --> event number
//     %PAT%    --> pattern passed in call to view
//
// This is an example of a stateful tool, i.e. one that carries state beyond
// its configuration and configuration-derived data. This state is all held
// in the subclass State which is mutable so the view method can remain const.
// In a multithreaded environment, we will need to guard against access
// from multiple threads, e.g. may want a separate data thread.
//
// The tool creates a 2D plot of wire vs. drift coordinate and adds a point for
// each selected wire and its reflection on the cathode. It also adds points for
// each ADC wire-tick point with signal above threshold. The ADC data from the
// preceding calls are retained if the event (run, subrun and event number) are
// unchanged.

#ifndef AdcDetectorPlotter_H
#define AdcDetectorPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/Geometry/WireSelector.h"
#include <memory>

namespace geo {
  class GeometryCore;
}

class AdcChannelStringTool;

class AdcDetectorPlotter : public AdcChannelTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using TPadManipulatorPtr = std::unique_ptr<TPadManipulator>;

  class State {
  public:
    geo::GeometryCore* pgeo =nullptr;  // Geometry
    Index jobCount =0;                 // # calls for the job
    Index reportCount =0;              // # calls for the current report
    Index channelCount =0;             // # channels in the current report
    Index run =0;
    Index subrun =0;
    Index event =0;
    WireSelector sel;
    TPadManipulatorPtr ppad;
    std::string ofname;
  };

  using StatePtr = std::shared_ptr<State>;

  AdcDetectorPlotter(fhicl::ParameterSet const& ps);

  ~AdcDetectorPlotter() override =default;

  // AdcChannelTool interface.
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

  int addChannel(const AdcChannelData& acd, double xfac) const;

  // Return the state.
  // Shared pointer so we can make sure only one reference is out at a time.
  StatePtr getState() const { return m_state; }
  //StatePtr getState() const { return const_cast<StatePtr&>(m_state); }

  // Initialize the tool.
  void initialize();

private:

  // Configuration data.
  int            m_LogLevel;
  float          m_WireAngle;
  int            m_DataType;
  float          m_Tick0;
  float          m_DriftSpeed;
  float          m_XMin;
  float          m_XMax;
  float          m_ZMin;
  float          m_ZMax;
  float          m_SignalThreshold;
  Index          m_FirstTick;
  Index          m_LastTick;
  bool           m_ShowWires;
  bool           m_ShowCathode;
  IndexVector    m_ShowTpcSets;
  bool           m_ShowGrid;
  std::string    m_Title;
  std::string    m_FileName;

  StatePtr m_state;

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

};

DEFINE_ART_CLASS_TOOL(AdcDetectorPlotter)

#endif
