// AdcDataPlotter.h

// David Adams
// July 2017
//
// Tool to make channel vs. tick displays of data from an ADC channel data map.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   DataType - Which data to plot: 0=prepared, 1=raw-pedestal, 2=signal, 3=not signal
//   DataView - Which view to use: "" for top, xxx merges everything from xxx
//   TickRange - Name of the tick range used in the display
//               The name must be defined in the IndexRangeTool tickRanges
//               If blank or not defined, the full range is used.
//   TickRebin - If > 1, histo bins include this # ticks.
//   ChannelRanges - Names of channel ranges to display.
//                   Ranges are obtained from the tool channelRanges.
//                   Special name "" or "data" plots all channels in data with label "All data".
//                   If the list is empty, data are plotted.
//   ClockFactor - Clock to tick conversion factor (0.04 for protoDUNE).
//   ClockOffset - Clock offset between trigger and nominal tick 0.
//   FembTickOffsets - Tick offset for each FEMB. FEMB = (offline channel)/128
//                     Offset is zero for FEMBs beyond range.
//                     Values should be zero (empty array) for undistorted plots
//   OnlineChannelMapTool - Name of tool mapping channel # to online channel #.
//   MinSignal - Formula for min signal. If absent, -(max signal) is used.
//   MaxSignal - Formula for max signal. Displayed signal range is (min signal, max signal)
//   SkipBadChannels - If true, skip channels flagged as bad.
//   EmptyColor - If >=0, empty bins are drawn in this color (See TAttColor).
//                Otherwise empty bins are drawn with value zero.
//                Bins may be empty if a channel is nor processed, if a tick out of range
//                or a tick is not selected (outside ROI) for DataType 2.
//                EmptyColor is not used when rebinning.
//   ChannelLineModulus - Repeat spacing for horizontal lines
//   ChannelLinePattern - Pattern for horizontal lines
//   HistName - Histogram name (should be unique within Root file)
//   HistTitle - Histogram title (appears above histogram)
//   PlotTitle - Plot title (appears below histogram an only on plots)
//   PlotSizeX, PlotSizeY: Size in pixels of the plot file.
//                         Root default (700x500?) is used if either is zero.
//   PlotFileName - Name for output plot file.
//                  If blank, no file is written.
//                  Existing file with the same name is replaced.
//   RootFileName - Name for the output root file.
//                  If blank, histograms are not written out.
//                  Existing file with the same is updated.
// For the title and file names, substitutions are made with adcStringBuilder, e.g.
//     %RUN%    --> run number
//     %SUBRUN% --> subrun number
//     %EVENT%  --> event number
//     %CHAN1%  --> First channel number 
//     %CHAN2%  --> Last channel number 
// Drawings may include horizontal lines intended to show boundaries of APAs,
// FEMBs, wire planes, etc.
//
// Lines are draw at N*ChannelLineModulus + ChannelLinePattern[i] for any
// integer N and any value if i in range of the array which are within
// the drawn channel range.
// If ChannelLineModulus is zero, then lines are drawn for the channels in
// ChannelLinePattern.
//
// If FirstChannel < LastChannel, then only channels in that range are displayed
// and no histogram is produced if the passed data has no channels in the range.
//
// If ClockFactor > 0, then tick + ClockFactor*(channelClock - triggerClock + ClockOffset)
// is used in place of tick.

#ifndef AdcDataPlotter_H
#define AdcDataPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Utility/ParFormula.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <vector>
#include <memory>

class AdcChannelStringTool;
class IndexMapTool;
class IndexRangeTool;
namespace lariov {
  class ChannelStatusProvider;
}
class RunDataTool;

class AdcDataPlotter : TpcDataTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IntVector = std::vector<int>;
  using IndexRangeVector = std::vector<IndexRange>;
  using Name = std::string;
  using NameVector = std::vector<Name>;

  AdcDataPlotter(fhicl::ParameterSet const& ps);

  ~AdcDataPlotter() override =default;

  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

private:

  // Configuration data.
  int            m_LogLevel;
  int            m_DataType;
  Name           m_DataView;
  Name           m_TickRange;
  Index          m_TickRebin;
  NameVector     m_ChannelRanges;
  float          m_ClockFactor;
  float          m_ClockOffset;
  IntVector      m_FembTickOffsets;
  Name           m_OnlineChannelMapTool;
  ParFormula*    m_MinSignal;
  ParFormula*    m_MaxSignal;
  bool           m_SkipBadChannels;
  Index          m_EmptyColor;
  Index          m_ChannelLineModulus;
  IndexVector    m_ChannelLinePattern;
  int            m_Palette;
  Name           m_HistName;
  Name           m_HistTitle;
  Name           m_PlotTitle;
  Index          m_PlotSizeX;
  Index          m_PlotSizeY;
  Name           m_PlotFileName;
  Name           m_RootFileName;

  // Derived configuration data.
  IndexRange m_tickRange;
  bool m_needRunData;

  // Channel ranges.
  IndexRangeVector m_crs;

  // Client tools and services.
  const AdcChannelStringTool* m_adcStringBuilder;
  const IndexMapTool* m_pOnlineChannelMapTool;
  const lariov::ChannelStatusProvider* m_pChannelStatusProvider;
  const RunDataTool* m_prdtool;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, const IndexRange& ran) const;

};


#endif
