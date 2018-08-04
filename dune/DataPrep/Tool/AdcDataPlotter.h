// AdcDataPlotter.h

// David Adams
// July 2017
//
// Tool to make event displays of prepared data from an ADC channel data map.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   DataType - Which data to plot: 0=prepared, 1=raw-pedestal, 2=signal
//   TickRange - Name of the tick range used in the display
//               The name must be defined in the IndexRangeTool tickRanges
//               If blank or not defined, the full range is used.
//   FirstTick - First tick number to display
//   LastTick - Last+1 tick number to display
//   ChannelRanges - Names of channel ranges to display.
//                   Ranges are obtained from the tool channelRanges.
//                   Special name "" or "data" plots all channels in data with label "All data".
//                   If the list is empty, data are plotted.
//   FembTickOffsets - Tick offset for each FEMB. FEMB = (offline channel)/128
//                     Offset is zero for FEMBs beyond range.
//                     Values should be zero (empty array) for undistorted plots
//   OnlineChannelMapTool - Name of tool mapping channel # to online channel #.
//   MaxSignal - Displayed signal range is (-MaxSignal, MaxSignal)
//   ChannelLineModulus - Repeat spacing for horizontal lines
//   ChannelLinePattern - Pattern for horizontal lines
//   HistName - Histogram name (should be unique within Root file)
//   HistTitle - Histogram title
//   PlotSizeX, PlotSizeY: Size in pixels of the plot file.
//                         Root default (700x500?) is used if either is zero.
//   PlotFileName - Name for output plot file.
//                  If blank, no file is written.
//                  Existing file with the same name is replaced.
//   RootFileName - Name for the output root file.
//                  If blank, histograms are not written out.
//                  Existing file with the same is updated.
// For the title and file names, the following sustitutions are made:
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

#ifndef AdcDataPlotter_H
#define AdcDataPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <vector>
#include <memory>

class AdcChannelStringTool;
class IndexMapTool;
class IndexRangeTool;

class AdcDataPlotter : AdcChannelTool {

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
  std::string    m_TickRange;
  NameVector     m_ChannelRanges;
  IntVector      m_FembTickOffsets;
  std::string    m_OnlineChannelMapTool;
  double         m_MaxSignal;
  Index          m_ChannelLineModulus;
  IndexVector    m_ChannelLinePattern;
  int            m_Palette;
  std::string    m_HistName;
  std::string    m_HistTitle;
  Index          m_PlotSizeX;
  Index          m_PlotSizeY;
  std::string    m_PlotFileName;
  std::string    m_RootFileName;

  // Derived configuration data.
  IndexRange m_tickRange;

  // Channel ranges.
  IndexRangeVector m_crs;

  // Client tools.
  const AdcChannelStringTool* m_adcStringBuilder;
  const IndexMapTool* m_pOnlineChannelMapTool;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, const IndexRange& ran) const;

};

DEFINE_ART_CLASS_TOOL(AdcDataPlotter)

#endif
