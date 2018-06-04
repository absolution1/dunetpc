// AdcDataPlotter.h

// David Adams
// July 2017
//
// Tool to make event displays of prepared data from an ADC channel data map.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   DataType - Which data to plot: 0=prepared, 1=raw-pedestal, 2=signal
//   FirstTick - First tick number to display
//   LastTick - Last+1 tick number to display
//   FirstChannel - First channel to display
//   LastChannel - Last+1 channel to display
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
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <vector>

class AdcChannelStringTool;

class AdcDataPlotter : AdcChannelTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;

  AdcDataPlotter(fhicl::ParameterSet const& ps);

  ~AdcDataPlotter() override =default;

  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

private:

  // Configuration data.
  int            m_LogLevel;
  int            m_DataType;
  unsigned long  m_FirstTick;
  unsigned long  m_LastTick;
  Index          m_FirstChannel;
  Index          m_LastChannel;
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

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

};

DEFINE_ART_CLASS_TOOL(AdcDataPlotter)

#endif
