// AdcDataPlotter.h

// David Adams
// July 2017
//
// Tool to make event displays of prepared data from an ADC channel data map.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   DataType - Which data to plot: 0=prepared, 1=raw-pedestal
//   FirstTick - First tick number to display
//   LastTick - Last+1 tick number to display
//   FirstChannel - First channel to display
//   LastChannel - Last+1 channel to display
//   MaxSignal - Displayed signal range is (-MaxSignal, MaxSignal)
//   ChannelLineSpacing - Spacing between horizontal lines (start at 0)
//   HistName - Histogram name (should be unique within Root file)
//   HistTitle - Histogram title
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

#ifndef AdcDataPlotter_H
#define AdcDataPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AdcDataPlotter : AdcChannelTool {

public:

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
  unsigned int   m_FirstChannel;
  unsigned int   m_LastChannel;
  double         m_MaxSignal;
  unsigned int   m_ChannelLineSpacing;
  int            m_Palette;
  std::string    m_HistName;
  std::string    m_HistTitle;
  std::string    m_PlotFileName;
  std::string    m_RootFileName;

};

DEFINE_ART_CLASS_TOOL(AdcDataPlotter)

#endif
