// AdcChannelMetric.h

// David Adams
// April 2018
//
// Tool to evalute metrics for single ADC channel and make histograms
// for multiple channels. Subclasses may be used to extend the list of
// metrics (names and algorithms).
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   Metric - Name of the metric to plot:
//              pedestal 
//              pedestalRms
//   FirstChannel - First channel to display
//   LastChannel - Last+1 channel to display
//   MetricMin - Minimum for the metric axis.
//   MetricMax - Maximum for the metric axis.
//   ChannelLineModulus - Repeat spacing for horizontal lines
//   ChannelLinePattern - Pattern for horizontal lines
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
// Drawings may include vertical lines intended to show boundaries of APAs,
// FEMBs, wire planes, etc.
//
// Lines are draw at N*ChannelLineModulus + ChannelLinePattern[i] for any
// integer N and any value if i in range of the array which are within
// the drawn channel range.
// If ChannelLineModulus is zero, then lines are drawn for the channels in
// ChannelLinePattern.

#ifndef AdcChannelMetric_H
#define AdcChannelMetric_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <vector>

class AdcChannelMetric : AdcChannelTool {

public:

  using Name = std::string;
  using Index = unsigned int;
  using IndexVector = std::vector<Index>;

  AdcChannelMetric(fhicl::ParameterSet const& ps);

  ~AdcChannelMetric() override =default;

  // Tool interface.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

  // Local method that directly returns the metric value and units.
  // Subclasse may overrride this to add metrics. They are expected to
  // call the fcl ctor of this class.
  virtual int
  getMetric(const AdcChannelData& acd, float& metricValue, Name& metricUnits) const;

private:

  // Configuration data.
  int            m_LogLevel;
  Name           m_Metric;
  Index          m_FirstChannel;
  Index          m_LastChannel;
  float          m_MetricMin;
  float          m_MetricMax;
  Index          m_ChannelLineModulus;
  IndexVector    m_ChannelLinePattern;
  Name           m_HistName;
  Name           m_HistTitle;
  Name           m_PlotFileName;
  Name           m_RootFileName;

};

DEFINE_ART_CLASS_TOOL(AdcChannelMetric)

#endif
