// AdcChannelMetric.h

// David Adams
// April 2018
//
// Tool to evalute metrics for single ADC channel and make histograms
// of metric vs. channel for ranges of channels.
//
// If plots are made, graphs are shown instead of histograms.
// If a plot range is specified then values outside the range are
// shown at the nearest range limit.
//
// Subclasses may be used to extend the list of
// metrics (names and algorithms).
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   Metric - Name of the plotted metric. This can be the name of any
//            metadata field or any of the following:
//              pedestal 
//              pedestalDiff - pedestal - (pedestal for first event)
//              pedestalRms - pedestal noise from the ADC channel data (e.g. filled by pedestal finder)
//              fembID [0, 120) in protoDUNE
//              apaFembID - FEMB number in the APA [0, 20)
//              fembChannel - channel # in the FEMB [0, 128)
//              rawRms - RMS of (ADC - pedestal)
//              rawTailFraction - Fraction of ticks with |raw - ped| > 3*noise
//   ChannelRanges - Channel ranges for which metric channel histograms and plots are made.
//                   Ranges are obtained from the tool channelRanges.
//                   Special name "all" or "" plots all channels with label "All".
//                   If the list is empty, all are plotted.
//   MetricMin - Minimum for the metric axis.
//   MetricMax - Maximum for the metric axis.
//   MetricBins - If nonzero, # channels vs metric is plotted with this binning instead of
//                metric vs channel.
//   ChannelLineModulus - Repeat spacing for horizontal lines
//   ChannelLinePattern - Pattern for horizontal lines
//   HistName - Histogram name (should be unique within Root file)
//              If the HistName name does not include "EVENT%", then only summary histogram
//              and plot are created.
//   HistTitle - Histogram title
//   MetricLabel - Histogram label for the metric axis
//   PlotSizeX, PlotSizeY: Size in pixels of the plot file.
//                         Root default (700x500?) is used if either is zero.
//   PlotFileName - Name for output plot file.
//                  If blank, no file is written.
//                  Existing file with the same name is replaced.
//   PlotUsesStatus - If nonzero plot colors indicate channel status (good, bad noisy).
//   RootFileName - Name for the output root file.
//                  If blank, histograms are not written out.
//                  Existing file with the same is updated.
// For the title and file names, the following sustitutions are made:
//     %RUN%    --> run number
//     %SUBRUN% --> subrun number
//     %EVENT%  --> event number
//     %CHAN1%  --> First channel number 
//     %CHAN2%  --> Last channel number 
//     %CRNAME%  --> Channel range name
//     %CRLABEL%  --> Channel range label
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
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <vector>

class AdcChannelStringTool;
namespace lariov {
  class ChannelStatusProvider;
}

class AdcChannelMetric : AdcChannelTool {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexRangeVector = std::vector<IndexRange>;

  AdcChannelMetric(fhicl::ParameterSet const& ps);

  ~AdcChannelMetric() override;

  // Initialize this tool.
  // We cache the status for each channel.
  // Does nothing if already called unlesss force is true.
  void initialize(bool force =false);

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
  NameVector     m_ChannelRanges;
  IndexVector    m_ChannelCounts;
  float          m_MetricMin;
  float          m_MetricMax;
  Index          m_MetricBins;
  Index          m_ChannelLineModulus;
  IndexVector    m_ChannelLinePattern;
  Name           m_HistName;
  Name           m_HistTitle;
  Name           m_MetricLabel;
  Index          m_PlotSizeX;
  Index          m_PlotSizeY;
  Name           m_PlotFileName;
  int            m_PlotUsesStatus;
  Name           m_RootFileName;

  // Channel ranges.
  IndexRangeVector m_crs;
  
  // Flag indicating separate plots should be made based on status.
  bool m_doSummary;

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

  // Channel status provider.
  const lariov::ChannelStatusProvider* m_pChannelStatusProvider;

  // Create the plot for one range.
  DataMap viewMapForOneRange(const AdcChannelDataMap& acds, const IndexRange& ran) const;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, const IndexRange& ran) const;

  // Summary data for one channel.
  class MetricSummary {
  public:
    Index count = 0;
    double sum = 0.0;
    double sumsq = 0.0;
    // Add an entry.
    void add(double val) { ++count; sum+=val; sumsq+=val*val; }
    double mean() const { return count ? sum/count : 0.0; }
    double meansq() const { return count ? sumsq/count : 0.0; }
    double rms() const {
      double valm = mean();
      double arg = meansq() - valm*valm;
      return arg > 0 ? sqrt(arg) : 0.0;
    }
    double dmean() const { return count ? rms()/sqrt(double(count)) : 0.0; }
  };

  class Metric {
  public:
    float value =0.0;
    float error =0.0;
    void setValue(float a_value) { value = a_value; }
    void setError(float a_error) { error = a_error; }
    bool operator<(const Metric& rhs) const { return value < rhs.value; }
  };
  using MetricMap = std::map<Index, Metric>;
  using MetricSummaryVector = std::vector<MetricSummary>;
  using MetricSummaryMap = std::map<IndexRange, MetricSummaryVector>;

  // This subclass carries the state for this tool, i.e. data that can change
  // after initialization.
  class State {
  public:
    Index initCount = 0;
    Index callCount = 0;
    Index firstRun =0;
    Index lastRun =0;
    Index firstEvent =0;
    Index lastEvent =0;
    Index eventCount =0;
    Index runCount =0;
    MetricSummaryMap crsums;    // Summary for each channel.
    MetricMap pedRefs;          // Reference pedestal for each channel;
    void update(Index run, Index event);
    IndexVector channelStatuses;
  };

  // Shared pointer so we can make sure only one reference is out at a time.
  using StatePtr = std::shared_ptr<State>;
  StatePtr m_state;

  // Return the state.
  State& getState() const { return *m_state; }

  // Local method that fills the metric histogram and creates plots
  // for the provided range and data.
  void processMetricsForOneRange(const IndexRange& ran, const MetricMap& mets, TH1* ph,
                                 Name ofname, Name ofrname, bool useErrors) const;

  // Local method to create the histogram.
  TH1* createHisto(const AdcChannelData& acd, const IndexRange& ran) const;

  // Return the status for a channel.
  //   0 - unknown (no status provider).
  //   1 - got (not bad and not noisy)
  //   2 - bad
  //   3 - good
  Index channelStatus(Index icha) const;

};

DEFINE_ART_CLASS_TOOL(AdcChannelMetric)

#endif
