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
//            metadata field or any of the following single values:
//              pedestal 
//              pedestalDiff - pedestal - (pedestal reference)
//              pedestalRms - pedestal noise from the ADC channel data (e.g. filled by pedestal finder)
//              fembID [0, 120) in protoDUNE
//              apaFembID - FEMB number in the APA [0, 20)
//              fembChannel - channel # in the FEMB [0, 128)
//            or any of the following calculated values
//              rawRms - RMS of (ADC - pedestal)
//              samRms - RMS of sample
//              samRmsNN - RMS of integration over NN contiguous samples (NN = 1, 2, ...)
//              rawTailFraction - Fraction of ticks with |raw - ped| > 3*noise
//              sigFrac: Fraction of samples that are signal.
//              sigRms: RMS of the signal samples.
//              nsgRms: RMS of the not-signal samples.
//              nsgRmsNN: RMS of integration over NN contiguous not-signal samples.
//            In the case a data view other than the top (DataView not blank), the single values
//            are taken from the first entry. The calculated values include all entries.
//   DataView - Name of the data view to use.
//   PedestalReference - Name of the FloatArrayTool that holds the pedestal reference values.
//                       If the value is "first", the pedestal for the first event is used.
//   MetricSummaryView - If not empty and a summary is requested, this specifies the view
//                       that is plotted, this view of the metric summary is plotted.
//                       The format is is VVV or VVV:EEE where VVV=position and EEE=error
//                       can be any of the following. Default is "mean:dmean".
//                  count - Number of values
//                  mean - Mean of the value
//                  dmean - error on the mean = rms/sqrt(count)
//                  rms - RMS of the values
//                  drms - error on the RMS
//                  sdev - RMS from the mean of the values
//                  min - Minimum value
//                  max - Maximum value
//                  center - 0.5*(min+max)
//                  range - max - min
//                  halfRange - 0.5*range
//   ChannelRanges - Channel ranges for which metric channel histograms and plots are made.
//                   Ranges are obtained from the tool channelRanges.
//                   Special name "all" or "" plots all channels with label "All".
//                   If the list is empty, all are plotted.
//   MetricMin - Formula for the minimum for the metric axis.
//   MetricMax - Formula for the maximum for the metric axis.
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
//   MetadataFlags - Vector of any of none of the following:
//                     write - Write value as ADC channel metadata.
//                     warnpresent - Warn if metadata field is present before write.
//                     read - read value from metadata (instead of calculation) if present
//                     warnabsent - Warn if requested metatdata field is not present
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
#include "dune/DuneInterface/Utility/ParFormula.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <vector>

class AdcChannelStringTool;
namespace lariov {
  class ChannelStatusProvider;
}

class FloatArrayTool;
class RunDataTool;

class AdcChannelMetric : TpcDataTool {

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
  DataMap update(AdcChannelData& acd) const override;
  DataMap updateMap(AdcChannelDataMap& acds) const override;

  // Local method that directly returns the metric value and units.
  // Subclasses may overrride this to add metrics. They are expected to
  // call the fcl ctor of this class.
  // The weight is used when averaging over events, e.g. the number
  // of samples or ROIs contributing to the metric value.
  virtual int
  getMetric(const AdcChannelData& acd, Name met, double& metricValue,
            Name& metricUnits, double& metricWeight) const;

private:

  // Configuration data.
  int            m_LogLevel;
  Name           m_Metric;
  Name           m_DataView;
  Name           m_PedestalReference;
  Name           m_MetricSummaryView;
  NameVector     m_ChannelRanges;
  IndexVector    m_ChannelCounts;
  ParFormula*    m_MetricMin;
  ParFormula*    m_MetricMax;
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
  NameVector     m_MetadataFlags;

  // Derived data.
  bool m_mdRead;
  bool m_mdWrite;
  bool m_mdWarnAbsent;
  bool m_mdWarnPresent;
  RunDataTool* m_prdtool;

  // Channel ranges.
  IndexRangeVector m_crs;
  
  // Summary info.
  bool m_doSummary;
  bool m_doSummaryError;
  Name m_summaryValue;
  Name m_summaryError;

  // Pedestal reference tool.
  const FloatArrayTool* m_pPedestalReference;

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

  // Channel status provider.
  const lariov::ChannelStatusProvider* m_pChannelStatusProvider;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, const IndexRange& ran) const;

  // Summary data for one channel.
  // Calculates mean, RMS, their uncertainties and more from accumulated data.
  // Note that the RMS may be more appropriate for RMS-like variables like noise
  // as it averages squares instead of values.
  //
  // The value for each event is added with a weight that is used in the calculation
  // of the mea, RMS and other values. Only the relative values of the weights
  // affects these calculations.
  //
  // Uncertainties are evaluated by dividing the variance by sqrt(Neff) where Neff
  // is the effective number independent measurements. Its value depends on
  // the the weight flag:
  //   0 - Neff = # events with weight > 0
  //   1 - Neff = sum of weights
  class MetricSummary {
  public:
    Index weightFlag = 0;
    Index eventCount = 0;
    Index weightedEventCount = 0;
    double weightSum;
    double sum = 0.0;
    double sumsq = 0.0;
    double minval = 0.0;
    double maxval = 0.0;
    // Add an entry.
    void add(double val, double weight) {
      ++eventCount;
      if ( weight <= 0.0 ) return;
      ++weightedEventCount;
      if ( weightSum == 0 || val < minval ) minval = val;
      if ( weightSum == 0 || val > maxval ) maxval = val;
      weightSum += weight;
      sum += weight*val;
      sumsq += weight*val*val;
    }
    double neff() const { return weightFlag ? weightSum : weightedEventCount; }
    double mean() const { return weightSum ? sum/weightSum : 0.0; }
    double dmean() const { return weightSum > 0.0 ? sdev()/sqrt(neff()) : 0.0; }
    double meansq() const { return weightSum ? sumsq/weightSum : 0.0; }
    double rms() const {
      return sqrt(meansq());
    }
    double drms() const {
      if ( weightSum <= 0.0 ) return 0.0;
      double rmsVal = rms();
      double rmsVar = meansq() + rmsVal*(rmsVal - 2.0*mean());
      return rmsVar > 0.0 ? sqrt(rmsVar/neff()) : 0.0;
    }
    double sdev() const {
      double valm = mean();
      double arg = meansq() - valm*valm;
      return arg > 0 ? sqrt(arg) : 0.0;
    }
    double center() const { return 0.5*(minval + maxval); }
    double range() const { return  maxval - minval; }
    // Return if the provided string is a value name.
    static bool isValueName(Name vnam) {
      const std::set<Name> sumVals =
        {"weightFlag", "eventCount", "weightedEventCount", "weightSum",
         "mean", "rms", "sdev", "min", "max", "dmean", "drms",
         "center", "range", "halfRange"};
      return sumVals.find(vnam) != sumVals.end();
    }
    // Return a value by name.
    double getValue(Name vnam) const {
      if ( vnam == "weightFlag" ) return weightFlag;
      if ( vnam == "weightedEventCount" ) return weightedEventCount;
      if ( vnam == "eventCount" ) return eventCount;
      if ( vnam == "weightSum" ) return weightSum;
      if ( vnam == "mean" ) return mean();
      if ( vnam == "rms" ) return rms();
      if ( vnam == "sdev" ) return sdev();
      if ( vnam == "min" ) return minval;
      if ( vnam == "max" ) return maxval;
      if ( vnam == "center" ) return center();
      if ( vnam == "range" ) return range();
      if ( vnam == "halfRange" ) return 0.5*range();
      if ( vnam == "dmean" ) return dmean();
      if ( vnam == "drms" ) return drms();
      return 0.0;
    }
  };

  class Metric {
  public:
    double value =0.0;
    double error =0.0;
    void setValue(double a_value) { value = a_value; }
    void setError(double a_error) { error = a_error; }
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
    bool update(Index run, Index event);
    IndexVector channelStatuses;
    float metricMin =0;
    float metricMax =100;
  };

  // Shared pointer so we can make sure only one reference is out at a time.
  using StatePtr = std::shared_ptr<State>;
  StatePtr m_state;

  // Return the state.
  State& getState() const { return *m_state; }

  // Tool interface plus map to hold results.
  // The update methods can attach the metrics to the channel data.
  // The view methods ignore those reults.
  DataMap viewMapLocal(const AdcChannelDataMap& acds, MetricMap& mets) const;

  // Create the plot for one range.
  DataMap viewMapForOneRange(const AdcChannelDataMap& acds, const IndexRange& ran,
                             MetricMap& mets) const;

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

  // Evaluate the metric formulas and store results in the state.
  void evaluateFormulas() const;

};


#endif
