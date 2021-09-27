// AdcRoiViewer.h
//
// David Adams
// October 2017
//
// Tool to extract information about the ROIs in an ADC channel.
//
// Configuration:
//            LogLevel - Logging level: 0=none, 1=init, 2=call, ...
//           SigThresh - if <0, then only keep ROIs with a tick below this value
//                       if >0, then only keep ROIs with a tick above this value
//          TickBorder - if > 0, only keep ROIs forwhich there are this many ticks or more
//                       before the start and after the end of the ROI.
//          RoiHistOpt - histo option:  0 - No histograms
//                                      1 - sample vs. tick
//                                      2 - raw vs. tick
//                                     10 + i - As above for i except vs. tick - tick0
//              FitOpt - ROI fitting option
//                         0 - no fit
//                         1 - fit with coldelecReponse
//          RoiPlotOpt - 0 = none, 1 = for separate for each event, 2 = multi-event
//           StartTime - Offset for time meaurements in sec since 1970.
//    PulserStepCharge - Charge per unit step in a pulser run
//     PulserDacOffset - Offset in pulser: Qin = PulserStepCharge*(DAC - PulserDacOffset)
//    PulserChargeUnit - Unit for the pulser charge (ke, fC, ...)
//         MaxRoiPlots - Maximum # of ROI plots (<0 is no limit, 0 is no plots)
//         RoiPlotPadX - Number of pad columns in ROI plots. No plots if 0.
//         RoiPlotPadY - Number of pad rows in ROI plots. No plots if 0.
//            SumHists - Array of summary histogram specifiers. See below.
//           SumNegate - If true, the following variable replacements are made for all sum hists:
//                         fitHeight --> fitHeightNeg
//         SumPlotPadX - Number of pad columns in summary plots.
//         SumPlotPadY - Number of pad rows in summary plots.
//       ChannelRanges - Ranges of channels for channel summary plots.
//                       Obtained from IndexRangeTool channelRanges.
//        ChanSumHists - Array of specifiers for the channel summary histograms.
//  ChannelLineModulus - Repeat spacing for horizontal lines in summary plots
//  ChannelLinePattern - Pattern for horizontal lines in summary plots
//         RunDataTool - Name for the run data tool. If found and pulser is on, then each
//                       ROI is assigned a charge corresponding to the pulser setting.
//      TickOffsetTool - Name of the tool that provides the tick offset.
//     RoiRootFileName - Name of file to which the ROI histograms are copied.
//     SumRootFileName - Name of file to which the summary histograms are copied.
// ChanSumRootFileName - Name of file to which the channel summary histograms are copied.
//          PlotLabels - Array of strings that may be used whe constructin plot titles.
//                       The are referencesd as %LAB0%, %LAB1%, ... %LAB9%
//
// Summary histograms
// ------------------
// Summary histograms show the number of ROIs in bins of some ROI variable such
// as tick or pulse height.
// A summary histogram specifier is a parameter set with the following fields:
//     var: Name of the variable to draw:
//            sigArea - Raw area
//            sigWidth - Raw width
//            fitHeight - Height from ROI fit
//            fitHeightNeg - Negative of height from ROI fit
//            fitHeightGain - Height/(pulser charge)
//            fitWidth - Shaping time from ROI fit
//            fitPos - Postion [ticks] from ROI fit
//            fitPosRem - remainder(fitPos, 1)
//            fitToffPulser - fmod(fitPos + toff, Tpulser) where
//              toff is the timing offset, i.e. the time for tick 0 and
//              Tpulser is the period of the pulser (e.g. 497 ticks)
//            fitToffPulserMod10 - fmod(fitToffPulser, 10)
//            fitChiSquare - Chi-square from ROI fit (as reported by TF1)
//            fitChiSquareDof - Chi-square/DOF from ROI fit (both from TF1)
//            fitCSNorm - Chi-square/(ped RMS)^2
//            fitCSNormDof - Chi-square/DOF/(ped RMS)^2
//            timingPhase_fitToffPulserMod10 - 2D plot of timing phase (0 to 1) vs offset tick
//            timeSec, timeHour, timeDay - DAQ time since StartTime
//            procEvent - event number filled once per event
//            procTimeSec, procTimeHour, procTimeDay - DAQ time since StartTime filled once/event
//    name: Name of the histogram. Include %CHAN% to get separate histos for each channel
//   title: Histogram title
//    nbin: # bins
//    xmin: Lower edge of the first bin
//    xmax: Upper edge of the last bin
//     fit: Name of function used to fit the distribution, e.g. "gaus"
//    plot: Name of file where plot should be created. E.g. myvar%CHAN%.png
//    pwid: Plot only includes region of this width around the peak.
// If xmin < xmax and xmin > 0, the range will have width xmin centered on the median
// value for the first event. If xmax > 0, the lower edge is rounded to that value.
// If xmin <= xmax otherwise (e.g. xmin = xmax = 0), Root will do autoscaling of the axis.
// E.g.: {var:fitHeight name:"hfh_chan%0CHAN" title:"Fit height for channel %CHAN%"
//        nbin:50 xmin:0.0 xmax:5.0}
// Unless otherwise noted, the histograms are filled once for each ROI.
// For variable proc*, the histogram is filled once each time a channel is processed.
//
// Channel summary histograms
// --------------------------
// Channel summary histograms hold a metric for each channel derived from a summary histogram.
// The following fields specify a channel summary histogram:
//     name - name for the summary histogram (%CRNAME% is replaced with the channel range name)
//    title - histogram title  (substitutions for %CRLABEL%, %RUN%, ...)
//  valHist - Name of the summary histogram template from which the metric is derived (should include %CHAN%)
//  valType - Specifies the metric to be extracted and used to set the bin content for each channe:
//             entries - Root GetEntries() (Includes under and overflow. The following do not.)
//               count - Root Integral()
//                mean - Root GetMean()
//                peak - Root GetMaximumBin()
//                 rms - Root GetRMS()
//               rmsFF - Root GetRMS()
//              fitXXX - Parameter XXX from the fit made to the summary histogram, e.g. Mean for gaus.
//              fitratXXX - Ratio of parameter XXX from the fit to the mean from the fit.
//  errType - Specifies the metric used to set the bin error for each channel. Any of the value options or:
//                none - Do not set error
//                 rms - Root GetRMS()
//           meanError - Root GetMeanError()
//            rmsError - Root GetRMSError()
//               rmsFF - max(Root GetRMS(), FF)
//                zero - Set the error to zero
//     bins - if > 0,  plot # channels vs. variable in nbins bins
//            if 0, plot variable vs channel
//     pran - Range of y axis: ymin:ymax:yscal
//            yscal = pamp: Multiply range by pulserAmplitude
//            yscal = pampg14: Multiply range by pulserAmplitude*pulserGain/14.0
//     plot - Name of the file where the histogram should be plotted.
//            The histogram name is substituted for %HNAME%.
//       cr - Name of the channel range to plot. If "list", each value in ChannelRanges.
//
// Output data map for view:
//           int              roiRun - Run number
//           int           roiSubRun - SubRun number
//           int            roiEvent - Event number
//           int          roiChannel - Channel number
//           int         roiRawCount - # ROI (nROI) before selection
//           int            roiCount - # ROI (nROI)
//           int     roiNTickChannel - # ticks in the ADC channel
//     int[nROI]           roiTick0s - First tick for each ROI
//     int[nROI]           roiNTicks - # ticks for each ROI
//     int[nROI]      roiNUnderflows - # bins with ADC underflow in each ROI
//     int[nROI]       roiNOverflows - # bins with ADC overflow in each ROI
//     int[nROI]         roiTickMins - tick-tick0 for signal minimum for each ROI
//     int[nROI]         roiTickMaxs - tick-tick0 for signal maximum for each ROI
//   float[nROI]          roiSigMins - Signal minimum for each ROI
//   float[nROI]          roiSigMaxs - Signal maximum for each ROI
//   float[nROI]         roiSigAreas - Signal area for each ROI
//    TH1*[nROI]            roiHists - Histogram of samples or raw for each ROI
// If fit is done:
//   float[nROI]       roiFitHeights - Height in signal units from fit
//   float[nROI]        roiFitWidths - Shaping time in ticks from fit
//   float[nROI]     roiFitPositions - T0 in ticks from fit
//   float[nROI]    roiFitChiSquares - Chi-square from fit
//   float[nROI] roiFitChiSquareDofs - Chi-square/DOF from fit
//     int[nROI]         roiFitStats - Return status from fit
// If run data is found:
//         float          tpcAmpGain - Nominal preamp gain setting [mV/fC]
//         float      tpcShapingTime - Nominal preamp shaping time [us]
//           int         pulserIndex - Pulser gain setting (0, 1, ..., 63)
//           int      pulserInternal - 1 if pulser on ADC ASIC is used, 0 for FEMB pulser
// If fit, run data ane non-zero pulser index:
//  float(nRoi]
//
// Output data map for viewMap:
//           int        roiCount - # ROI (nROI)
//           int roiChannelCount - # channels
//     int roiFailedChannelCount - # failed channels
//   int[nfail]   failedChannels - List of failed channels
//
// Lines are draw at N*ChannelLineModulus + ChannelLinePattern[i] for any
// integer N and any value if i in range of the array which are within
// the drawn channel range.
// If ChannelLineModulus is zero, then lines are drawn for the channels in
// ChannelLinePattern.


#ifndef AdcRoiViewer_H
#define AdcRoiViewer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Data/RunData.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include <iostream>

class AdcChannelStringTool;
class RunDataTool;
class TimeOffsetTool;
class IndexRangeTool;
class TH1;

class AdcRoiViewer : TpcDataTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using HistVector = std::vector<TH1*>;
  using HistMap = std::map<Name, TH1*>;
  using HistVectorMap = std::map<Name, HistVector>;
  using NameMap = std::map<Name, Name>;
  using NameVectorMap = std::map<Name, NameVector>;
  using ChannelRange = IndexRange;
  using ChannelRangeMap = std::map<Name, ChannelRange>;
  using FloatMap = std::map<Name, float>;
  using IndexByIndexMap = std::map<Index, Index>;
  using IndexByNameMap = std::map<Name, Index>;
  using TpmPtr = std::unique_ptr<TPadManipulator>;
  using TpmMap = std::map<Index, TpmPtr>;
  using TpmNameMap = std::map<Index, Name>;
  using TpmCountMap = std::map<Index, Index>;

  // Subclass that associates a variable name with a histogram.
  //  vary != "" ==> 2D histo
  class HistInfo {
  public:
    TH1* ph = nullptr;
    Name varx;
    Name vary;
    Name fitName;
    Name plotName;
    float plotWidth;
  };

  using HistInfoMap = std::map<Name, HistInfo>;

  // This subclass carries the state for this tool, i.e. data that can change
  // after initialization.
  class State {
  public:
    // ROI plots.
    TpmMap roiPads;
    TpmNameMap roiPadNames;
    TpmCountMap roiPadCounts;
    // Summary histogram templates.
    HistInfoMap sumHistTemplates;
    // Summary histograms.
    HistMap sumHists;            // Histograms indexed by histogram name.
    HistVectorMap sumPlotHists;  // Histograms for each plot indexed by plot template name.
    NameMap sumFitNames;         // Fit name for each plotted histogram indexed by hist name.
    NameMap sumPlotNames;        // File names for each plotted histogram indexed by hist name.
                                 // The first file name is used for plots with multiple hists.
    FloatMap sumPlotWidths;      // Plot width for each plotted histogram indexed by hist name.
    IndexByNameMap sumHistChannels; // Channel for each sum histogram
    // Channel summary histograms.
    HistMap chanSumHists;
    NameMap chanSumHistTemplateNames;   // Sum template name indexed by chansum name
    NameMap chanSumHistVariableTypes;   // Variable type indexed by chansum name.
    NameMap chanSumHistErrorTypes;      // Error type indexed by chansum name.
    IndexByNameMap chanSumHistTypes;    // Type (0=var vs. chan, 1=#ROI vs var)
    NameMap chanSumPlotNames;           // Plot name indexed by chansum name
    FloatMap chanSumPlotYMins;          // Min value of y for plot.
    FloatMap chanSumPlotYMaxs;          // Max value of y for plot.
    NameMap chanSumPlotYOpts;           // Y scaling option for plot ("", "pamp")
    IndexByNameMap chanSumChaBegin;     // First channel for each histogram.
    IndexByNameMap chanSumChaEnd;       // Last channel for each histogram.
    ~State();
    // Fetch properties indexed by a histogram name.
    TH1* getSumHist(Name hnam);
    Name getSumFitName(Name hnam) const;
    Name getSumPlotName(Name hnam) const;
    float getSumPlotWidth(Name hnam) const;
    Name getChanSumHistTemplateName(Name hnam) const;
    Name getChanSumHistVariableType(Name hnam) const;
    Name getChanSumHistErrorType(Name hnam) const;
    Name getChanSumPlotName(Name hnam) const;
    Index getChannelStatus(Index icha) const;
    Index getChannelStatus(Name hnam) const;  // Argument is a chansum histogram name
    Index cachedRunCount = 0;  // Increment each time run number changes.
    Index cachedRun = AdcChannelData::badIndex();
    Name cachedSampleUnit;
    Index nRoiPlot =0;
    IndexByIndexMap channelStatuses;     // Status indexed by channel number
    // Run data.
    RunData runData;
    // Count doView calls.
    Index callCount =0;                   // Total
    IndexByIndexMap eventCallCount;       // For each event. Size of this is the event count.
    Index closeCount = 0;
  };

  using StatePtr = std::shared_ptr<State>;

  AdcRoiViewer(fhicl::ParameterSet const& ps);

  ~AdcRoiViewer() override;

  // AdcChannelTool methods.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }
  DataMap close(const DataMap* dmin) override;

  // Internal methods where most of the work is done.

  // Read ROIs, fit them and put results in data map.
  using DataMapVector = std::vector<DataMap>;
  int doView(const AdcChannelData& acd, int dbg, DataMap& dm) const;

  // Save the ROI histograms to the ROI Root file.
  void writeRoiHists(const DataMap& res, int dbg) const;
  void writeRoiHists(const DataMapVector& res, int dbg) const;

  // Plot the ROIs for an event/channel.
  // ADC channel data is used to build plot names.
  void writeRoiPlots(const HistVector& hists, const AdcChannelData& acd) const;

  // Return the state.
  State& getState() const { return *m_state; }

  // Fill the summary histograms for one channel.
  void fillSumHists(const AdcChannelData& acd, const DataMap& dm) const;

  // Fit the summary histograms to the summary Root file.
  void fitSumHists() const;

  // Fill the channel summary histograms.
  void fillChanSumHists() const;

  // Write the summary histograms to the summary Root files
  // and summary plots to the plot files.
  void writeSumHists() const;
  void writeSumPlots(const DataMap* pdmin) const;
  void writeChanSumHists() const;
  void writeChanSumPlots() const;

  // Replace %LABX% with the corresponding plot label.
  void setPlotLabels(Name& sttl) const;

private:

  // Configuration data.
  int m_LogLevel;
  float m_SigThresh;
  Index m_TickBorder;
  int m_RoiHistOpt;
  int m_FitOpt;
  Index m_RoiPlotOpt;
  int m_MaxRoiPlots;
  Index m_RoiPlotPadX;
  Index m_RoiPlotPadY;
  time_t m_StartTime;
  float m_PulserStepCharge;
  float m_PulserDacOffset;
  Name m_PulserChargeUnit;
  bool m_SumNegate;
  Index m_SumPlotPadX;
  Index m_SumPlotPadY;
  Index          m_ChannelLineModulus;
  IndexVector    m_ChannelLinePattern;
  Name m_RunDataTool;
  Name m_TickOffsetTool;
  Name m_ChannelRangeTool ="channelRanges";
  Name m_RoiRootFileName;
  Name m_SumRootFileName;
  Name m_ChanSumRootFileName;
  NameVector m_ChannelRanges;
  NameVector m_PlotLabels;

  // Derived from configuration.
  ChannelRangeMap m_crmap; 
  NameMap m_plotLabelSubs;

  // Shared pointer so we can make sure only one reference is out at a time.
  StatePtr m_state;

  // Tools.
  const AdcChannelStringTool* m_adcStringBuilder  =nullptr;
  const RunDataTool*          m_pRunDataTool      =nullptr;
  const TimeOffsetTool*       m_pTickOffsetTool   =nullptr;
  const IndexRangeTool*       m_pChannelRangeTool =nullptr;

};


#endif
