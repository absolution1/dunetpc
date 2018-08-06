// AdcRoiViewer.h
//
// David Adams
// October 2017
//
// Tool to extract information about the ROIs in an ADC channel.
//
// Configuration:
//           LogLevel - Logging level: 0=none, 1=init, 2=call, ...
//          SigThresh - if <0, then only keep ROIs with a tick below this value
//                      if >0, then only keep ROIs with a tick above this value
//         TickBorder - if > 0, only keep ROIs forwhich there are this many ticks or more
//                      before the start and after the end of the ROI.
//         RoiHistOpt - histo option:  0 - No histograms
//                                     1 - sample vs. tick
//                                     2 - raw vs. tick
//                                    10 + i - As above for i except vs. tick - tick0
//             FitOpt - ROI fitting option
//                        0 - no fit
//                        1 - fit with coldelecReponse
//   PulserStepCharge - Charge per unit step in a pulser run
//    PulserDacOffset - Offset in pulser: Qin = PulserStepCharge*(DAC - PulserDacOffset)
//   PulserChargeUnit - Unit for the pulser charge (ke, fC, ...)
//           SumHists - Array of summary histogram specifiers. See below.
//      ChannelRanges - Ranges of channels for channel summary plots.
//       ChanSumHists - Array of specifiers for the channel summary histograms.
//        RunDataTool - Name for the run data tool. If found and pulser is on, then each
//                      ROI is assigned a charge corresponding to the pulser setting.
//     TickOffsetTool - Name of the tool that provides the tick offset.
//    RoiRootFileName - Name of file to which the ROI histograms are copied.
//    SumRootFileName - Name of file to which the evaluated parameter histograms are copied.
//
// Summary histograms
// ------------------
// Summary histograms show the number of ROIs is bins of some ROI variable such
// as tick or pulse height.
// A summary histogram specifier is a parameter set with the following fields:
//     var: Name of the variable to draw:
//            fitHeight - Height from ROI fit
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
//    name: Name of the histogram. Include %CHAN% to get separate histos for each channel
//   title: Histogram title
//    nbin: # bins
//    xmin: Lower edge of the first bin
//    xmax: Upper edge of the last bin
//     fit: Name of function used to fit the distribution, e.g. "gaus"
// If xmin < xmax and xmin > 0, the range will have width xmin centered on the median
// value for the first event. If xmax > 0, the lower edge is rounded to that value.
// If xmin <= xmax otherwise (e.g. xmin = xmax = 0), Root will do autoscaling of the axis.
// E.g.: {var:fitHeight name:"hfh_chan%0CHAN" title:"Fit height for channel %CHAN%"
//        nbin:50 xmin:0.0 xmax:5.0}
//
// Channel summary histograms
// --------------------------
// Channel summary histograms hold a metric for each channel derived from a summary histogram.
// The following fileds specify a channel summary histogram:
//     name - name for the summary histogram (%CRNAME% is replaced with the channel range name)
//    title - histogram title  (substitutions for %CRLABEL%, %RUN%, ...)
//  valHist - Name of the summary histogram template from which the metric is derived (should include %CHAN%)
//  valType - Specifies the metric to be extracted and used to set the bin content for each channe:
//                mean - Root GetMean()
//                 rms - Root GetRMS()
//              fitXXX - Parameter XXX from the fit made to the summary histogram, e.g. Mean for gaus.
//  errType - Specifies the metric used to set the bin error for each channel. Any of the value options or:
//                none - Do not set error
//                zero - Set the error to zero
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

#ifndef AdcRoiViewer_H
#define AdcRoiViewer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include <iostream>

class AdcChannelStringTool;
class RunDataTool;
class TimeOffsetTool;

class AdcRoiViewer : AdcChannelTool {

public:

  using Index = unsigned int;
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using HistVector = std::vector<TH1*>;
  using HistMap = std::map<Name, TH1*>;
  using NameMap = std::map<Name, Name>;
  using ChannelRange = IndexRange;
  using ChannelRangeMap = std::map<Name, ChannelRange>;

  // Subclass that associates a variable name with a histogram.
  //  vary != "" ==> 2D histo
  class HistInfo {
  public:
    TH1* ph = nullptr;
    Name varx;
    Name vary;
  };

  using HistInfoMap = std::map<Name, HistInfo>;

  // This subclass carries the state for this tool, i.e. data that can change
  // after initialization.
  class State {
  public:
    // Summary histogram templates.
    HistInfoMap sumHistTemplates;
    // Summary histograms.
    HistMap sumHists;
    // Channel summary histograms.
    HistMap chanSumHists;
    NameMap chanSumHistTemplateNames;  // Sum template name indexed by chansum name
    NameMap chanSumHistVariableTypes;  // Variable type indexed by chansum name.
    NameMap chanSumHistErrorTypes;     // Error type indexed by chansum name.
    ~State();
    // Fetch the summary histogram for a histogram name.
    TH1* getSumHist(Name hname);
    Name getChanSumHistTemplateName(Name hnam) const;
    Name getChanSumHistVariableType(Name hnam) const;
    Name getChanSumHistErrorType(Name hnam) const;
    Index cachedRunCount = 0;  // Increment each time run number changes.
    Index cachedRun = AdcChannelData::badIndex;
    Name cachedSampleUnit;
  };

  using StatePtr = std::shared_ptr<State>;

  AdcRoiViewer(fhicl::ParameterSet const& ps);

  ~AdcRoiViewer() override;

  // AdcChannelTool methods.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

  // Internal methods where most of the work is done.

  // Read ROIs, fit them and put results in data map.
  using DataMapVector = std::vector<DataMap>;
  int doView(const AdcChannelData& acd, int dbg, DataMap& dm) const;

  // Save the ROI histograms to the ROI Root file.
  void writeRoiHists(const DataMap& res, int dbg) const;
  void writeRoiHists(const DataMapVector& res, int dbg) const;

  // Return the state.
  State& getState() const { return *m_state; }

  // Fill the summary histograms for one channel.
  void fillSumHists(const AdcChannelData acd, const DataMap& dm) const;

  // Fit the summary histograms to the summary Root file.
  void fitSumHists() const;

  // Fill the channel summary histograms.
  void fillChanSumHists() const;

  // Write the summary histograms to the summary Root file.
  void writeSumHists() const;
  void writeChanSumHists() const;

private:

  // Configuration data.
  int m_LogLevel;
  float m_SigThresh;
  Index m_TickBorder;
  int m_RoiHistOpt;
  int m_FitOpt;
  float m_PulserStepCharge;
  float m_PulserDacOffset;
  Name m_PulserChargeUnit;
  Name m_RunDataTool;
  Name m_TickOffsetTool;
  Name m_RoiRootFileName;
  Name m_SumRootFileName;
  Name m_ChanSumRootFileName;
  ChannelRangeMap m_ChannelRanges;

  // Derived from configuration.

  // Shared pointer so we can make sure only one reference is out at a time.
  StatePtr m_state;

  // Tools.
  const AdcChannelStringTool* m_adcStringBuilder =nullptr;
  const RunDataTool*          m_pRunDataTool     =nullptr;
  const TimeOffsetTool*       m_pTickOffsetTool  =nullptr;

};

DEFINE_ART_CLASS_TOOL(AdcRoiViewer)

#endif
