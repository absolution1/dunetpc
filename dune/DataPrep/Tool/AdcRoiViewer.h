// AdcRoiViewer.h

// David Adams
// October 2017
//
// Tool to extract information about the ROIs in an ADC channel.
//
// Configuration:
//         LogLevel - Logging level: 0=none, 1=init, 2=call, ...
//     SigMinThresh - if <0, then only keep ROIS with a tick below this value
//     SigMaxThresh - if >0, then only keep ROIS with a tick above this value
//       RoiHistOpt - histo option:  0 - No histograms
//                                   1 - sample vs. tick
//                                   2 - raw vs. tick
//                                  10 + i - As above for i except vs. tick - tick0
//           FitOpt - ROI fitting option
//                      0 - no fit
//                      1 - fit with coldelecReponse
//         SumHists - Array of summary histogram specifiers. See below.
//  RoiRootFileName - Name of file to which the ROI histograms are copied.
//  SumRootFileName - Name of file to which the evaluated parameter histograms are copied.
//
// A summary histogram specifier is a paramter set with the following fields:
//     var: Name of the variable to draw:
//            fitHeight
//            fitWidth
//            fitPosition
//            fitTickRem - Fractional part of tick
//            fitStat
//            fitChiSquare
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
// Output data map for view:
//           int           roiRun - Run number
//           int        roiSubRun - SubRun number
//           int         roiEvent - Event number
//           int       roiChannel - Channel number
//           int      roiRawCount - # ROI (nROI) before selection
//           int         roiCount - # ROI (nROI)
//           int  roiNTickChannel - # ticks in the ADC channel
//     int[nROI]        roiTick0s - First tick for each ROI
//     int[nROI]        roiNTicks - # ticks for each ROI
//     int[nROI]   roiNUnderflows - # bins with ADC underflow in each ROI
//     int[nROI]    roiNOverflows - # bins with ADC overflow in each ROI
//     int[nROI]      roiTickMins - tick-tick0 for signal minimum for each ROI
//     int[nROI]      roiTickMaxs - tick-tick0 for signal maximum for each ROI
//   float[nROI]       roiSigMins - Signal minimum for each ROI
//   float[nROI]       roiSigMaxs - Signal maximum for each ROI
//   float[nROI]      roiSigAreas - Signal area for each ROI
//    TH1*[nROI]         roiHists - Histogram of samples or raw for each ROI
// If fit is done:
//   float[nROI]    roiFitHeights - Height in signal units from fit
//   float[nROI]     roiFitWidths - Shaping time in ticks from fit
//   float[nROI]  roiFitPositions - T0 in ticks from fit
//   float[nROI] roiFitChiSquares - Chi-square from fit
//     int[nROI]      roiFitStats - Return status from fit
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
#include <iostream>

class AdcChannelStringTool;

class AdcRoiViewer : AdcChannelTool {

public:

  using Index = unsigned int;
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using HistVector = std::vector<TH1*>;
  using HistMap = std::map<Name, TH1*>;
  using NameMap = std::map<Name, Name>;

  // This subclass carries the state for this tool, i.e. data that can change
  // after initialization.
  class State {
  public:
    // Variable for each histogram.
    NameVector vars;
    // Fit function for each histogram.
    NameVector fits;
    // Template for each histogram.
    HistVector histTemplates;
    // Filled histograms.
    HistMap hists;
    ~State();
    // Fetch the histogram for a histogram name.
    TH1* getHist(Name hname);
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

  // Write the summary histograms to the summary Root file.
  void writeSumHists() const;

private:

  // Configuration data.
  int m_LogLevel;
  float m_SigMinThresh;
  float m_SigMaxThresh;
  int m_RoiHistOpt;
  int m_FitOpt;
  std::string m_RoiRootFileName;
  std::string m_SumRootFileName;

  // Fixed configuration data.
  int m_TickPeriod = 497;

  // Shared pointer so we can make sure only one reference is out at a time.
  StatePtr m_state;

  // ADC string tools.
  const AdcChannelStringTool* m_adcStringBuilder;

};

DEFINE_ART_CLASS_TOOL(AdcRoiViewer)

#endif
