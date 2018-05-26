// AdcRoiViewer.h

// David Adams
// October 2017
//
// Tool to extract information about the ROIs in an ADC channel.
//
// Configuration:
//   LogLevel - Login level: 0=none, 1=init, 2=call, ...
//   HistOpt - histo option:  0 - No histograms
//                            1 - sample vs. tick
//                            2 - raw vs. tick
//                           10 + i - As above for i except vs. tick - tick0
//           FitOpt - Fit option
//                      0 - no fit
//                      1 - fit with coldelecReponse
//  RoiRootFileName - Name of file to which the ROI histograms are copied.
//  SumRootFileName - Name of file to which the evaluated parameter histograms are copied.
//
// Output data map for view:
//           int          roiRun - Run number
//           int       roiSubRun - SubRun number
//           int        roiEvent - Event number
//           int      roiChannel - Channel number
//           int        roiCount - # ROI (nROI)
//           int roiNTickChannel - # ticks in the ADC channel
//     int[nROI]       roiTick0s - First tick for each ROI
//     int[nROI]       roiNTicks - # ticks for each ROI
//     int[nROI]  roiNUnderflows - # bins with ADC underflow in each ROI
//     int[nROI]   roiNOverflows - # bins with ADC overflow in each ROI
//     int[nROI]     roiTickMins - tick-tick0 for signal minimum for each ROI
//     int[nROI]     roiTickMaxs - tick-tick0 for signal maximum for each ROI
//   float[nROI]      roiSigMins - Signal minimum for each ROI
//   float[nROI]      roiSigMaxs - Signal maximum for each ROI
//   float[nROI]     roiSigAreas - Signal area for each ROI
//    TH1*[nROI]        roiHists - Histogram of samples or raw for each ROI
// If fit is done:
//   float[nROI]   roiFitHeights - Height in signal units from fit
//   float[nROI]    roiFitWidths - Shaping time in ticks from fit
//   float[nROI] roiFitPositions - T0 in ticks from fit
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
using HistVector = std::vector<TH1*>;
using HistMap = std::map<std::string, TH1*>;

class AdcRoiViewer : AdcChannelTool {

public:

  using Index = unsigned int;
  using Name = std::string;

  class State {
  public:
    HistMap hists;
    ~State();
    // Fetch the histogram for a histogram name.
    TH1* getHist(Name hname);
  };

  using StatePtr = std::shared_ptr<State>;

  AdcRoiViewer(fhicl::ParameterSet const& ps);

  ~AdcRoiViewer() override =default;

  // AdcChannelTool methods.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

  // Internal methods where most of the work is done.
  using DataMapVector = std::vector<DataMap>;
  int doView(const AdcChannelData& acd, int dbg, DataMap& dm) const;
  void doSave(const DataMap& res, int dbg) const;
  void doSave(const DataMapVector& res, int dbg) const;

  // Return the state.
  StatePtr getState() const { return m_state; }

  // Name for dist var histogram from variable name and channel data.
  Name getDistHistName(Name vname, const AdcChannelData& acd) const;

  // Write the dist hists.
  void writeDistHists(Index dbg) const;

private:

  // Configuration data.
  int m_LogLevel;
  int m_HistOpt;
  int m_FitOpt;
  std::string m_RoiRootFileName;
  std::string m_SumRootFileName;

  // Shared pointer so we can make sure only one reference is out at a time.
  StatePtr m_state;

  // ADC string tools.
  const AdcChannelStringTool* m_adcStringBuilder;

};

DEFINE_ART_CLASS_TOOL(AdcRoiViewer)

#endif
