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
//   FitOpt - Fit option
//              0 - no fit
//              1 - fit with coldelecReponse
//   RootFileName - Name of file to which histograms are copied.
//
// Output data map for view:
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

class AdcRoiViewer : AdcChannelTool {

public:

  AdcRoiViewer(fhicl::ParameterSet const& ps);

  ~AdcRoiViewer() override =default;

  // AdcChannelTool methods.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

  // Internal methods where most of the work is done.
  using DataMapVector = std::vector<DataMap>;
  int doView(const AdcChannelData& acd, int dbg, DataMap& dm) const;
  void doSave(const DataMap& res, std::string ofrname, int dbg) const;
  void doSave(const DataMapVector& res, std::string ofrname, int dbg) const;

private:

  // Configuration data.
  int m_LogLevel;
  int m_HistOpt;
  int m_FitOpt;
  std::string m_RootFileName;

  // ADC string tools.
  const AdcChannelStringTool* m_adcStringBuilder;

};

DEFINE_ART_CLASS_TOOL(AdcRoiViewer)

#endif
