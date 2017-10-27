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
//
// Output data map:
//           int        roiCount - # ROI (nROI)
//           int roiNTickChannel - # ticks in the ADC channel
//     int[nROI]       roiTick0s - First tick for each ROI
//     int[nROI]       roiNTicks - # ticks for each ROI
//     int[nROI]  roiNUnderflows - # bins with ADC underflow in each ROI
//     int[nROI]   roiNOverflows - # bins with ADC overflow in each ROI
//   float[nROI]      roiSigMins - Signal minimum for each ROI
//   float[nROI]      roiSigMaxs - Signal minimum for each ROI
//   float[nROI]     roiSigAreas - Signal area for each ROI
//    TH1*[nROI]        roiHists - Histogram of samples or raw for each ROI

#ifndef AdcRoiViewer_H
#define AdcRoiViewer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelViewer.h"
#include <iostream>

class AdcRoiViewer : AdcChannelViewer {

public:

  AdcRoiViewer(fhicl::ParameterSet const& ps);

  ~AdcRoiViewer() override =default;

  DataMap view(const AdcChannelData& acd) const override;

private:

  // Configuration data.
  int m_LogLevel;
  int m_HistOpt;

};

DEFINE_ART_CLASS_TOOL(AdcRoiViewer)

#endif
