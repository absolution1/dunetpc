// AdcRoiViewer.h

// David Adams
// October 2017
//
// Tool to extract information about the ROIs in an ADC channel.
//
// Configuration:
//   LogLevel - Login level: 0=none, 1=init, 2=call, ...
//
// Output data map:
//           int     nROI - # ROI
//    TH1*[nROI] roiHists - Histogram of samples for each ROI
//   float[nROI] roiSigMins - Signal minimum for each ROI

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

};

DEFINE_ART_CLASS_TOOL(AdcRoiViewer)

#endif
