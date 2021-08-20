// AdcRoiSlicer.h

// David Adams
// September 2019
//
// Tool to slice ADC data on ROI boundaries.
//
// Configuration parameters.
//         LogLevel - Logging level (0=none, 1=ctor only, ...)
//      OutViewName - Name for the view where the slices are recorded.
//         SliceOpt - 1 = ROIs, 2 = not ROIs, 3 = both
//          CopyRaw - If true, raw data is copied to the new view

#ifndef AdcRoiSlicer_H
#define AdcRoiSlicer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"

class AdcRoiSlicer
: public TpcDataTool {

public:

  using Name = std::string;

  AdcRoiSlicer(fhicl::ParameterSet const& ps);

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int m_LogLevel;
  Name m_OutViewName;
  int m_SliceOpt;
  bool m_CopyRaw;

};


#endif
