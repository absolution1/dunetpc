// AdcToRoi2d.h
//
// David Adams
// February 2021
//
// Tool to convert maps of ADC channel data into 2D ROIs.
//
// Configuration:
//            LogLevel - Logging level: 0=none, 1=init, 2=call, ...
//              Option - Action to take:
//                         0 - none except logging
//                         1 - One ROI for each channel map
//                         2 - All channel maps merged into a singe ROI.
// Returned data map:
//   a2r_nroi - # ROIS created                        
//   a2r_nsams - # samples for each ROI                

#ifndef AdcToRoi2d_H
#define AdcToRoi2d_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Data/RunData.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include <iostream>

class AdcToRoi2d : TpcDataTool {

public:

  using Index = unsigned int;

  AdcToRoi2d(fhicl::ParameterSet const& ps);

  ~AdcToRoi2d() override;

  // TpcDataTool methods.
  DataMap updateTpcData(TpcData& tpd) const override;

private:

  // Configuration data.
  int m_LogLevel;
  Index m_Option;

};


#endif
