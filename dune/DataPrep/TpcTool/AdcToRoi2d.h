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
//                         2 - All channel maps merged into a single ROI.
//        InputAdcMaps - Indices of the ADC maps that should be processed.
//                       If empty, all maps are processed.
//         OutputNames - Names for the subdirectories where the 2D ROIs are written.
//                       If blank, they are written into the top level.
//                       Otherwise there should be one name for each entry in the ADC vector.
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
  using IndexVector = std::vector<Index>;
  using Name = std::string;
  using NameVector = std::vector<Name>;

  AdcToRoi2d(fhicl::ParameterSet const& ps);

  ~AdcToRoi2d() override;

  // TpcDataTool methods.
  DataMap updateTpcData(TpcData& tpd) const override;

private:

  // Configuration data.
  int m_LogLevel;
  Index m_Option;
  IndexVector m_InputAdcMaps;
  NameVector m_OutputNames;

};


#endif
