// Roi2dToAdc.h
//
// David Adams
// February 2021
//
// Tool to copy 2D ROIs to maps of ADC channel data.
//
// The samples in the ADC channel data are zeroed and then replaced
// with ROI data where present.Updated samples are flagges as signal.
//
// Configuration:
//            LogLevel - Logging level: 0=none, 1=init, 2=call, ...
//
// Returned data map:
//   r2a_nroi - # ROIS copied                        
//   r2a_nchaZeroed - # channels zeroed in ADC data                
//   r2a_nchaFilled - # channels filled in ADC data                
//   r2a_nsamZeroed - # samples zeroed in ADC data                
//   r2a_nsamFilled - # samples filled in ADC data                
//   r2a_nsig - # samples flags as signal             

#ifndef Roi2dToAdc_H
#define Roi2dToAdc_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneInterface/Data/IndexRange.h"
#include "dune/DuneInterface/Data/RunData.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include <iostream>

class Roi2dToAdc : TpcDataTool {

public:

  using Index = unsigned int;

  Roi2dToAdc(fhicl::ParameterSet const& ps);

  ~Roi2dToAdc() override;

  // TpcDataTool methods.
  DataMap updateTpcData(TpcData& tpd) const override;

private:

  // Configuration data.
  int m_LogLevel;

};


#endif
