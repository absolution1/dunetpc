// PddpVintageDeconvoluter.h
//
// Tool to preform deconvolution of data from an ADC.
// This is just a copy of VintangeDeconvoluter tool 
// but with ProtoDUNE DP signal shaping service. 
//
// This code should be deleted point when different signal 
// shaping services could be loaded via a single interface
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more

#ifndef PddpVintageDeconvoluter_H
#define PddpVintageDeconvoluter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <vector>

class PddpVintageDeconvoluter : TpcDataTool {

public:

  PddpVintageDeconvoluter(fhicl::ParameterSet const& ps);

  ~PddpVintageDeconvoluter() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int            m_LogLevel;

};


#endif
