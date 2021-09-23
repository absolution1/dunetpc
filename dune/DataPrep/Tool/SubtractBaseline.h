////////////////////////////////////////////////////////////////////////
// SubtractBaseline.h
//
// Tool to preform subtract baseline using linear interpolation between 
// regions defined by the datasize and fBaseSampleBins
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   BaseSampleBins - 
//   BaseVarCut -
//
/////////////////////////////////////////////////////////////////////////
#ifndef SubtractBaseline_H
#define SubtractBaseline_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <vector>

class SubtractBaseline : TpcDataTool {

public:

  SubtractBaseline(fhicl::ParameterSet const& ps);

  ~SubtractBaseline() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int            m_LogLevel;
  int            m_BaseSampleBins;
  float          m_BaseVarCut;

};


#endif
