// AdcRangeSampleScaler.h
//
// David Adams
// Augest 2020
//
// Tool to scale the samples in AdcSampleData according to channel range.
//
// Configuration:
//   LogLevel: 0=silent, 1=init, 2=each event, >2=more
//   RangeLimits: Channel limits
//   ScaleFactors: Factor by which samples are scaled.
//
// All samples for channel ICHA are scaled by ScaleFactors[IRAN] where
// IRAN is the first range for which ICHA < Rangelimits[IRAN].
// If RangeModulus > 0, then ICHA%RangeModulus is used in place of ICHA.

#ifndef AdcRangeSampleScaler_H
#define AdcRangeSampleScaler_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"

class AdcRangeSampleScaler : TpcDataTool {

public:

  AdcRangeSampleScaler(fhicl::ParameterSet const& ps);

  ~AdcRangeSampleScaler() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using FloatVector = std::vector<float>;

  // Configuration data.
  int                m_LogLevel;
  IndexVector        m_RangeLimits;
  Index              m_RangeModulus;
  FloatVector        m_ScaleFactors;

};


#endif
