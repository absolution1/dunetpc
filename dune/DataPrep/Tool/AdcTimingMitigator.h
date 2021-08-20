// AdcTimingMitigator.h

// David Adams
// October 2017
//
// Tool to mitigate the timing in ADC samples to chave the sampling
// frequency and offset.
//
// If i is the old tick index and j is the new value:
//   t = j T2 = (i + ioff) T1
//   i = (T2/T1) j - ioff
//     = SamplingRatio j + SamplingOffset
//
// The value for each new sample is obtained by interpolating linearly
// between the surrounding old samples. The new waveform ends when
// the end of the old waveform is reached.
//
// If the offset is positive, the original waveform is missing the data
// required for interpolation at the start of the new waveform and the
// first value in the old waveform is used.
//
// Configuration parameters.
//         LogLevel - Logging level (0=none, 1=ctor only, ...)
//    SamplingRatio - Sampling period is increased by this factor.
//   SamplingOffset - Wave form is offset by this many (original) ticks
//            FEMBs - FEMB indices for which correction should be applied.

#ifndef AdcTimingMitigator_H
#define AdcTimingMitigator_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <string>
#include <vector>
#include <set>

class HistogramManager;
class TH1;

class AdcTimingMitigator
: public TpcDataTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;

  AdcTimingMitigator(fhicl::ParameterSet const& ps);

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int m_LogLevel;
  float m_SamplingRatio;
  float m_SamplingOffset;
  IndexVector m_FEMBs;

  IndexSet m_fembSet;

};


#endif
