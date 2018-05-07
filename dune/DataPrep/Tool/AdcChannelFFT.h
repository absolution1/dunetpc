// AdcChannelFFT.h

// David Adams
// April 2018
//
// Tool to perform FFT of the prepared data in an ADC channel data map.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   FirstTick - First tick used in transform.
//   NTick - # ticks to use in transform. To end if 0.
//   NormOpt - Normalization:
//     0 - Root default: Transform followed by inverse will be N times larger than original
//     1 - Fixed power: Transform followed by inverse gives back original
//                      Normalization of 1/sqrt(N) applied at both stages.
//                      Power for original and transform are the same
//     2 - Compensating: Normalization is 1/N.
//                       Transform with option 0 followed by inverse with option 2 or
//                       Transform with option 2 followed by inverse with option 0
//                       results in the original samples.

#ifndef AdcChannelFFT_H
#define AdcChannelFFT_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AdcChannelFFT : AdcChannelTool {

public:

  using Index = unsigned int;

  AdcChannelFFT(fhicl::ParameterSet const& ps);

  ~AdcChannelFFT() override =default;

  DataMap view(const AdcChannelData& acd) const override;
  bool updateWithView() const override { return true; }

private:

  // Configuration data.
  int    m_LogLevel;
  Index  m_FirstTick;
  Index  m_NTick;
  Index  m_NormOpt;

};

DEFINE_ART_CLASS_TOOL(AdcChannelFFT)

#endif
