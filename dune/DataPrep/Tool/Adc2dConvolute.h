// Adc2dConvolute.h
//
// Tool to perform 2D (time and channel) convolution of AdcChannelData.
//
// The transformation is
//
//   d(m;i) = SUM(n;j) R(m,n;i,j) q(n,j)
//
// where q denotes the input charge with NT (BinsPerTick) bins in time (i, j) and
// NW (BinsPerWire) bins in channel (m, n).
//
// The response matrix is constructed assuming translational invariance (except for edges)
// in both time and channel, and reflection symmetry in channel, and so takes the form
//
//   R(m,n;i,j) = r(S, T)
//              = V(A, B, T)
//   T = NT*i - (j - Joff)
//   S = 0, 1, 2, ...
// 
// where
//
//   A = abs(n/NW - m)
//   B = n % NW,            for n - NW m < N/2 or
//       NW - 1 - n % NW    otherwise
//
// The response is provided as a vector r (ResponseVectors) which is a vector over channel
// bins of vectors over time bins with
//   r = { r(0, (NW-1)/2), ..., r(0, NW-1)],
//         r(1, 0), ..., r(1, NW-1),
//         r(2, 0), ..., r(2, NW-1),
//         ... }
//
// where v(A,B) is the response of a channel to unit charge in bin B of the channel R to
// the right, i.e. the second half of the diagonal reponse vector followed by the vectors
// or the right neighbors. Left elements are obtained by symmetry as shown above.
//
// The input charge is taken from binSamples or samples according to BinsPerTick.
// Output charge is stored in samples.
// If the reponse extends over multiple channels, then contiguous numbering from MinChannel to
// MaxChannel, inclusive is used.
//
// Configuration:
//   LogLevel: 0=silent, 1=init, 2=each event, >2=more
//   BinsPerTick - # input bins for each tick. Value 0 means 1 with input from samples.
//   BinsPerWire - # bins for each input channel
//   ResponseVectors: Vector of discrete sequences representing the response function.
//   ReponseCenter: Center tick for the response function (Joff above).
//   MinChannel: First channel
//   MaxChannel: Last channel
//

#ifndef Adc2dConvolute_H
#define Adc2dConvolute_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"

class Adc2dConvolute : TpcDataTool {

public:

  using Index = unsigned int;

  Adc2dConvolute(fhicl::ParameterSet const& ps);

  ~Adc2dConvolute() override =default;

  DataMap updateMap(AdcChannelDataMap& acd) const override;

private:

  using Name = std::string;
  using IndexVector = std::vector<Index>;
  using ResponseVector = AdcSignalVector;   // Same type as AdcChannel data
  using ResponseVectorVector = std::vector<ResponseVector>;

  // Configuration data.
  int                  m_LogLevel;
  Index                m_BinsPerTick;
  Index                m_BinsPerWire;
  ResponseVectorVector m_ResponseVectors;
  Index                m_ResponseCenter;
  Index                m_MinChannel;
  Index                m_MaxChannel;

};


#endif
