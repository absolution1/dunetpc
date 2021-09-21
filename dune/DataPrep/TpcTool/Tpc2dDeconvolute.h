// Tpc2dDeconvolute.h
//
// David Adams
// August 2021
//
// Tool to perform 2D (time and channel) deconvolution of the 2D ROIs in TpcData.
//
// The observed charge qd is assumed to be related to the true charge by
//   qd(m,i) = SUM(n,j) r(m,n;i,j) q0(n,j)
// where the indices repectively denote channel and time bin (tick) and r
// is the unbinned response matrix.
//
// Discrete Fourier transforms (DFTs) are used to express this in channel (kn, km)
// and time (ki, kj) wave numbers:
//
//   Qd(km,ki) = R(km-kn,ki-kj) Q0(kn,kj)
//
// with charge and response assumed periodic in channel and time.
// This expression is then easily inverted to obtain and estimate Qe of Q0:
//
//   Qe(kn,kj) = Fc(kn) Ft(kj) Qd(km,ki) /R(km-kn,ki-kj)
//
// HF filter functions Fc and Ft are introduced to suppress large values that can
// arise where R approaches zero and noise in the measurement of Qd keeps
// that value finite.
// The filter functions are Gaussian derived from user-supplied Gaussian sigmas in
// sample (tick) and channel.
//
// Low-frequency filtering is provided in sample only by the function
//   1/[ 1 + 1/(W N k)^P ]
// N is the number of samples and and k is the wave number.
// The width W and power P are configuration parameters.
//   W < 0 disables the filter
//   W = 0 filters out k = 0 only
//   W > 0 removes k = 0 and applies the above for other values
//
// Other options may later be added for the response or filter.
//
// The response matrix is derived from an input reponse matrix ri(m,i) read from a
// channel-indexed vector of tick-indexed reponse vectors. This input reponse is
// the signal induced in the mth neighbor channel by a charge deposited at tick0.
// The tick-indexed reponse vectors can have varying lenght.
//
// The channel-tick dimensions of the response matrix r used in the deconvolution are
// those of the input data which must be rectangular in that space.
// The reponse matrix is derived from the input reponse assuming reflection symmetry
// in channel. Missing entries are filled with zeroes and those byond the data
// dimensions are ignored.
//
// Configuration:
//   LogLevel: 0=silent, 1=init, 2=each event, >2=more
//   ResponseVectors: Vector of discrete sequences giving the 2D response (ri above).
//   ReponseCenter: Center tick for the response function (tick0 above).
//   FftSize: Maximum FFT size, roughly ncha*nsam.
//   InPath: Path to the input 2D ROIs.
//   OutPath: Path to the output 2D ROIs. If "" or ".", the input ROIs are updated.
//   SampleSigma: Tick sigma used to construct the sample filter. Zero means no filter.
//   ChannelSigma: Channel sigma used to construct the sample filter. Zero means no filter.
//   LowFilterWidth: - Width [tick] for the low-frequency sample filter. Zero means only k=0.
//   LowFilterPower: - Power for the LF sample filter. Value <=0 disables the filter.
//

#ifndef Tpc2dDeconvolute_H
#define Tpc2dDeconvolute_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "dune/DuneCommon/Utility/Fw2dFFT.h"

class Tpc2dDeconvolute : TpcDataTool {

public:

  using Index = unsigned int;

  Tpc2dDeconvolute(fhicl::ParameterSet const& ps);

  ~Tpc2dDeconvolute() override =default;

  // Deconvolute.
  DataMap updateTpcData(TpcData&) const override;

  // Viewing fails (for now).
  DataMap viewTpcData(const TpcData&) const override;

  // Return the FFT transform (eventually for this thread).
  Fw2dFFT& fft() const;

private:

  using Name = std::string;
  using IndexVector = std::vector<Index>;
  using ResponseVector = AdcSignalVector;   // Same type as AdcChannel data
  using ResponseVectorVector = std::vector<ResponseVector>;

  // Configuration data.
  int                  m_LogLevel;
  ResponseVectorVector m_ResponseVectors;
  int                  m_ResponseCenter;
  Index                m_FftSize;
  Name                 m_InPath;
  Name                 m_OutPath;
  float                m_SampleSigma;
  float                m_ChannelSigma;
  float                m_LowFilterWidth;
  float                m_LowFilterPower;

  // Transform.
  std::unique_ptr<Fw2dFFT> m_pfft;

};

DEFINE_ART_CLASS_TOOL(Tpc2dDeconvolute)

#endif
