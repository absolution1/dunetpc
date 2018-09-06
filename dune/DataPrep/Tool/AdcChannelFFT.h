// AdcChannelFFT.h

// David Adams
// April 2018
//
// Tool to perform FFT of the prepared data in an ADC channel data map.
//
// For AdcChannelData acd, the waveform (WF) is held in acd.samples, the DFT magnitudes
// in acd.dftmags and DFT phases in acd.fftphases. For a WF of length N, the length of the
// magnitudes vector is (N+2)/2 and the length of the phases is (N+1)/2, i.e. the same for
// N odd and one less for N even. I.e. N = N_mag + N_phase - 1. The first entry in the phase
// vector is always zero or pi. The magnitudes are all positive except, when N is even, there
// the last (Nyquist) magnitude is signed because there is no corresponding phase.
//
// The DFT magnitudes are stored in a float vector of lenght NMag normalized
// so that the relative power for each term is the square of its value.
// I.e. the DFT magnitudes only extend up to the Nyquist frequency and are normalized to
// include the power for all alias (e.g. folded) frequencies. To convert to a complex array
// of nsam terms all but the first and Nyquiset (last only for even # samples) must be scaled
// by 1/sqrt(2).
//
// The DFT phases
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   FirstTick - First tick used in transform.
//   NTick - # ticks to use in transform. To end if 0.
//   NormOpt - Normalization:
//     0 - Unnormalized ==> integrated raw power of transform is N times larger than input.
//         Normalization of 1/N is applied on inverse to recover original waveform.
//     1 - Normalized ==> Integrated power of transform or inverse is the same as its input.
//         I.e. normalization of 1/sqrt(N) is applied both forward and backward.
//     2 - Bin normalized ==> Integrated power of transform is 1/N times input.
//         Normalization of 1/N on forward transform and one on inverse.
//   Action -  0 - Use existing DFT
//             1 - Evaluate DFT from WF. Do not save result.
//             2 - Evaluate DFT if it is not already present. Do not save result.
//             3 - Evaluate DFT and save result
//             4 - Evaluate DFT if it is not already present and save result.
//            10 - Use existing WF.
//            11 - Evaluate WF from DFT. Do not save result.
//            12 - Evaluate WF if it is not already present. Do not save result.
//            13 - Evaluate WF and save result
//            14 - Evaluate WF if it is not already present and save result.
//   ReturnOpt - Controls how much data is written to the returned data map (see below)
//
// Only consistent results can be saved. I.e. FirstTick must be zero and Ntick zero
// or the same size as samples.
//
// In addition to the status, the returned data map includes:
//    For forward transforms (Action = 0-4):
//      if ReturnOpt%10 >= 1:
//        fftTick0 - First tick used to evaluate the DFT.
//        fftNTick - # ticks used to evaluate the DFT.
//        fftNMag - # DFT magnitudes
//        fftNPhase - # DFT phases (first is zero)
//        fftNSample- # samples = fftNMag + fftNPhase - 1
//      if ReturnOpt%10 >= 2:
//        fftMags - FFT magnitudes
//        fftPhases - FFT phases
//      if ReturnOpt%10 >= 3:
//        fftReals - Real part of each term (NTick entries)
//        fftImags - Imaginary part of each term (NTick entries)
//      if ReturnOpt = 10-13:
//        fftSamples - Waveform

#ifndef AdcChannelFFT_H
#define AdcChannelFFT_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AdcChannelFFT : AdcChannelTool {

public:

  using Index = unsigned int;
  using FloatVector = AdcSignalVector;

  AdcChannelFFT(fhicl::ParameterSet const& ps);

  ~AdcChannelFFT() override =default;

  // AdcChannelTool methods.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap update(AdcChannelData& acd) const override;

  // These methods do the work of calling Root interface to FFTW and converting to conventions here.
  // They could be moved to a utility class.
  int fftForward(Index normOpt, Index ntick, const float* psam,
                 FloatVector& xres, FloatVector& xims, FloatVector& mags, FloatVector& phases) const;
  int fftInverse(Index normOpt, const FloatVector& mags, const FloatVector& phases,
                 FloatVector& xres, FloatVector& xims, FloatVector& sams) const;

  // This does all the work of view:
  //   deciding what action to take
  //   calling one of the above  to carry out the action
  //   constructing the result
  void internalView(const AdcChannelData& acd, FloatVector& sams, FloatVector& mags, FloatVector& phas, DataMap& ret) const;

private:

  // Configuration data.
  int    m_LogLevel;
  Index  m_FirstTick;
  Index  m_NTick;
  Index  m_NormOpt;
  Index  m_Action;
  Index  m_ReturnOpt;

};

DEFINE_ART_CLASS_TOOL(AdcChannelFFT)

#endif
