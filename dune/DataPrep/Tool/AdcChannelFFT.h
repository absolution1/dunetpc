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
// of nsam terms all but the first and Nyquist (last only for even # samples) must be scaled
// by 1/sqrt(2).
//
// This is the Consistent-Power normalization described in DuneCommon/Utility/RealDftData.h.
//
// The DFT phases
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   FirstTick - First tick used in transform.
//   NTick - # ticks to use in transform. To end if 0.
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
//   DataView - Action is taken for each entry on this view. Blank means the top.
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
// If DatView is not blank, then the data map instead holds
//   fftNproc - # AdcChannelData objects processed
//   fftNfail - # objects for which the processing failed

#ifndef AdcChannelFFT_H
#define AdcChannelFFT_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneCommon/Utility/DuneFFT.h"

class AdcChannelFFT : TpcDataTool {

public:

  using Index = unsigned int;
  using FloatVector = AdcSignalVector;
  using DFT = DuneFFT::DFT;
  using Name = std::string;

  AdcChannelFFT(fhicl::ParameterSet const& ps);

  ~AdcChannelFFT() override =default;

  // AdcChannelTool methods.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap update(AdcChannelData& acd) const override;

private:

  // AdcChannelTool that act on the top of the data view.
  DataMap viewTop(const AdcChannelData& acd) const;
  DataMap updateTop(AdcChannelData& acd) const;

  // This does all the work of view:
  //   deciding what action to take
  //   calling one of the above  to carry out the action
  //   constructing the result
  void internalView(const AdcChannelData& acd, FloatVector& sams, FloatVector& amps, FloatVector& phas, DataMap& ret) const;

private:

  // Configuration data.
  int    m_LogLevel;
  Index  m_FirstTick;
  Index  m_NTick;
  Index  m_Action;
  Index  m_ReturnOpt;
  Name   m_DataView;

};


#endif
