// AdcDeconvoluteFFT.h
//
// Tool to perform deconvolution or convolution on AdcChannelData.
// The tool is configured with a response function and a filter function.
// It updates the samples in a passed AdcChannelData as specified by action:
//   0 - No change in samples
//   1 - Deconvolute with the reponse and convolute with filter.
//   2 - Convolute with the response using FFT. Filter is ignored.
//   3 - Convolute with the filter using FFT. Response is ignored.
//   4 - Convolute with the response directly. Filter is ignored.
//   5 - Convolute with the filter directly. Response is ignored.
// Note that an attempt to deconvolute with insufficient smearing, e.g. a delta
// function filter, is likely to result in large oscillations due to noise
// or numerical rounding.
//
// The response function is configured with an explicit vector that is padded
// zeroes to the length of the passed data. It is is implictly periodic.
//
// The filter function is Gaussian with time-domain sigma provided in
// the configuration.
//
// Other options may later be added for the reponse or filter.
//
// Configuration:
//   LogLevel: 0=silent, 1=init, 2=each event, >2=more
//   Action: Option for action described above.
//   ResponseVector: Discrete sequence representing the response function.
//   GausFilterRms - Time-domain sigma [Tick] for the Gaussian filter.

#ifndef AdcDeconvoluteFFT_H
#define AdcDeconvoluteFFT_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AdcDeconvoluteFFT : AdcChannelTool {

public:

  using Index = unsigned int;

  AdcDeconvoluteFFT(fhicl::ParameterSet const& ps);

  ~AdcDeconvoluteFFT() override =default;

  DataMap update(AdcChannelData& acd) const override;

  // Internal methods.
  // Get the
private:

  // Configuration data.
  int             m_LogLevel;
  Index           m_Action;
  AdcSignalVector m_ResponseVector;
  float           m_GausFilterSigma;

  // Derived from configuration.
  bool m_useResponse;
  bool m_useFilter;
  bool m_doFFTConvolute;
  bool m_doDirectConvolute;
  bool m_doDeconvolute;

};

DEFINE_ART_CLASS_TOOL(AdcDeconvoluteFFT)

#endif
