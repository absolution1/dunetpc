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
//   6-9 - Options for direct deconvolution. Most successful is option 8 which
//         does chi-square minimization inluding a smoothing term.
// Note that an attempt to deconvolute with insufficient smearing, e.g. a delta
// function filter, is likely to result in large oscillations due to noise
// or numerical rounding.
//
// The response function is configured with an explicit vector that is padded
// zeroes to the length of the passed data. It is is implictly periodic.
// If empty, no action is taken.
//
// High frequncy filtering is provided by a Gaussian with time-domain sigma
// specified by a configuration parameter.
//
// Low-frequency filtering is provided by the function
//   1/[ 1 + 1/(W N k)^P ]
// N is the number of samples and and k is the wave number.
// The width W and power P are configuration parameters.
//   W < 0 disables the filter
//   W = 0 filters out k = 0 only
//   W > 0 removes k = 0 and applies the above for other values
//
// Other options may later be added for the response or filter.
//
// If an index map tool is provided, then it is used with the channel number  to choose
// the reponse vector and filter. Otherwise, the first entry is used.
//
// Configuration:
//   LogLevel: 0=silent, 1=init, 2=each event, >2=more
//   Action: Option for action described above.
//   ResponseVectors: Vector of discrete sequences representing the response function.
//   ResponseCenters: Vector specifying the zero tick for each of the response fundtions.
//   SmoothVectors: Vector o discrete sequnces holding smoothing functions.
//   SmoothScales: Scale factors for the smoothing functions: 1 if missing.
//   GausFilterSigmas: - Vector of time-domain sigmas [Tick] for the Gaussian filter.
//   LowFilterWidth: - Vector of widths [tick] [Tick] for the low-frequency filter.
//   LowFilterPower: - Vector of powers time-domain sigma [Tick] for the low-frequncy filter.
//   IndexMapTool - Name of the tool mapping channel to response index.

#ifndef AdcDeconvoluteFFT_H
#define AdcDeconvoluteFFT_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"

class AdcDeconvoluteFFT : TpcDataTool {

public:

  using Index = unsigned int;

  AdcDeconvoluteFFT(fhicl::ParameterSet const& ps);

  ~AdcDeconvoluteFFT() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Name = std::string;
  using IndexVector = std::vector<Index>;
  using ResponseVector = AdcSignalVector;   // Same type as AdcChannel data
  using ResponseVectorVector = std::vector<ResponseVector>;
  using FloatVector = std::vector<float>;
  using IndexMapToolPtr = std::unique_ptr<IndexMapTool>;

  // Configuration data.
  int                  m_LogLevel;
  Index                m_Action;
  ResponseVectorVector m_ResponseVectors;
  IndexVector          m_ResponseCenters;
  ResponseVectorVector m_SmoothVectors;
  FloatVector          m_SmoothScales;
  FloatVector          m_GausFilterSigmas;
  FloatVector          m_LowFilterWidths;
  FloatVector          m_LowFilterPowers;
  Name                 m_IndexMapTool;

  // Derived from configuration.
  bool m_useResponse;
  bool m_useFilter;
  bool m_doFFTConvolute;
  bool m_doDirectConvolute;
  bool m_doDeconvolute;
  int m_directDeconvolute;
  IndexMapToolPtr m_channelToIndex;

};


#endif
