// AdcChannelTrimmer.h

// David Adams
// May 2019
//
// Tool to trim or pad the ADC samples to a fixed length. If longer,
// later samples are dropped. If shorter, samples at the start are
// appended to the end.
//
// If the # samples to trim or pad exceeds a threshold, no action is taken
// and an error message is broadcast.
//
// Configuration:
//   LogLevel: 0=errors only, 1=show config, 2=message for each trim/pad
//   Length - Output length.
//   MaxTrim - Max # samples to trim or pad.
//
// Output datamap:
//   int trimAction: 0 - Input # samples matched the requested value
//                   1 - Samples vector trimmed to requested length
//                  -1 - From call to view. Vector would have been trimmed if
//                       update was called.
//                   2 - Required trim exceeded threshold.
//                   3 - Input data had no samples.
//   int trimLength: (# input samples) - (requested length)

#ifndef AdcChannelTrimmer_H
#define AdcChannelTrimmer_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "TH1.h"
#include <memory>

class AdcChannelStringTool;
class TPadManipulator;

class AdcChannelTrimmer : public TpcDataTool {

public:

  using Index = unsigned int;

  AdcChannelTrimmer(fhicl::ParameterSet const& ps);

  ~AdcChannelTrimmer() override =default;

  // Inherited methods.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int    m_LogLevel;
  Index  m_Length;
  Index  m_MaxTrim;

};


#endif
