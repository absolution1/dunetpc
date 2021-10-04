// AdcChannelSplitter.h

// David Adams
// September 2019
//
// Tool to split ADC samples to a specified length.
//
// Splits only raw if there are no samples, only samples if there is no
// raw or both if they have the same length.
//
// Configuration:
//   LogLevel: 0=errors only, 1=show config, 2=message for each trim/pad
//   Length - Output length.
//   DataPath: Data view path where the split is applied.
//   DataView: View where the split is written.
//
// Output datamap:
//   int splitInputCount: Number input objects.
//   int splitOutputCount: Number of created objects.
//   int splitRawCopiedCount: Number of copied raw sample ticks.
//   int splitSampleCopiedCount: Number of copied sample ticks.

#ifndef AdcChannelSplitter_H
#define AdcChannelSplitter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "TH1.h"
#include <memory>

class AdcChannelStringTool;
class TPadManipulator;

class AdcChannelSplitter : public TpcDataTool {

public:

  using Index = unsigned int;
  using Name = std::string;

  AdcChannelSplitter(fhicl::ParameterSet const& ps);

  ~AdcChannelSplitter() override =default;

  // Inherited methods.
  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int    m_LogLevel;
  Index  m_Length;
  Name   m_DataPath;
  Name   m_DataView;

};


#endif
