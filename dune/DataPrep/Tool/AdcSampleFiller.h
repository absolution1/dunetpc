// AdcSampleFiller.h

// David Adams
// October 2017
//
// Tool to convert raw ADC values into samples.
//
// Subtracts the pedestal and records the sample in ADC units.
//
// If AdcOverflow > AdcUnderflow, then
//   Underflow flag is set for ADC value at or below AdcUnderflow
//   Overflow flag is set for ADC value at or above AdcOverflow
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   AdcUnderflow - value for underflow
//   AdcOverflow - value for underflow
//
// The output results holds:
//   int nUnderflow - # ticks at or below AdcUnderflow
//   int nOverflow - # ticks at or below AdcUnderflow
//   int nOutOfRange - # ticks below AdcUnderflow or above AdcOverflow
//
// Reads: raw
// Writes: samples, flags

#ifndef AdcSampleFiller_H
#define AdcSampleFiller_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <string>
#include <vector>

class HistogramManager;
class TH1;

class AdcSampleFiller
: public AdcChannelTool {

public:

  AdcSampleFiller(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int m_LogLevel;
  AdcIndex m_AdcUnderflow;
  AdcIndex m_AdcOverflow;

};

DEFINE_ART_CLASS_TOOL(AdcSampleFiller)

#endif
