// TimingRawDecoderOffsetTool.h
//
// David Adams
// June 2018
//
// Tool that reads the time offset from the text file written by
// the protoDUNE TimingRawDecoder module.
// 
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   TpcTickPhase - See below.
//   Unit - Unit in which the timing offset should be returns (daq, tick, ns).
//
// For protoDUNE, daq is us/50, i.e 20 ns.
//
// The TPC tick is
//
//  t_TPC = (t_daq + TpcTickPhase)/TpcTickPeriod
//
// where t_daq is the count from the DAQ clock, TpcTickPeriod is (for now)
// fixed to 25 (50 MHz/2 MHz) and TpcTickPhase is a configuration parameter.
// Sensible values for the latter are 0, 1, 2, ..., TpcTickPeriod-1.

#ifndef TimingRawDecoderOffsetTool_H
#define TimingRawDecoderOffsetTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TimeOffsetTool.h"
#include <vector>

class TimingRawDecoderOffsetTool : public TimeOffsetTool {

public:

  using Name = std::string;

  // Ctor.
  TimingRawDecoderOffsetTool(fhicl::ParameterSet const& ps);

  // Dtor.
  ~TimingRawDecoderOffsetTool() override =default;

  // Return run data.
  Offset offset(const Data& dat) const override;

private:

  // Parameters.
  Index m_LogLevel;
  Index m_TpcTickPhase;
  Name m_Unit;

};

DEFINE_ART_CLASS_TOOL(TimingRawDecoderOffsetTool)

#endif
