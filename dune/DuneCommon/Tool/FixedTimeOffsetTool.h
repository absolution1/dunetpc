// FixedTimeOffsetTool.h
//
// David Adams
// June 2018
//
// Tool that returns a fixed time offset.
// 
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   Value - Integral part of the time offset
//   Rem - Floating part of the time offset
//   Unit - Unit for the time offset (e.g. daq, tick, ns).
//
// For protoDUNE, daq is us/50, i.e 20 ns.

#ifndef FixedTimeOffsetTool_H
#define FixedTimeOffsetTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TimeOffsetTool.h"
#include <vector>

class FixedTimeOffsetTool : public TimeOffsetTool {

public:

  using Name = std::string;

  // Ctor.
  FixedTimeOffsetTool(fhicl::ParameterSet const& ps);

  // Dtor.
  ~FixedTimeOffsetTool() override =default;

  // Return run data.
  Offset offset(const Data& dat) const override;

private:

  // Parameters.
  Index m_LogLevel;
  long m_Value;
  double m_Rem;
  Name m_Unit;

};

DEFINE_ART_CLASS_TOOL(FixedTimeOffsetTool)

#endif
