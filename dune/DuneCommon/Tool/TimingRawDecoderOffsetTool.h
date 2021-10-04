// TimingRawDecoderOffsetTool.h
//
// David Adams
// June 2018
//
// Tool that calculates time offset using the timing clock received in an
// argumment.
//
// An earlier version of this tool read the timing clock from a text file
// written by the protoDUNE TimingRawDecoder module.
// 
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   TpcTickPhase - See below.
//   Unit - Unit in which the timing offset should be returned (daq, tick, ns).
//   FembScaleIds - IDs of the FEMBs fo which the offset should be scaled.
//   FembScaleValues - Scale factors for the scaled FEMBs.
//   RunDataTool - Name of the tool used to fetch the run data.
//
// For protoDUNE, daq is us/50, i.e 20 ns.
//
// The TPC tick is
//
//  t_TPC = scale*(t_daq + TpcTickPhase + runPhase)/TpcTickPeriod
//
// where t_daq is the count from the DAQ clock, TpcTickPeriod is (for now)
// fixed to 25 (50 MHz/2 MHz) and TpcTickPhase is a configuration parameter.
// Sensible values for the latter are 0, 1, 2, ..., TpcTickPeriod-1.
//
// The factor scale is otained from the corresponding position in FembScaleValues
// if the FEMB ID
//
// If the run data tool is found and has phaseGroup, then
//  runPhase = rdat.phases[igrp]

#ifndef TimingRawDecoderOffsetTool_H
#define TimingRawDecoderOffsetTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TimeOffsetTool.h"
#include <vector>
#include <map>

class RunDataTool;

class TimingRawDecoderOffsetTool : public TimeOffsetTool {

public:

  using Name = std::string;
  using IndexVector = std::vector<Index>;
  using DoubleVector = std::vector<double>;
  using ScaleMap = std::map<Index, double>;

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
  IndexVector m_FembScaleIds;
  DoubleVector m_FembScaleValues;
  Name m_RunDataTool;

  const RunDataTool*          m_pRunDataTool     =nullptr;

  // Derived from configuration.
  ScaleMap m_fembScales;    // Indexed by FEMB ID

};


#endif
