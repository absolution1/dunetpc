// FloatArrayGainCalibration.h

// David Adams
// February 2019
//
// Applies a gain calibration to ADC channel data.
// The assigned sample charge is
//   q = Gains[icha]*(adc - ped)
// where adc is the raw value, ped is the pedestal and Gains is taken from
// a FloatArrayTool.
//
//   Underflow flag is set for ADC value at or below AdcUnderflow
//   Overflow flag is set for ADC value at or above AdcOverflow
//
// Parameters:
//   LogLevel - 0: silent
//              1: display configuration in ctor
//              2: display when channel is out of range
//   Unit - units for the calibrated samples ("fC", "ke", "mV", ...)
//          If "derived", then /(ADC count) is removed from the gain unit.
//   GainDefault - Gain used if input gain is <=0 or missing.
//                 If this value is negative, then the default value of
//                 GainTool is used instead (Feb 2021).
//   AdcUnderflowDefault - value for underflow
//   AdcOverflowDefault - value for overflow
//   GainTool - Name of the FloatArrayTool holding the gains.
//   ScaleFactor - Scale factor formula for the gain (tool or default)
//
// Output:
//   calibSampleCount    - # of samples calibrated (size of raw data)
//   calibUnderflowCount - # samples flagged as underflow (ADC <= AdcMin)
//   calibOverflowCount  - # samples flagged as underflow (ADC >= AdcMax)
//   calibAdcMin         - min ADC for the channel
//   calibAdcMax         - max ADC for the channel
//   calibGain           - applied gain
//
// In the case of zero or negative gain, only the last is returned.
//   

#ifndef FloatArrayGainCalibration_H
#define FloatArrayGainCalibration_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Utility/ParFormula.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"

class FloatArrayTool;
class RunDataTool;

class FloatArrayGainCalibration
: public TpcDataTool {

public:

  using Name = std::string;
  using FloatArrayPtr = const FloatArrayTool*;

  FloatArrayGainCalibration(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration parameters.
  int m_LogLevel;
  Name m_Unit;
  ParFormula* m_GainDefault;
  AdcIndex m_AdcUnderflowDefault;
  AdcIndex m_AdcOverflowDefault;
  Name m_GainTool;
  ParFormula* m_ScaleFactor;

  // Derived parameters.
  FloatArrayPtr m_pgains;
  const RunDataTool* m_prdtool; // Run data tool.

};


#endif
