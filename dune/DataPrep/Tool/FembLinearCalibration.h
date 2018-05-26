// FembLinearCalibration.h

// David Adams
// November 2017
//
// Applies a linear FEMB-based calibration to ADC channel data.
// The sample charge is
//   q = Gains[icha]*(adc - ped)
// where adc is the raw value, ped is the pedestal and Gains is a calibration
// parameter.
// The sample flag is set to underflow if the raw ADC value is at or below
// the calibration parameter AdcMins[icha].
// The sample flag is set to overflow if the raw ADC value is at or above
// the calibration parameter AdcMaxs[icha].
// If AdcMins has length zero, then AdcMin is used for all channels.
// If AdcMaxs has length zero, then AdcMax is used for all channels.
//
// If the gain is zero or negative, then no calibration is performed and the
// prepared data vector is cleared.
//
// Calibration parameters are read from fcl files with the name FclNameBaseFFF
// where FFF is the FEMB ID.
//
// Parameters:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   Units - units for the calibrated samples ("fC", "ke", "mV", ...)
//   Gains - gain for each channel
//   AdcMin - min ADC if AdcMins is empty
//   AdcMins - min ADC for each channel
//   AdcMax - max ADC if AdcMaxs is empty
//   AdcMaxs - max ADC for each channel
//
// Output:
//   calibSampleCount    - # of samples calibrated (size of raw data)
//   calibUnderflowCount - # samples flagged as underflow (ADC <= AdcMin)
//   calibOverflowCount  - # samples flagged as underflow (ADC >= AdcMax)
//   calibAdcMin         - min ADC for the channel
//   calibAdcMax         - max ADC for the channel
//   calibGain           - applied gain
//
// Int he case of zero or negative gain, only the last is returned.
//   

#ifndef FembLinearCalibration_H
#define FembLinearCalibration_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class FembLinearCalibration
: public AdcChannelTool {

public:

  FembLinearCalibration(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Parameters.
  int m_LogLevel;
  std::string m_Units;
  AdcSignalVector m_Gains;
  AdcCount m_AdcMin;
  AdcCountVector m_AdcMins;
  AdcCount m_AdcMax;
  AdcCountVector m_AdcMaxs;

};

DEFINE_ART_CLASS_TOOL(FembLinearCalibration)

#endif
