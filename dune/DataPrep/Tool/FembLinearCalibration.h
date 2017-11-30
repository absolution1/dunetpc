// FembLinearCalibration.h

// David Adams
// November, 2017
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
//   nSample    - # of samples calibrated (size of raw data)
//   nUnderflow - # samples flagged as underflow (ADC <= AdcMin)
//   

#ifndef FembLinearCalibration_H
#define FembLinearCalibration_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelDataModifier.h"

class FembLinearCalibration
: public AdcChannelDataModifier {

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
