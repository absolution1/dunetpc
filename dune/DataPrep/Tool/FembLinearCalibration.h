// FembLinearCalibration.h

// David Adams
// November, 2017
//
// Applies a linear FEMB-based calibration to ADC channel data.
// The sample charge is
//   q = Gain*(adc - ped)
// where adc is the raw value, ped is the pedestal and Gain is a calibration
// parameter.
// The sample flag is set to underflow if the raw ADC value is at or below
// the calibration parameter AdcMin.
//
// Calibration parameters are read from fcl files with the name FclNameBaseFFF
// where FFF is the FEMB ID.
//
// Parameters:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   Gains - gain for each channel
//   AdcMin - min ADC for each channel
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
  AdcSignalVector m_Gains;
  AdcCountVector m_AdcMins;

};

DEFINE_ART_CLASS_TOOL(FembLinearCalibration)

#endif
