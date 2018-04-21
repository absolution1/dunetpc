// AdcUnderflowFlagger.h

// David Adams
// October 2017
//
// Tool to flag underflows in ADC data.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=minor errors each event, >2=more
//   AdcThresholds - vector with threshold for each channel
//   DefaultThreshold - If this is >= 0, then this value is used for channels
//                      that are not set in (negative) or beyond the range of
//                      the AdcThresholds vector.
//                      Otherwise no action is taken and an error is returned.
//
// Output in DataMap:
//   status: 0 - success
//           1 - channel not set in ADC data
//           2,3 - threshold not set for the requested channel
//   int channel - channel number
//   int nUnderflow - # ticks flagged as underflow
//   int nModify - # calls to modify methods

#ifndef AdcUnderflowFlagger_H
#define AdcUnderflowFlagger_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <string>
#include <vector>

class HistogramManager;
class TH1;

class AdcUnderflowFlagger
: public AdcChannelTool {

public:

  using IntVector = std::vector<short>;

  AdcUnderflowFlagger(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

  // Allow user to modify thresholds on the fly.
  AdcCountVector& adcThreshold() { ++m_modifyCount; return m_AdcThresholds; }
  const AdcCountVector& adcThreshold() const { return m_AdcThresholds; }
  int& defaultThreshold() { ++m_modifyCount; return m_DefaultThreshold; }
  AdcCount defaultThreshold() const { return m_DefaultThreshold; }

private:

  // Configuration data.
  int m_LogLevel;
  IntVector m_AdcThresholds;
  int m_DefaultThreshold;

  // Counter for access to modifiers.
  AdcIndex m_modifyCount;


};

DEFINE_ART_CLASS_TOOL(AdcUnderflowFlagger)

#endif
