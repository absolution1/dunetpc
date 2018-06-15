// AdcResultDumper.h

// David Adams
// May 2018
//
// ADC channel tool that calls another ADC channel tool and dumps its result to
// the log and then returns that same result.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   Tool - Name of the called tool.

#ifndef AdcResultDumper_H
#define AdcResultDumper_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <memory>

class AdcResultDumper : public AdcChannelTool {

public:

  AdcResultDumper(fhicl::ParameterSet const& ps);

  ~AdcResultDumper() override =default;

  // AdcChannelTool interface.
  DataMap view(const AdcChannelData& acd) const override;
  DataMap update(AdcChannelData& acd) const override;
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  DataMap updateMap(AdcChannelDataMap& acds) const override;

private:

  // Configuration data.
  int            m_LogLevel;
  std::string    m_Tool;

  std::unique_ptr<AdcChannelTool> m_ptool;

};

DEFINE_ART_CLASS_TOOL(AdcResultDumper)

#endif
