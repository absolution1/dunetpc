// ToolBasedRawDigitPrepService.h
//
// David Adams
// April 2018
//
// Implementation of service that prepares TPC raw digits for reconstruction.
// It receives an ADC channel data map, applies a sequence of ADC channel
// tools and then constructs wires.

#ifndef ToolBasedRawDigitPrepService_H
#define ToolBasedRawDigitPrepService_H

#include "dune/DuneInterface/RawDigitPrepService.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <map>

class AdcWireBuildingService;
class AdcChannelDataCopyService;

class ToolBasedRawDigitPrepService : public RawDigitPrepService {

public:

  using AdcChannelToolPtr = std::unique_ptr<AdcChannelTool>;
  struct NamedTool {
    std::string name;
    const AdcChannelTool* tool;
    NamedTool(std::string a_name ="", const AdcChannelTool* a_tool =nullptr) : name(a_name), tool(a_tool) { }
  };
  using AdcChannelToolVector = std::vector<AdcChannelToolPtr>;
  using AdcChannelNamedToolVector = std::vector<NamedTool>;

  ToolBasedRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int prepare(AdcChannelDataMap& prepdigs,
              std::vector<recob::Wire>* pwires,
              WiredAdcChannelDataMap* pintStates) const override;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const override;

private:

  // Configuration parameters.
  int m_LogLevel;
  bool m_DoWires;
  std::vector<std::string> m_AdcChannelToolNames;

  AdcChannelToolVector m_AdcChannelTools;
  AdcChannelNamedToolVector m_AdcChannelNamedTools;
  const AdcWireBuildingService* m_pWireBuildingService;


};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ToolBasedRawDigitPrepService, RawDigitPrepService, LEGACY)

#endif
