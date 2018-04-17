// ToolBasedRawDigitPrepService.h
//
// David Adams
// April 2018
//
// Implementation of service that prepares TPC raw digits for reconstruction.
// It receives an ADC channel data map, applies a sequence of ADC channel
// tools (calling updateMap) and then constructs wires.
//
// Confguration parameters.
//   LogLevel - logging level
//              0 - errors only
//              1 - initilization information including configuration
//              2 - one line per processed map
//              3 - one line for each step (tool call)
//              4 - display result from each step
//   DoWires - If true, the wire building service is called after processing.
//   AdcChannelToolNames - Names of the ADC channel tools.

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
