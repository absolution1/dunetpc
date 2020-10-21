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
#include <chrono>

class AdcWireBuildingService;
class AdcChannelDataCopyService;

class ToolBasedRawDigitPrepService : public RawDigitPrepService {

public:

  using Index = unsigned int;
  using AdcChannelToolPtr = std::unique_ptr<AdcChannelTool>;
  struct NamedTool {
    std::string name;
    const AdcChannelTool* tool;
    NamedTool(std::string a_name ="", const AdcChannelTool* a_tool =nullptr) : name(a_name), tool(a_tool) { }
  };
  using AdcChannelToolVector = std::vector<AdcChannelToolPtr>;
  using AdcChannelNamedToolVector = std::vector<NamedTool>;

  ToolBasedRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);
  ~ToolBasedRawDigitPrepService();

  // Called at begin and end of event processing.
  // Calls the same method for each tool.
  int beginEvent(const DuneEventInfo& devt) const override;
  int endEvent(const DuneEventInfo& devt) const override;

  int prepare(detinfo::DetectorClocksData const& clockData,
              AdcChannelDataMap& prepdigs,
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

  using Clock = std::chrono::steady_clock;
  using Duration = std::chrono::duration<double>;
  class State {
  public:
    // Timing.
    std::vector<Duration> toolTimes;
  };
  std::unique_ptr<State> m_pstate;
  State& state() const { return *m_pstate; }

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ToolBasedRawDigitPrepService, RawDigitPrepService, LEGACY)

#endif
