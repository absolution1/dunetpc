// FclStickyCodeFlagger.h

// David Adams
// October 2018
//
// Tool to flag sticky codes (i.e. set acd.flags) from fcl.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more

#ifndef FclStickyCodeFlagger_H
#define FclStickyCodeFlagger_H

#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include <vector>
#include <map>

class FclStickyCodeFlagger : AdcChannelTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexMap = std::map<Index, IndexVector>;
  using Name = std::string;

  FclStickyCodeFlagger(fhicl::ParameterSet const& ps);

  ~FclStickyCodeFlagger() override =default;

  DataMap update(AdcChannelData& acds) const override;

private:

  // Configuration data.
  int            m_LogLevel;
  AdcFlag        m_StickyCode;

  // Configuration derived data.
  IndexMap m_stickyCodes;

};

DEFINE_ART_CLASS_TOOL(FclStickyCodeFlagger)

#endif
