// FclStickyCodeFlagger.h

// David Adams
// October 2018
//
// Tool to flag sticky codes (i.e. set acd.flags) from fcl.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   StickyCodes - Array fcl-mapped vectors of sticky codes.
//   StickyRanges - Array of inclusive code ranges.
//
// Example syntax for the code map:
//   StickyCode: {
//     chan0001: [2002, 2004, 226, 2008]
//     chan0016: [2016]
//   }
//   StickyRange: {
//     chan0003: [0, 300]
//     chan0024: [1, 50]
//     chan0024: [3960, 4095]
//   }
//
// Channel numbers are encoded in the fcl element names and are obtained
// by stripping the leading "chan" followed by an arbitrary number of '0'
// and stripping an arbitrary number of trailing 'x'. The latter allow
// multiple entries for a single channel.

#ifndef FclStickyCodeFlagger_H
#define FclStickyCodeFlagger_H

#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include <vector>
#include <map>

class FclStickyCodeFlagger : TpcDataTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexVectorMap = std::map<Index, IndexVector>;
  using IndexPair = std::pair<Index, Index>;
  using IndexPairMap = std::multimap<Index, IndexPair>;
  using Name = std::string;

  FclStickyCodeFlagger(fhicl::ParameterSet const& ps);

  ~FclStickyCodeFlagger() override =default;

  DataMap update(AdcChannelData& acds) const override;

private:

  // Configuration data.
  int            m_LogLevel;
  AdcFlag        m_StickyCode;

  // Configuration derived data.
  IndexVectorMap m_stickyCodes;
  IndexPairMap m_stickyRanges;

};


#endif
