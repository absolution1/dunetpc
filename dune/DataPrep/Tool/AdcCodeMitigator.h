// AdcCodeMitigator.h

// David Adams
// October 2018
//
// Tool to mitigate sticky codes (i.e. set acd.flags) from fcl.
//
// If the flag for a sample is in FixFlags, then the sample is set to zero.
// If the flag is in InterpolateFlags, then the sample is obtained by
// linear interpolation between the first non-sticky samples before and after.
// Here non-sticky mean not in InterpolateFlags or SkipFlags.
//
// The fix mitigation is done before the interpolation.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   FixFlags - Samples with these flags are set to zero
//   InterpolateFlags - Samples with these flags are interpolated from neighbors
//   SkipFlags - Samples with these flags are not used for interpolation.

#ifndef AdcCodeMitigator_H
#define AdcCodeMitigator_H

#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include <vector>
#include <set>
#include <map>

class AdcCodeMitigator : AdcChannelTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;
  using Name = std::string;

  AdcCodeMitigator(fhicl::ParameterSet const& ps);

  ~AdcCodeMitigator() override =default;

  DataMap update(AdcChannelData& acds) const override;

private:

  // Configuration data.
  int            m_LogLevel;
  IndexVector    m_FixFlags;
  IndexVector    m_InterpolateFlags;
  IndexVector    m_SkipFlags;

  // Configuration derived data.
  IndexSet m_fixSet;
  IndexSet m_interpolateSet;
  IndexSet m_skipSet;

};

DEFINE_ART_CLASS_TOOL(AdcCodeMitigator)

#endif
