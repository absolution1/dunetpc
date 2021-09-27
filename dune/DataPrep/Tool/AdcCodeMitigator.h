// AdcCodeMitigator.h

// David Adams
// October 2018
//
// Tool to mitigate sticky codes (i.e. set acd.flags) from fcl.
//
// If the flag for a sample is in FixFlags, then the sample is set to zero.
// If the flag is in InterpolateFlags, then the sample is obtained by
// interpolating between the first non-sticky samples before and after.
//
// The algorithm for interpolation depends on the data and the the
// curvature threshold FCT:
//
// 1. If FCT > 0 and two non-sticky codes are found on at least one side
// of the mitigated samples, the signal jump is evaluted on each side
// difference between two nearest non-sticky samples. A constant-curvature
// interpolation is done if the magnitude of either of these jumps exceeds FCT.
//
// The interpolation passes through the closest point on each side and equally
// far from the second point on each side.
//
// 2. Otherwise and if one non-sticky code is found on either side of the
// mitigated samples, linear interpolate is performed between the nearest
// samples on each side.
//
// 3. Otherwise if a non-sticky sample is found on one side, that value is used.
//
// 4. Otherwise the mitigated samples are set to zero.
//
// Here non-sticky means not in InterpolateFlags or SkipFlags.
//
// The fix mitigation is done before the interpolation.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   FixFlags - Samples with these flags are set to zero
//   InterpolateFlags - Samples with these flags are interpolated from neighbors
//   SkipFlags - Samples with these flags are not used for interpolation.
//   FixedCurvThresh - The parameter FCT discussed above.

#ifndef AdcCodeMitigator_H
#define AdcCodeMitigator_H

#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include <vector>
#include <set>
#include <map>

class AdcCodeMitigator : TpcDataTool {

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
  double         m_FixedCurvThresh;

  // Configuration derived data.
  IndexSet m_fixSet;
  IndexSet m_interpolateSet;
  IndexSet m_skipSet;

};


#endif
