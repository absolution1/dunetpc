// ChannelStatusConfigTool.h
//
// David Adams
// April 2020
//
// Tool to return channel status from fcl configuration that
// is eaily copied from that used in the channel status service.
// Example configuration:
//
// cstool.IndexVectors: [
//   [],
//   @local::services.ChannelStatusService.BadChannels,
//   @local::services.ChannelStatusService.NoisyChannels,
//   myListOfOtherBadChannels
// ]
//
// Parameters:
//   LogLevel - Message logging level:
//              0: none
//              1: show config
//              2: check for duplicates
//   DefaultIndex - Index for values that appear in none of the index vector lists.
//   IndexVectors - Vector of index vectors. Values for the ith entry
//                  are mapped to index i.
//
// A channel appearing in multiple lists is assigned to the first in which
// it appears.

#ifndef ChannelStatusConfigTool_H
#define ChannelStatusConfigTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include <vector>

class ChannelStatusConfigTool : public IndexMapTool {

public:

  using Index = IndexMapTool::Index;
  using IndexVector = std::vector<Index>;
  using IndexVectorVector = std::vector<IndexVector>;

  // Ctor.
  ChannelStatusConfigTool(fhicl::ParameterSet const& ps);

  // Dtor.
  ~ChannelStatusConfigTool() override =default;

  // Return the channel status.
  Index get(Index icha) const override;

private:

  // Parameters.
  Index m_LogLevel;
  Index m_DefaultIndex;
  IndexVectorVector m_IndexVectors;

  // Derived from configuration.
  IndexVector m_vals;

};


#endif
