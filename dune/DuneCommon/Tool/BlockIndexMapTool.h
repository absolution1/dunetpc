// BlockIndexMapTool.h
//
// David Adams
// November 2020
//
// Tool to map a block of indices (e.g. a channel) to an index (e.g. wire plane index).
//
// Parameters:
//   LogLevel - Message logging level:
//              0: none
//              1: show config
//              2: check for duplicates
//   Starts: Index for which each block starts
//   Values: Value to which each blokc of indices is mapped
//   Period: Mapping is repeated with this period
//
// The starts and values do not need to be ordered.
// If any starts are duplicated, the value for last one takes precedence.
// If there are fewer values than starts, value zer is used for the extra starts.
// If Period is nonzero the  mapping is done with the remainder after dividing by that value.

#ifndef BlockIndexMapTool_H
#define BlockIndexMapTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include <vector>
#include <map>

class BlockIndexMapTool : public IndexMapTool {

public:

  using Index = IndexMapTool::Index;

  // Ctor.
  BlockIndexMapTool(fhicl::ParameterSet const& ps);

  // Dtor.
  ~BlockIndexMapTool() override =default;

  // Return the channel status.
  Index get(Index icha) const override;

private:

  using IndexVector = std::vector<Index>;
  using IndexMap = std::map<Index, Index>;

  // Parameters.
  Index m_LogLevel;
  IndexVector m_Starts;
  IndexVector m_Values;
  Index       m_Period;

  // Derived data.
  IndexMap m_map;

};


#endif
