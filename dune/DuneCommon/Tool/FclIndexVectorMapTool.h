// FclIndexVectorMapTool.h
//
// David Adams
// July 2018
//
// Fcl-based tool to return Index vectors indexed by an index, e.g.
// event lists indexed by run.
//
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   Entries - vectors of indices in format [RRR, EE1, EE2, ...]
//             RRR is the map index
//             [EE1, EE2, ...] is the mapped vector
//             If a map index appeared ealier, the new entries are appended.

#ifndef FclIndexVectorMapTool_H
#define FclIndexVectorMapTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexVectorMapTool.h"
#include <map>

class FclIndexVectorMapTool : public IndexVectorMapTool {

public:

  using Index = IndexVectorMapTool::Index;
  using IndexVector = std::vector<Index>;
  using IndexVectorVector = std::vector<IndexVector>;
  using IndexVectorMap = std::map<Index, IndexVector>;

  // Ctor.
  FclIndexVectorMapTool(fhicl::ParameterSet const& ps);

  // Dtor.
  ~FclIndexVectorMapTool() override =default;

  // Return the vector for an index.
  IndexVector get(Index ient) const override;

private:

  // Parameters.
  Index m_LogLevel;
  IndexVectorVector m_Entries;
  IndexVectorMap m_entryMap;

};


#endif
