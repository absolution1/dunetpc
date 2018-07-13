// FclIndexRangeTool.h
//
// David Adams
// July 2018
//
// Fcl-based tool to return IndexRange objects.
//
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   Ranges - Vector fcl blocks each describing a range:
//              name - name for the range
//             label - label
//             begin - first index
//               end - last+1 index

#ifndef FclIndexRangeTool_H
#define FclIndexRangeTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include <map>

class FclIndexRangeTool : public IndexRangeTool {

public:

  using Index = IndexRange::Index;
  using IndexRangeMap = std::map<Name, IndexRange>;

  // Ctor.
  FclIndexRangeTool(fhicl::ParameterSet const& ps);

  // Dtor.
  ~FclIndexRangeTool() override =default;

  // Return a range.
  IndexRange get(Name nam) const override;

private:

  // Parameters.
  Index m_LogLevel;
  IndexRangeMap m_Ranges;

};

DEFINE_ART_CLASS_TOOL(FclIndexRangeTool)

#endif
