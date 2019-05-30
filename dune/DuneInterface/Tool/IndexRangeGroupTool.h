// IndexRangeGroupTool.h

// David Adams
// May 2019
//
// Interface for tools that that provide access IndexRangeGroup objects.

#ifndef IndexRangeGroupTool_H
#define IndexRangeGroupTool_H

#include "dune/DuneInterface/Data/IndexRangeGroup.h"

class IndexRangeGroupTool {

public:

  using Name = IndexRangeGroup::Name;

  virtual ~IndexRangeGroupTool() =default;

  virtual IndexRangeGroup get(Name nam) const =0;

};

#endif
