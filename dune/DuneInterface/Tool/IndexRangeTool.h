// IndexRangeTool.h

// David Adams
// July 2018
//
// Interface for tools that that provide access IndexRange objects.

#ifndef IndexRangeTool_H
#define IndexRangeTool_H

#include "dune/DuneInterface/Data/IndexRange.h"

class IndexRangeTool {

public:

  using Name = IndexRange::Name;

  virtual ~IndexRangeTool() =default;

  virtual IndexRange get(Name nam) const =0;

};

#endif
