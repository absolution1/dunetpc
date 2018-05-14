// RunDataTool.h

#ifndef RunDataTool_H
#define RunDataTool_H

// Interface for a tool providing access to run conditions data.

#include "dune/DuneInterface/Data/RunData.h"

class RunDataTool {

public:

  using Index = unsigned int;

  virtual ~RunDataTool() =default;

  virtual RunData runData(Index run, Index subRun =0) const =0;

};

#endif
