// AdcChannelStringTool.h

// David Adams
// April 2018
//
// Interface for a tool that builds a string from AdcChannelData and
// a string pattern and a count.

#ifndef AdcChannelTool_H
#define AdcChannelTool_H

#include "dune/DuneInterface/AdcChannelData.h"

class AdcChannelStringTool {

public:

  using Index = unsigned int;

  virtual ~AdcChannelStringTool() =default;

  virtual std::string
  build(const AdcChannelData& acd, std::string spat,
        Index count =AdcChannelData::badIndex) const =0;

};

#endif
