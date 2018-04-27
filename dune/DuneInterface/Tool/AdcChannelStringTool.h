// AdcChannelStringTool.h

// David Adams
// April 2018
//
// Interface for a tool that builds a string from AdcChannelData and
// a string pattern and a count.

#ifndef AdcChannelStringTool_H
#define AdcChannelStringTool_H

#include "dune/DuneInterface/AdcChannelData.h"

class AdcChannelStringTool {

public:

  using Index = unsigned int;

  // Helper function that calls tool from a pointer and returns the
  // input string if that pointer is null.
  static std::string
  build(const AdcChannelStringTool* ptool,
        const AdcChannelData& acd, std::string spat,
        Index count =AdcChannelData::badIndex) {
    return ptool == nullptr ? spat : ptool->build(acd, spat, count);
  }
                        
  virtual ~AdcChannelStringTool() =default;

  virtual std::string
  build(const AdcChannelData& acd, std::string spat,
        Index count =AdcChannelData::badIndex) const =0;

};

#endif
