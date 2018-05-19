// AdcChannelStringTool.h

// David Adams
// April 2018
//
// Interface for a tool that builds a string from AdcChannelData, 
// a DataMap and a string pattern.

#ifndef AdcChannelStringTool_H
#define AdcChannelStringTool_H

#include "dune/DuneInterface/AdcChannelData.h"
#include "dune/DuneInterface/Data/DataMap.h"

class AdcChannelStringTool {

public:

  using Index = unsigned int;

  // Helper function that calls tool from a pointer and returns the
  // input string if that pointer is null.
  static std::string
  build(const AdcChannelStringTool* ptool,
        const AdcChannelData& acd,  const DataMap& dm, std::string spat) {
    return ptool == nullptr ? spat : ptool->build(acd, dm, spat);
  }

  // Same without data map.
  static std::string
  build(const AdcChannelStringTool* ptool,
        const AdcChannelData& acd, std::string spat) {
    return ptool == nullptr ? spat : ptool->build(acd, DataMap(), spat);
  }
                        
  virtual ~AdcChannelStringTool() =default;

  virtual std::string
  build(const AdcChannelData& acd, const DataMap& dm, std::string spat) const =0;

};

#endif
