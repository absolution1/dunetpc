// AdcChannelData.h
//
// David Adams
// June 2016
//
// Struct to hold the prepared time samples for single TPC channel.

#ifndef AdcChannelData_H
#define AdcChannelData_H

#include <map>
#include "dune/DuneInterface/AdcTypes.h"

namespace raw {
  class RawDigit;
}

class AdcChannelData {

public:

  AdcChannel channel =-1;
  AdcSignalVector samples;
  AdcFlagVector flags;
  const raw::RawDigit* digit =nullptr;

};

typedef std::map<AdcChannel, AdcChannelData> AdcChannelDataMap;

#endif
