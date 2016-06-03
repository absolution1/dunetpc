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

class AdcChannelData {

public:

  AdcSignalVector samples;
  AdcFlagVector flags;

};

typedef std::map<AdcChannel, AdcChannelData> AdcChannelDataMap;

#endif
