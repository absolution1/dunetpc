// AdcChannelData.h
//
// David Adams
// June 2016
//
// Struct to hold the prepared time samples for single TPC channel.
//
//    channel - Offline channel number
//   pedestal - Pedestal subtracted from the raw count
//    samples - Array holding the prepared signal value for each tick
//      flags - Array holding the status flag for each tick
//      digit - Corresponding raw digit

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
  AdcSignal pedestal = -99999;
  AdcSignalVector samples;
  AdcFlagVector flags;
  const raw::RawDigit* digit =nullptr;

};

typedef std::map<AdcChannel, AdcChannelData> AdcChannelDataMap;

#endif
