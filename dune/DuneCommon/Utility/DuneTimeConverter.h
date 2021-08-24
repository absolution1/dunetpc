// DuneTimeConverter

#ifndef DuneTimeConverter_H
#define DuneTimeConverter_H

// David Adams
// March 2016
//
// Utility that converts between NOvA time and DUNE time.
//
// DUNE uses art::Timestamp with
//   Unix time in sec in the high word
//   Remainder in ns in the low word
// NOvA time is the # ticks since 2010 with 1 tick = 1/(64 MHz),
// i.e. there are 64 ticks/us.

#include <stdint.h>
#include <string>
#include "canvas/Persistency/Provenance/Timestamp.h"

class DuneTimeConverter {
public:

  // Unix time when NOvA time is zero.
  static uint64_t novaT0Sec() {
    return 1262304000;  // 01Jan2010 UTC in unix time
  }
    
  // Number of NOvA ticks per microsecond.
  static uint64_t novaTicksPerUsec() {
    return 64;
  }
    
  // Number of NOvA ticks per second.
  static uint64_t novaTicksPerSec() {
    return 1000000*novaTicksPerUsec();
  }
    
  // Convert NOvA time to DUNE time.
  static art::Timestamp fromNova(uint64_t tnova);

  // Convert NOvA time to dune time.
  static uint64_t toNova(art::Timestamp tart);

  // Construct timestamp from the time in sec and remainder in ns
  static art::Timestamp makeTimestamp(uint32_t tsec, uint32_t trem);

  // Convert to a string with format "SSS.NNNNNNNNN" where SSS is the truncated
  // time in seconds and NNNNNNNNN is the remainder in ns.
  static std::string toString(art::Timestamp tart);

  // Convert from a string in the above format.
  static art::Timestamp fromString(std::string stime);

};

#endif
