// DuneEventInfo.h

// David Adams
// June 2019
//
// Event info for dataprep in DUNE.

#ifndef DuneEventInfo_H
#define DuneEventInfo_H

#include <sstream>

class DuneEventInfo {

public:

  using Index = unsigned int;
  using LongIndex = unsigned long;

  Index run =0;
  Index event = 0;
  Index subRun = 0;
  time_t time =0;
  int timerem =0;
  Index trigger =0;
  LongIndex triggerClock =0;
  LongIndex channelClock =0;

  // Reset to no event.
  void clear() {
    run = 0;
    event = 0;
    subRun = 0;
    time = 0;
    timerem = 0;
    trigger = 0;
    triggerClock = 0;
    triggerClock = 0;
  }

  // Is this valid info?
  bool isValid() const { return run > 0; }

  // Order two event info objects.
  bool operator<(const DuneEventInfo& rhs) {
    if ( run < rhs.run ) return true;
    if ( rhs.run < run ) return false;
    if ( subRun < rhs.subRun ) return true;
    if ( rhs.subRun < subRun ) return false;
    if ( event < rhs.event ) return true;
    if ( rhs.event < event ) return false;
    if ( triggerClock < rhs.triggerClock ) return true;
    if ( rhs.triggerClock < triggerClock ) return false;
    return false;
  }

  // Compare two event info objects.
  bool operator==(const DuneEventInfo& rhs) const {
    if ( run != rhs.run ) return false;
    if ( subRun != rhs.subRun ) return false;
    if ( event != rhs.event ) return false;
    if ( time != rhs.time ) return false;
    if ( timerem != rhs.timerem ) return false;
    if ( trigger != rhs.trigger ) return false;
    if ( triggerClock != rhs.triggerClock ) return false;
    if ( channelClock != rhs.channelClock ) return false;
    return true;
  }

  bool operator!=(const DuneEventInfo& rhs) const {
    return ! (*this == rhs);
  }

  // Return the event number as a string.
  std::string eventString() const {
    std::ostringstream ssout;
    ssout << event;
    return ssout.str();
  }

  // For run RRR and subrun SSS:
  //   opt = 0:  Return RRR.
  //         1:  Return RRR-SSS
  //         2:  if SSS==0, return  RRR, otherwise RRR-SSS.
  std::string runString(Index opt =2) const {
    std::ostringstream ssout;
    ssout << run;
    if ( opt == 1 || (opt ==2 && subRun > 0) ) {
      ssout << "-" << subRun;
    }
    return ssout.str();
  }

};

#endif
