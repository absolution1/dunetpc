// DuneEventInfo.h

// David Adams
// June 2019
//
// Event info for dataprep in DUNE.

#ifndef DuneEventInfo_H
#define DuneEventInfo_H

#include <sstream>
#include <memory>

class DuneEventInfo {

public:

  using Index = unsigned int;
  using LongIndex = unsigned long;

  Index run =badIndex();
  Index event =badIndex();
  Index subRun =badIndex();
  time_t time =0;
  int timerem =0;
  Index trigger =badIndex();
  LongIndex triggerClock =badLongIndex();
  LongIndex triggerTick0 =badLongIndex();

  static Index badIndex() { return -1u; }
  static LongIndex badLongIndex() { return -1ul; }
  static const DuneEventInfo& badEventInfo() {
    static DuneEventInfo bei;
    return bei;
  }

  // Reset to no event.
  void clear() {
    run = badIndex();
    event = badIndex();
    subRun = badIndex();
    time = 0;
    timerem = 0;
    trigger = badIndex();
    triggerClock = badLongIndex();
    triggerTick0 = badLongIndex();
  }

  // Default ctor.
  DuneEventInfo() = default;

  // Ctor from components.
  DuneEventInfo(Index a_run, Index a_event, Index a_subRun =badIndex(),
                time_t a_time =0, int a_timerem =0,
                Index a_trigger =badIndex(),
                LongIndex a_triggerClock =badLongIndex(),
                LongIndex a_triggerTick0 =badLongIndex())
  : run(a_run), event(a_event), subRun(a_subRun),
    time(a_time), timerem(a_timerem),
    trigger(a_trigger),
    triggerClock(a_triggerClock),
    triggerTick0(a_triggerTick0) { }

  // Is this valid info?
  bool isValid() const { return run != badIndex(); }

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
    if ( triggerTick0 < rhs.triggerTick0 ) return true;
    if ( rhs.triggerTick0 < triggerTick0 ) return false;
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
    if ( triggerTick0 != rhs.triggerTick0 ) return false;
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
    if ( opt == 1 || (opt ==2 && subRun > 0 && subRun!=badIndex()) ) {
      ssout << "-" << subRun;
    }
    return ssout.str();
  }

};

#endif
