//File: Timer.cpp
//Brief: A CRT::Timer approximates the time represented by a CRT timestamp and the elapsed time since an initial time it was given. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

//Include header
#include "Timer.h"

//Include definition of CRT::Trigger for interface
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

namespace CRT
{
  Timer::Timer(const CRT::Trigger& trigger): Timer(trigger.Timestamp())
  {
  }

  //TODO: Fill in constructor
  Timer::Timer(const uint64_t timestamp): fLastUNIXSecond(Timestamp(timestamp).upper), fNsOfLastUNIXSecond(Timestamp(timestamp).lower), 
                                          fNsOfLastSync(0), fPrevNs(Timestamp(timestamp).lower), fFirstTime(GetSeconds(timestamp))
  {
  }

  double Timer::seconds(const CRT::Trigger& trigger)
  {
    return seconds(trigger.Timestamp());
  }

  double Timer::seconds(const uint64_t timestamp)
  {
    return GetSeconds(timestamp);
  }

  double Timer::elapsed(const CRT::Trigger& trigger)
  {
    return elapsed(trigger.Timestamp());
  }

  double Timer::elapsed(const uint64_t timestamp)
  {
    return GetSeconds(timestamp) - fFirstTime;
  }

  double Timer::GetSeconds(const uint64_t timestamp)
  {
    //Extract upper and lower parts of timestamp
    Timestamp time(timestamp);

    //If a sync pulse arrived between this timestamp and the previous one seen.  This could skip Sync pulses entirely.
    if(time.lower < fPrevNs)
    {
      fNsOfLastSync = fPrevNs; //Estimate time in ns of previous Sync as previous timestamp seen
      fPrevNs = time.lower;
      return time.upper+(time.lower+fNsOfLastSync)*1e-9;
    }

    //If the UNIX timestamp has moved to a new value
    if(time.upper > fLastUNIXSecond)
    {
      fLastUNIXSecond = time.upper;
      fNsOfLastUNIXSecond = time.lower+fNsOfLastSync; 
    }

    return time.upper+(time.lower-fNsOfLastUNIXSecond)*1e-9; 
  }

  double Timer::elapsed()
  {
    const uint64_t lastTimestamp = ((uint64_t)fLastUNIXSecond << 32) | fPrevNs;
    return elapsed(lastTimestamp);
  }
} 
