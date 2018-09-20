//File: Timer.h
//Brief: A CRT::Timer computes the current and/or elapsed time since its' construction.  CRT timestamps have 
//       a UNIX timestamp that is not synchronized with the clock timestamp.  So, a CRT::Timer uses the times 
//       that it is asked to convert to infer when the lower 32 bits of the timestamp have reset. This means 
//       that you will get less accurate time values (on the order of single seconds) if you don't "show" 
//       a CRT::Timer all of the timestamps in an event.  
//Author: Andrew Olivier aolivier@ur.rochester.edu 

#ifndef CRT_TIMER_H
#define CRT_TIMER_H

//c++ includes
#include <cstdint>

namespace CRT
{
  class Trigger;
  
  //Represent a CRT timestamp in two parts
  struct Timestamp
  {
    //Create a Timestamp from the format it is saved in in a CRT::Trigger
    Timestamp(const uint64_t timestamp)
    {
      lower = timestamp; //Let the compiler strip off the upper 32 bits
      upper = timestamp >> 32; //Push the lower 32 bits off the end of the number, then let the compiler 
                               //strip off the remaining 32 bits
    }

    //Two parts of a CRT timestamp: 
    uint32_t lower; //Time in ns since last Sync
    uint32_t upper; //UNIX timestamp in seconds
  };

  class Timer
  {
    public:
      //TODO: using statment to make sure timestamp type is whatever CRT::Trigger uses as a timestamp? 
      //Construct a CRT::Timer from the timestamp relative to which it will report elapsed time
      Timer(const CRT::Trigger& trigger);
      Timer(const uint64_t timestamp); 

      //Get the time in seconds represented by a trigger or timestamp
      double seconds(const CRT::Trigger& trigger);
      double seconds(const uint64_t timestamp);

      //Get the elapsed time since the timestamp with which this Timer was constructed
      double elapsed(const CRT::Trigger& trigger);
      double elapsed(const uint64_t timestamp);

      //Get the elapsed time between the timestamp with which this Timer was constructed and the last timestamp it saw
      double elapsed();

    private:
      //State that Timer uses to estimate the real time of a CRT timestamp
      uint32_t fLastUNIXSecond; //Last UNIX second seen by Timer 
      uint32_t fNsOfLastUNIXSecond; //Inferred ns since the UNIX timestamp last went to a new second
      uint32_t fNsOfLastSync; //Inferred number of ns when the last time a sync signal was sent to the CRT boards and they 
                                    //reset their 
      uint32_t fPrevNs; //Previous ns component of time 
      const double fFirstTime; //Time of the timestamp with which this CRT::Timer was constructed

      //Turn a 64 bit timestamp into seconds
      //TODO: Move this to the Timestamp struct below as a cast to double?
      double GetSeconds(const uint64_t timestamp);
  };
}

#endif //CRT_TIMER_H
