////////////////////////////////////////////////////////////////////////
//
// basic timer
//
//
////////////////////////////////////////////////////////////////////////
#ifndef __TIMER_H__
#define __TIMER_H__

#include <sys/time.h>
#include <cstdio>
#include <cstdlib>

class Timer
{
 private:
  timeval startTime;
  timeval lastTime;
  
  double duration(timeval &tstart, timeval &tend)
  {
    long sec  = tend.tv_sec  - tstart.tv_sec;
    long usec = tend.tv_usec - tstart.tv_usec;
    double tdiff = sec + usec/1.0E+6;
    return tdiff;
  }

  Timer() { start(); }
    
  Timer(const Timer&);
  Timer& operator=(const Timer&);
  ~Timer(){;}

 public:
  static Timer& GetTimer()
    {
      static Timer inst;
      return inst;
    }

  void start()
  {
    gettimeofday(&startTime, NULL);
    lastTime = startTime;
  }
  
  double splittime(bool _save, bool _echo = true)
  {
    timeval curTime;
    gettimeofday(&curTime, NULL);
    
    double tdiff = duration(lastTime, curTime);
    if(_save) lastTime = curTime;
    if(_echo) printf("Time from last %.6f s\n", tdiff);
    return tdiff;
  }
  
  double stop(bool _echo = true)
  {
    timeval curTime;
    gettimeofday(&curTime, NULL);
    
    double tdiff = duration(startTime, curTime);
    if(_echo) printf("Time from start %.6f s\n", tdiff);
    return tdiff;
  }
};


#endif
