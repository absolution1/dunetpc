////////////////////////////////////////////////////////////////////////
// \file SlowControls.h
//
// \brief pure virtual base interface for slow controls data
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef SLOW_CONTROLS_H
#define SLOW_CONTROLS_H

#include <string>

namespace slowctrls {
  
  class SlowControls {
    
  public:
    
    SlowControls(const SlowControls &) = delete;
    SlowControls(SlowControls &&) = delete;
    SlowControls& operator = (const SlowControls &) = delete;
    SlowControls& operator = (SlowControls &&) = delete;
    virtual ~SlowControls() = default;

    virtual double GetValue(std::string& chan, float t) = 0;

  protected:
    SlowControls() = default;

  }; // class SlowControls
} //namespace slowctrls
#endif // SLOW_CONTROLS_H
