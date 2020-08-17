////////////////////////////////////////////////////////////////////////
// \file LifetimeCalib.h
//
// \brief pure virtual base interface for lifetime calibration
//
// \author tjyang@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef CALIB_LIFETIMECALIB_H
#define CALIB_LIFETIMECALIB_H

namespace calib {
  
  class LifetimeCalib {
    
  public:
    
    LifetimeCalib(const LifetimeCalib &) = delete;
    LifetimeCalib(LifetimeCalib &&) = delete;
    LifetimeCalib& operator = (const LifetimeCalib &) = delete;
    LifetimeCalib& operator = (LifetimeCalib &&) = delete;
    virtual ~LifetimeCalib() = default;

    virtual double GetLifetime() = 0;
    virtual double GetLifetimeLow() = 0;
    virtual double GetLifetimeHigh() = 0;

  protected:
    LifetimeCalib() = default;

  }; // class LifetimeCalib
} //namespace calib
#endif // CALIB_LIFETIMECALIB_H
