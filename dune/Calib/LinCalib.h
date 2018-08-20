////////////////////////////////////////////////////////////////////////
// \file LinCalib.h
//
// \brief pure virtual base interface for linearity calibrations
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef CALIB_LINEARITYCALIB_H
#define CALIB_LINEARITYCALIB_H

namespace calib {
  
  class LinCalib {
    
  public:
    
    LinCalib(const LinCalib &) = delete;
    LinCalib(LinCalib &&) = delete;
    LinCalib& operator = (const LinCalib &) = delete;
    LinCalib& operator = (LinCalib &&) = delete;
    virtual ~LinCalib() = default;

  protected:
    LinCalib() = default;

  }; // class LinCalib
} //namespace calib
#endif // CALIB_LINCALIB_H
