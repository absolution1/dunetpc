////////////////////////////////////////////////////////////////////////
// \file XYZCalib.h
//
// \brief pure virtual base interface for (x,y,z) calibrations
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef CALIB_XYZCALIB_H
#define CALIB_XYZCALIB_H

namespace calib {
  
  class XYZCalib {
    
  public:
    
    XYZCalib(const XYZCalib &) = delete;
    XYZCalib(XYZCalib &&) = delete;
    XYZCalib& operator = (const XYZCalib &) = delete;
    XYZCalib& operator = (XYZCalib &&) = delete;
    virtual ~XYZCalib() = default;

    virtual double GetNormCorr(int plane) = 0;
    virtual double GetXCorr(int plane, double x) = 0;
    virtual double GetYZCorr(int plane, int side, double x, double y) = 0;

  protected:
    XYZCalib() = default;

  }; // class XYZCalib
} //namespace calib
#endif // CALIB_XYZCALIB_H
