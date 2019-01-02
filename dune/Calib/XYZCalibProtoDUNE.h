////////////////////////////////////////////////////////////////////////
// \file CalibProtoDUNE.h
//
// \brief header of class for accessing calibration data for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef CALIB_XYZCALIBPROTODUNE_H
#define CALIB_XYZCALIBPROTODUNE_H


// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

// dunetpc includes
#include "dune/Calib/XYZCalib.h"

namespace calib {

  typedef struct {
    double corr;
    double corr_err;
  } NormCorr_t;

  typedef struct {
    double x;
    double dx;
    double corr;
    double corr_err;
  } XCorr_t;

  bool operator<(const XCorr_t& a, const XCorr_t& b)
  {
    return (a.x < b.x);
  }

  typedef struct {
    double y;
    double dy;
    double z;
    double dz;
    double corr;
    double corr_err;
  } YZCorr_t;
  
  bool operator<(const YZCorr_t& a, const YZCorr_t& b)
  {
    double dy = a.y - b.y;
    if (dy < -1.e-5) return true;
    else if (dy > 1.e-5) return false;
    else return (a.z < b.z);
  }

  class XYZCalibProtoDUNE : public XYZCalib {
    
  public:

    XYZCalibProtoDUNE();
    XYZCalibProtoDUNE(fhicl::ParameterSet const& pset);
    XYZCalibProtoDUNE(XYZCalibProtoDUNE const&) = delete;
    virtual ~XYZCalibProtoDUNE() = default;
      
    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(uint64_t ts=0);
    
    virtual double GetNormCorr(int plane) override;
    virtual double GetXCorr(int plane, double x) override;
    virtual double GetYZCorr(int plane, int side, double y, double z) override;
    
    void SetIsMC(bool v) { fIsMC = v; }
    void SetNormCorr(int plane, double val);
    void SetXCorr(int plane, double x, double dx, double val);
    void SetYZCorr(int plane, int side, double y, double dx, double val);
    void SetUseCondb(bool v) { fUseCondbXYZCorr = v; }
    void SetXCorrFileName(std::string f) { fXCorrFileName=f; }
    void SetYZCorrFileName(std::string f) { fYZCorrFileName=f; }
    void SetNormCorrFileName(std::string f) { fNormCorrFileName=f; }

    void SetInterpolate(bool v) { fInterpolate = v; }

  protected:
      bool LoadNormCorr();
      bool LoadXCorr();
      bool LoadYZCorr();

  protected:
      bool fUseCondbXYZCorr;
      bool fNormCorrLoaded;
      bool fXCorrLoaded;
      bool fYZCorrLoaded;
      bool fIsMC;
      bool fInterpolate;
      uint64_t fCurrentTS;

      std::string fXCorrFileName;
      std::string fYZCorrFileName;
      std::string fNormCorrFileName;
      std::string fXCorrDBTag;
      std::string fYZCorrDBTag;
      std::string fNormCorrDBTag;

      std::map<int,NormCorr_t> fNormCorr;
      std::map<int,TH1F> fXCorrHist;
      std::map<int,TH2F> fYZCorrHist;

  }; // class XYZCalibProtoDUNE
} //namespace calib
#endif // CALIB_XYZCALIBPROTODUNE_H
