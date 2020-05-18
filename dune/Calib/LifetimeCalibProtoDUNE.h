////////////////////////////////////////////////////////////////////////
// \file CalibProtoDUNE.h
//
// \brief header of class for accessing calibration data for ProtoDUNE
//
// \author tjyang@fnal.gov
//         wwu@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef CALIB_LIFETIMECALIBPROTODUNE_H
#define CALIB_LIFETIMECALIBPROTODUNE_H


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
#include "dune/Calib/LifetimeCalib.h"

namespace calib {

  class LifetimeCalibProtoDUNE : public LifetimeCalib {
    
  public:

    LifetimeCalibProtoDUNE();
    LifetimeCalibProtoDUNE(fhicl::ParameterSet const& pset);
    LifetimeCalibProtoDUNE(LifetimeCalibProtoDUNE const&) = delete;
    virtual ~LifetimeCalibProtoDUNE() = default;
      
    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(uint64_t ts=0);
    
    virtual double GetLifetime() override;
    virtual double GetLifetimeLow() override;
    virtual double GetLifetimeHigh() override;
    
    void SetIsMC(bool v) { fIsMC = v; }
    void SetLifetime(double val);
    void SetLifetimeLow(double val);
    void SetLifetimeHigh(double val);
    void SetUseCondb(bool v) { fUseCondbLifetime = v; }
    void SetLifetimeFileName(std::string f) { fLifetimeFileName=f; }
    void SetInterpolate(bool v) { fInterpolate = v; }

  protected:
      bool LoadLifetime();

      bool fUseCondbLifetime;
      bool fLifetimeLoaded;
      bool fIsMC;
      bool fInterpolate;
      uint64_t fCurrentTS;

      std::string fLifetimeFileName;
      std::string fLifetimeDBTag;

  }; // class LifetimeCalibProtoDUNE
} //namespace calib
#endif // CALIB_LIFETIMECALIBPROTODUNE_H
