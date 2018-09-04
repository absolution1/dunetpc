////////////////////////////////////////////////////////////////////////
// \file LinCalibProtoDUNE.h
//
// \brief header of class for accessing linearity calibration 
//        data for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef CALIB_LINCALIBPROTODUNE_H
#define CALIB_LINCALIBPROTODUNE_H


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
#include "dune/Calib/LinCalib.h"

namespace calib {

  typedef struct {
    int status;
    float gain;
    float offset;
    float shape;
    float chi2;
    int adc_low;
    int adc_high;
    float nl[20];
  } LinConsts_t;

  class LinCalibProtoDUNE : public LinCalib {
    
  public:

    LinCalibProtoDUNE();
    LinCalibProtoDUNE(fhicl::ParameterSet const& pset);
    LinCalibProtoDUNE(LinCalibProtoDUNE const&) = delete;
    virtual ~LinCalibProtoDUNE() = default;
      
    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(uint64_t ts=0);
    
    LinConsts_t GetLinConsts(int chanId); 
    
    void SetIsMC(bool v) { fIsMC = v; }

    void SetUseCondb(bool v) { fUseCondb = v; }
    void SetCSVFileName(std::string f) { fCSVFileName=f; }
    void SetTag(std::string t) { fDBTag = t; }

  protected:
      bool LoadConsts();

      bool fUseCondb;
      bool fConstsLoaded;
      bool fIsMC;
      uint64_t fCurrentTS;

      std::string fCSVFileName;
      std::string fDBTag;

      std::map<int,LinConsts_t> fLinConsts;

  }; // class LinCalibProtoDUNE
} //namespace calib
#endif // CALIB_LINCALIBPROTODUNE_H
