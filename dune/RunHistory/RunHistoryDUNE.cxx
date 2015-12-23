////////////////////////////////////////////////////////////////////////
//
//  RunHistoryDUNE
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// C++ language includes
#include <iostream>

// LArSoft includes
#include "RunHistoryDUNE.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

#include "IFDatabase/Table.h"

namespace dune {
  //-----------------------------------------------
  RunHistoryDUNE::RunHistoryDUNE(int detid, int run)
  {
    if (detid > 0 && detid < kNDUNEDetectors && run > 0) {
      fDetId = detid;
      fRun=run;
      switch(fDetId) {
      case(k35t):
	detName = "dune35t";
	break;
      case(kProtoDUNE):
	detName = "protoDUNE";
	break;	
      case(kFarDet):
	detName = "FarDet";
	break;
      case(kNearDet):
	detName = "NearDet";
	break;
      default:
	detName = "";
	break;
      }
    }
    else 
      abort();
    
    if (!Update(run)) abort();
  }
  
  //------------------------------------------------
  RunHistoryDUNE::~RunHistoryDUNE()
  {
  }

  //------------------------------------------------
  bool RunHistoryDUNE::Update(uint64_t run) 
  {
    if (run == 0) return false;

    std::string tableName = "run_summary";
    auto * t = new nutools::dbi::Table(detName,tableName,nutools::dbi::kGenericTable);
    int cfgLabelIdx = t->AddCol("configuration_label", "integer");
    int runTypeIdx = t->AddCol("run_type", "text");
    int compListIdx = t->AddCol("component_list", "text");
    int startIdx = t->AddCol("start", "timestamp");
    int stopIdx = t->AddCol("stop", "timestamp");

    t->SetValidityRange("run",run);

    t->Load();

    if (t->NRow() != 1)
      abort();

    nutools::dbi::Row* row = t->GetRow(0);
    
    
    return true;
    
  }

  //------------------------------------------------
  std::string RunHistoryDUNE::RunTypeAsString() const
  {
    switch(fRunType) {
    case(detinfo::kProductionRun):
      return std::string("Production");
    case(detinfo::kCommissioningRun):
      return std::string("Commissioning");
    case(detinfo::kTestRun):
      return std::string("Test");
    case(detinfo::kPedestalRun):
      return std::string("Pedestal");
    case(detinfo::kCalibrationRun):
      return std::string("Calibration");
    case(detinfo::kUnknownRunType):
    default:
      return std::string("Uknown");
    }
  }
}
