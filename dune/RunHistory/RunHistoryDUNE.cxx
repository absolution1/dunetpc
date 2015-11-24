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
#include "IFDatabase/Util.h"
#include <boost/tokenizer.hpp>
#include <ctime>

namespace dune {
  //-----------------------------------------------
  RunHistoryDUNE::RunHistoryDUNE(int detid, int run)
  {
    if (detid > 0 && detid < kNDUNEDetectors && run > 0) {
      fDetId = detid;
      fRun=run;
      switch(fDetId) {
      case(k35t):
	fDetName = "dune35t";
	break;
      case(kProtoDUNE):
	fDetName = "protoDUNE";
	break;	
      case(kFarDet):
	fDetName = "FarDet";
	break;
      case(kNearDet):
	fDetName = "NearDet";
	break;
      default:
	fDetName = "";
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
    nutools::dbi::Table t;
    
    t.SetDetector(fDetName);
    t.SetTableName(tableName);
    t.SetTableType(nutools::dbi::kGenericTable);

    int runIdx = t.AddCol("run", "integer");
    int cfgLabelIdx = t.AddCol("configuration_label", "integer");
    int runTypeIdx = t.AddCol("run_type", "text");
    int compListIdx = t.AddCol("component_list", "text");
    int startIdx = t.AddCol("start", "timestamp");
    int stopIdx = t.AddCol("stop", "timestamp");

    t.SetValidityRange("run",run);

    t.SetVerbosity(100);
    t.Load();

    if (t.NRow() != 1)
      abort();

    std::string runTypeStr, compList;
    int runNum=0;
    
    nutools::dbi::Row* row = t.GetRow(0);
    std::cout << *row << std::endl;
    row->Col(runIdx).Get(runNum);
    row->Col(cfgLabelIdx).Get(fCfgLabel);
    row->Col(compListIdx).Get(compList);
    row->Col(runTypeIdx).Get(runTypeStr);
    row->Col(startIdx).Get(fTStartStr);
    row->Col(stopIdx).Get(fTStopStr);

    if (runNum != int(run)) abort();

    boost::tokenizer< boost::escaped_list_separator<char> > tok(compList);
    fComponents.assign(tok.begin(),tok.end());

    fTStart = fTStop = 0;
    time_t tval;
    nutools::dbi::Util::TimeAsStringToTime_t(fTStartStr,tval);
    fTStart = tval;
    if (fTStopStr != "" && fTStopStr != "None") {
      nutools::dbi::Util::TimeAsStringToTime_t(fTStopStr,tval);
      fTStop = tval;
    }
    
    std::cout << "Run " << runNum
      	      << "\nCfg: " << fCfgLabel
      	      << "\nComponents: ";
    for (size_t i=0; i<fComponents.size(); ++i)
      std::cout << fComponents[i] << " ";

    
    std::cout << "\nrunType: " << runTypeStr 
      	      << "\nStart time: " << fTStartStr << "(" << fTStart << ")"
      	      << "\nStop time: " << fTStopStr << "(" << fTStop << ")"
	      << std::endl;
    if (fTStop > fTStart)
      std::cout << "Duration: " << (fTStop - fTStart) << std::endl;
    
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
