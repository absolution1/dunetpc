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
#include "cetlib_except/exception.h"

#include "nutools/IFDatabase/Util.h"
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
    fSCChanMap.clear();
    fSCInvChanMap.clear();
    fSCDataTable.reset();
  }
  
  //------------------------------------------------
  RunHistoryDUNE::~RunHistoryDUNE()
  {
  }

  //------------------------------------------------
  bool RunHistoryDUNE::LoadSCChanMap()
  {
    if (!fSCChanMap.empty()) return true;
    
    nutools::dbi::Table t;

    if (fDetName.empty()) {
      std::cout << __PRETTY_FUNCTION__ << ": Error in line " << __LINE__ << ".  Detector name is undefined." << std::endl;
      return false; //std::abort();
    }
    
    t.SetDetector(fDetName);
    t.SetTableName("daq_sc_chanmap");
    t.SetTableType(nutools::dbi::kConditionsTable);
    t.SetDataTypeMask(nutools::dbi::kDataOnly);

    int chanNameIdx = t.AddCol("chan_name","text");
    
    t.SetMinTSVld(1);
    t.SetMaxTSVld(1);

    t.SetVerbosity(100);
    if (! t.Load()) {
      std::cout << "Error in " << __PRETTY_FUNCTION__ << ", line " << __LINE__ << ".  Load from database failed." << std::endl;
      return false; //std::abort();
    }
    
    if (t.NRow() == 0) {
      std::cout << "Error in " << __PRETTY_FUNCTION__ << ", line " << __LINE__ << ".  Number of rows in table is 0.  This should never be the case!" << std::endl;
      return false;
    }

    nutools::dbi::Row* row;
    std::string chanName;
    uint64_t chan;
    for (int i=0; i<t.NRow(); ++i) {
      row = t.GetRow(i);      
      chan = row->Channel();
      row->Col(chanNameIdx).Get(chanName);
      fSCChanMap[chanName] = chan;
      fSCInvChanMap[chan] = chanName;
    }    

    return true;
  }
    
  //------------------------------------------------
  void RunHistoryDUNE::DumpSCData()
  {
    LoadSCData();
    LoadSCChanMap();
    
    fSCDataTable->FillChanRowMap();
    int rvIdx = fSCDataTable->GetColIndex("rvalue");
    float rv;
    
    std::vector<uint64_t> chanList = fSCDataTable->VldChannels();

    for (size_t ichan=0; ichan<chanList.size(); ++ichan) {
      std::vector<nutools::dbi::Row*> vldRow = fSCDataTable->GetVldRows(chanList[ichan]);
      std::cout << fSCInvChanMap[chanList[ichan]];
      
      for (size_t irow=0; irow<vldRow.size(); ++irow) {
	vldRow[irow]->Col(rvIdx).Get(rv);
	std::cout << ", (" << rv << "," << vldRow[irow]->VldTime()-fTStart << ")";
      }
      std::cout << std::endl;
    }
  }
  
  //------------------------------------------------
  void RunHistoryDUNE::DumpASICSettings()
  {
    LoadASICSettings();

    for (auto i : fASICSettingsMap) {
      int chan = i.first%100000;
      int asic = i.first/100000;
      std::cout << "ASIC " << asic << ", " << chan << ": " << i.second.gain << " mV/fC, " << i.second.shape << " us, " << i.second.base << " mV" << std::endl;
    }
  }

  //------------------------------------------------
  bool RunHistoryDUNE::LoadASICSettings()
  {
    if (fASICSettingsTable.get() != nullptr) return true;
    
    fASICSettingsTable.reset(new nutools::dbi::Table);
    fASICSettingsTable->SetDetector(fDetName);

    fASICSettingsTable->SetTableName("asic_settings_by_run");
    fASICSettingsTable->SetTableType(nutools::dbi::kConditionsTable);
    fASICSettingsTable->SetDataTypeMask(nutools::dbi::kDataOnly);
    
    int gainIdx = fASICSettingsTable->AddCol("gain","float");
    int shapeIdx = fASICSettingsTable->AddCol("shape","float");
    int baseIdx = fASICSettingsTable->AddCol("base","int");
    
    fASICSettingsTable->SetMinTSVld(fRun);
    fASICSettingsTable->SetMaxTSVld(fRun);

    fASICSettingsTable->SetVerbosity(100);
    if (! fASICSettingsTable->Load()) {
      std::cout << "Error in " << __PRETTY_FUNCTION__ << ", line " << __LINE__ << ".  Load from database failed." << std::endl;
      return false; 
    }

    if (!fASICSettingsTable->NRow()) {
      std::cout << "WARNING! " << __PRETTY_FUNCTION__ << ", no ASIC settings found for this run." << std::endl;
      return false;
    }
    
    fASICSettingsMap.clear();
    float g=0.;
    float s=0.;
    int b=0;
    for (int irow=0; irow<fASICSettingsTable->NRow(); ++irow) {
      nutools::dbi::Row* r = fASICSettingsTable->GetRow(irow);
      int chan = r->Channel();
      r->Col(gainIdx).Get(g);
      r->Col(shapeIdx).Get(s);
      r->Col(baseIdx).Get(b);
      fASICSettingsMap[chan] = ASICSetting(g,s,b);
    }
    return true;
  }  
  
  //------------------------------------------------
  bool RunHistoryDUNE::LoadSCData() 
  {
    if (fSCDataTable.get() != nullptr) return true;
    
    fSCDataTable.reset(new nutools::dbi::Table);
    fSCDataTable->SetDetector(fDetName);

    fSCDataTable->SetTableName("daq_slowcontrols");
    fSCDataTable->SetTableType(nutools::dbi::kConditionsTable);
    fSCDataTable->SetDataTypeMask(nutools::dbi::kNone);

    fSCDataTable->AddCol("rvalue","float");
    
    fSCDataTable->SetMinTSVld(fTStart);
    fSCDataTable->SetMaxTSVld(fTStop);

    fSCDataTable->SetVerbosity(100);
    if (! fSCDataTable->Load()) {
      std::cout << "Error in " << __PRETTY_FUNCTION__ << ", line " << __LINE__ << ".  Load from database failed." << std::endl;
      return false; //std::abort();
    }

    std::cout << "Read in " << fSCDataTable->NRow() << " rows of slow control data for run " << fRun << std::endl;
    
    return true;
  }
  
  //------------------------------------------------
  bool RunHistoryDUNE::Update(uint64_t run) 
  {
    if (run == 0) return false;

    fSCDataTable.reset();
    fASICSettingsTable.reset();
    
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

    //    t.SetVerbosity(100);
    t.SetTimeQueries(false);
    t.SetTimeParsing(false);
    t.Load();

    if (t.NRow() != 1)
      abort();

    std::string runTypeStr, compList;
    int runNum=0;
    
    nutools::dbi::Row* row = t.GetRow(0);
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
