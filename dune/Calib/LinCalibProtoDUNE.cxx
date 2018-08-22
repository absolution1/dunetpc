////////////////////////////////////////////////////////////////////////
// \file LinCalibProtoDUNE.cxx
//
// \brief implementation of class for accessing linearity calibration 
//        constants for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

// LArSoft includes
#include "dune/Calib/LinCalibProtoDUNE.h"

// nutools includes
#include "nutools/IFDatabase/Table.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
calib::LinCalibProtoDUNE::LinCalibProtoDUNE()
{
  fIsMC = true;
  fConstsLoaded = false;
  fCurrentTS = 0;
  fCSVFileName="";
  fDBTag="";
}


//-----------------------------------------------
calib::LinCalibProtoDUNE::LinCalibProtoDUNE(
  fhicl::ParameterSet const& pset
)
{
  fIsMC = true;
  fConstsLoaded = false;
  fCurrentTS = 0;
  fCSVFileName="";
  fDBTag="";
}

//------------------------------------------------
bool calib::LinCalibProtoDUNE::Configure(fhicl::ParameterSet const& pset)
{  
  fUseCondb      = pset.get<bool>("UseCondb");
  fCSVFileName   = pset.get<std::string>("CSVFileName");
  fDBTag         = pset.get<std::string>("DBTag");

  return true;
}

//------------------------------------------------
bool calib::LinCalibProtoDUNE::Update(uint64_t ts) 
{

  if (fConstsLoaded && ts != fCurrentTS) {
    fLinConsts.clear();
    fConstsLoaded = false;
  }

  fCurrentTS = ts;
  // all done! 

  return true;
}

//------------------------------------------------
calib::LinConsts_t calib::LinCalibProtoDUNE::GetLinConsts(int chanId) 
{
  if (!fConstsLoaded) this->LoadConsts();

  if (fLinConsts.find(chanId) == fLinConsts.end()) {
    mf::LogError("LinCalibProtoDUNE") << "Channel " << chanId 
				      << "not found!";
    std::abort();
  }

  return fLinConsts[chanId];
}

//------------------------------------------------
bool calib::LinCalibProtoDUNE::LoadConsts()
{
  if (!fUseCondb) return true;

  if (fConstsLoaded) return true;

  nutools::dbi::Table t;
  
  t.SetDetector("pdunesp");
  t.SetTableName("linconsts");
  t.SetTableType(nutools::dbi::kConditionsTable);
  t.SetDataTypeMask(nutools::dbi::kDataOnly);
  if (fIsMC)
    t.SetDataTypeMask(nutools::dbi::kMCOnly);
  
  int statusIdx = t.AddCol("status","int");
  int gainIdx   = t.AddCol("gain","float");
  int offsetIdx = t.AddCol("offset","float");
  int shapeIdx  = t.AddCol("shape","float");
  int chi2Idx   = t.AddCol("chi2","float");
  int adcLowIdx = t.AddCol("adc_low","int");
  int adcHiIdx  = t.AddCol("adc_high","int");
  int nlIdx[20];
  char buff[64];
  for (int i=0; i<20; ++i) {
    sprintf(buff,"nl%d",i);
    nlIdx[i] = t.AddCol(buff,"float");
  }
  
  t.SetMinTSVld(fCurrentTS);
  t.SetMaxTSVld(fCurrentTS);
  t.SetTag(fDBTag);

  t.SetVerbosity(100);

  bool readOk = false;
  if (!fCSVFileName.empty()) 
    readOk = t.LoadFromCSV(fCSVFileName);
  else
    readOk = t.Load();

  if (! readOk) {
    mf::LogError("LinCalibProtoDUNE") << "Load from norm linconsts database table failed.";
    
    return false; //std::abort();

  }
  
  if (t.NRow() == 0) {
    mf::LogError("LinCalibProtoDUNE") << "Number of rows in linconsts calib table is 0.  This should never be the case!";
    return false;
  }
  
  nutools::dbi::Row* row;
  uint64_t chan;
  for (int i=0; i<t.NRow(); ++i) {
    LinConsts_t c;
    row = t.GetRow(i);      
    chan = row->Channel();
    row->Col(statusIdx).Get(c.status);
    row->Col(gainIdx).Get(c.gain);
    row->Col(offsetIdx).Get(c.offset);
    row->Col(shapeIdx).Get(c.shape);
    row->Col(chi2Idx).Get(c.chi2);
    row->Col(adcLowIdx).Get(c.adc_low);
    row->Col(adcHiIdx).Get(c.adc_high);
    for (int j=0; j<20; ++j) {
      row->Col(nlIdx[j]).Get(c.nl[j]);
    }
    
    fLinConsts[chan] = c;
  }    

  fConstsLoaded = true;
  return true;
}


