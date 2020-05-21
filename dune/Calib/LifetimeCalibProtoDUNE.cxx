////////////////////////////////////////////////////////////////////////
// \file LifetimeCalibProtoDUNE.cxx
//
// \brief implementation of class for accessing (x,y,z) calibration data for ProtoDUNE
//
// \author tjyang@fnal.gov
//         wwu@fnal.gov
//
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"
#include <iomanip> 

// LArSoft includes
#include "dune/Calib/LifetimeCalibProtoDUNE.h"

// nutools includes
#include "nuevdb/IFDatabase/Table.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
calib::LifetimeCalibProtoDUNE::LifetimeCalibProtoDUNE()
{
  fIsMC = true;
  fLifetimeLoaded = false;
  fInterpolate = false;
  fCurrentTS = 0;
  fLifetimeFileName="";
  fLifetimeDBTag="";
}


//-----------------------------------------------
calib::LifetimeCalibProtoDUNE::LifetimeCalibProtoDUNE(
  fhicl::ParameterSet const& pset
)
{
  fIsMC = true;
  fLifetimeLoaded = false;
  fInterpolate = false;
  fCurrentTS = 0;
  fLifetimeFileName="";
  fLifetimeDBTag="";
  Configure(pset);
}

//------------------------------------------------
bool calib::LifetimeCalibProtoDUNE::Configure(fhicl::ParameterSet const& pset)
{  
  fIsMC             = pset.get<bool>("IsMC");
  fUseCondbLifetime = pset.get<bool>("UseCondbLifetime");
  fInterpolate      = pset.get<bool>("Interpolate");
  fLifetimeFileName = pset.get<std::string>("LifetimeFileName");
  fLifetimeDBTag    = pset.get<std::string>("LifetimeDBTag");
  return true;
}

//------------------------------------------------
bool calib::LifetimeCalibProtoDUNE::Update(uint64_t ts) 
{
  
  if (fLifetimeLoaded && ts != fCurrentTS) {
    fLifetimePurMon.clear();
    fLifetimeLoaded = false;
  }

  fCurrentTS = ts;
  // all done! 

  return true;
}

//------------------------------------------------
double calib::LifetimeCalibProtoDUNE::GetLifetime() 
{
  if (!fLifetimeLoaded) this->LoadLifetime();

  if (fLifetimePurMon.find(1) == fLifetimePurMon.end()) {
    mf::LogError("LifetimeCalibProtoDUNE") << "channel not found!";
    return 0.;
  }

  return fLifetimePurMon[1].center;
}

//------------------------------------------------
double calib::LifetimeCalibProtoDUNE::GetLifetimeHigh() 
{
  if (!fLifetimeLoaded) this->LoadLifetime();

  if (fLifetimePurMon.find(1) == fLifetimePurMon.end()) {
    mf::LogError("LifetimeCalibProtoDUNE") << "channel not found!";
    return 0.;
  }

  return fLifetimePurMon[1].high;
}

//------------------------------------------------
double calib::LifetimeCalibProtoDUNE::GetLifetimeLow() 
{
  if (!fLifetimeLoaded) this->LoadLifetime();

  if (fLifetimePurMon.find(1) == fLifetimePurMon.end()) {
    mf::LogError("LifetimeCalibProtoDUNE") << "channel not found!";
    return 0.;
  }

  return fLifetimePurMon[1].low;
}

//------------------------------------------------
bool calib::LifetimeCalibProtoDUNE::LoadLifetime()
{
  if (!fUseCondbLifetime) return true;

  if (fLifetimeLoaded) return true;
  
  nutools::dbi::Table LifetimePurMonTable;
  LifetimePurMonTable.SetDetector("pdunesp");
  LifetimePurMonTable.SetTableName("lifetime_purmon");
  LifetimePurMonTable.SetTableType(nutools::dbi::kConditionsTable);
  LifetimePurMonTable.SetDataTypeMask(nutools::dbi::kDataOnly);
  if (fIsMC) LifetimePurMonTable.SetDataTypeMask(nutools::dbi::kMCOnly);

  int centerIdx = LifetimePurMonTable.AddCol("center", "double");
  int lowIdx = LifetimePurMonTable.AddCol("low", "double");
  int highIdx = LifetimePurMonTable.AddCol("high", "double");

  LifetimePurMonTable.SetMinTSVld(fCurrentTS);
  LifetimePurMonTable.SetMaxTSVld(fCurrentTS);
  LifetimePurMonTable.SetTag(fLifetimeDBTag);
  
  LifetimePurMonTable.SetVerbosity(100);

  bool readOk = false;
  if (!fLifetimeFileName.empty())
    readOk = LifetimePurMonTable.LoadFromCSV(fLifetimeFileName);
  else
    readOk = LifetimePurMonTable.Load();
  
  //std::cout << "www...  LifetimePurMonTable.GetMinTSVld(): "<< std::setprecision(12) << LifetimePurMonTable.GetMinTSVld() << std::endl;
  //std::cout << "www...  LifetimePurMonTable.GetMaxTSVld(): "<< std::setprecision(12)<< LifetimePurMonTable.GetMaxTSVld() << std::endl;
  //std::cout << "www...fCurrentTS: " << fCurrentTS << std::endl;

  if (! readOk) {
   mf::LogError("LifetimeCalibProtoDUNE") << "Load from lifetime calib database table failed.";

   return false; //std::abort();
  }

  if (LifetimePurMonTable.NRow() == 0) {
    mf::LogError("LifetimeCalibProtoDUNE") << "Number of rows in lifetime calib table is 0.  This should never be the case!";
    return false;
  }

  //std::cout << "www... LifetimePurMonTable.NRow(): " << LifetimePurMonTable.NRow() << std::endl;

  nutools::dbi::Row* row;
  uint64_t chan;
  for (int i=0; i<LifetimePurMonTable.NRow(); ++i) {
    LifetimePurMon_t lifetime;
    row = LifetimePurMonTable.GetRow(i);
    chan = row->Channel(); // One channel only, start with 1

    if (chan != 1) {
      mf::LogError("LifetimeCalibProtoDUNE") << "Channel numuber in lifetime calib table is not  1.  This should never be the case!";
      return false;
    }
    row->Col(centerIdx).Get(lifetime.center);
    row->Col(lowIdx).Get(lifetime.low);
    row->Col(highIdx).Get(lifetime.high);
    fLifetimePurMon[chan] = lifetime;
  }


  fLifetimeLoaded = true;
  return true;
}



