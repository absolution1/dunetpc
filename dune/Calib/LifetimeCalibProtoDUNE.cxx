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
  
  if (!fInterpolate) {
    LifetimePurMonTable.SetMinTSVld(fCurrentTS); // only load one previous lifetime
    LifetimePurMonTable.SetMaxTSVld(fCurrentTS);
  }
  else {
    // load lifetime IOV: runtime +- 1 days for interpolation
    LifetimePurMonTable.SetMinTSVld(fCurrentTS-1*86400.);
    LifetimePurMonTable.SetMaxTSVld(fCurrentTS+1*86400.);
  }
  LifetimePurMonTable.SetTag(fLifetimeDBTag);
  
  LifetimePurMonTable.SetVerbosity(100);

  bool readOk = false;
  if (!fLifetimeFileName.empty())
    readOk = LifetimePurMonTable.LoadFromCSV(fLifetimeFileName);
  else
    readOk = LifetimePurMonTable.Load();
  
  if (! readOk) {
   mf::LogError("LifetimeCalibProtoDUNE") << "Load from lifetime calib database table failed.";
   throw cet::exception("LifetimeCalibProtoDUNE") << "Failed to query lifetime database. Please check URL: https://dbdata0vm.fnal.gov:9443/dune_con_prod/app" << "\n";

   return false; //std::abort();
  }

  if (LifetimePurMonTable.NRow() == 0) {
    mf::LogError("LifetimeCalibProtoDUNE") << "Number of rows in lifetime calib table is 0.  This should never be the case!";
    return false;
  }

  std::vector<double> loaded_tv;
  std::vector<double> loaded_center;
  std::vector<double> loaded_low;
  std::vector<double> loaded_high;
  loaded_tv.clear();
  loaded_center.clear();
  loaded_low.clear();
  loaded_high.clear();

  if (LifetimePurMonTable.NRow() == 1) fInterpolate = false;

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

    loaded_tv.push_back(row->VldTime());
    loaded_center.push_back(lifetime.center);
    loaded_low.push_back(lifetime.low);
    loaded_high.push_back(lifetime.high);
    
    if (!fInterpolate) fLifetimePurMon[chan] = lifetime;
  }

  if (fInterpolate && loaded_tv.size()>1) {

    LifetimePurMon_t interpolate_lifetime;

    // find the previous and following lifetime for fCurrentTS
    int t0_idx = -1; 
    int t1_idx = -1;
    double temp_0 = 1e10;
    double temp_1 = 1e10; 
    for (size_t j=0; j<loaded_tv.size(); j++) {
      if ( fCurrentTS - loaded_tv[j] >= 0 && fCurrentTS - loaded_tv[j] < temp_0 ) {
        temp_0 = fCurrentTS- loaded_tv[j];
        t0_idx = j;
      }

      if ( loaded_tv[j]-fCurrentTS >= 0 && loaded_tv[j]-fCurrentTS < temp_1 ) {
        temp_1 = loaded_tv[j]-fCurrentTS;
        t1_idx = j;
      }
    }
    
    if (t0_idx == -1) { // only lifetime measurement after fCurrentTS during selected IOV
      // use the following one only
      interpolate_lifetime.center = loaded_center[t1_idx];
      interpolate_lifetime.low = loaded_low[t1_idx];
      interpolate_lifetime.high = loaded_high[t1_idx];
    }
    if (t1_idx == -1) { // only lifetime measurement before fCurrentTS during selected IOV
      // use the previous one only
      interpolate_lifetime.center = loaded_center[t0_idx];
      interpolate_lifetime.low = loaded_low[t0_idx];
      interpolate_lifetime.high = loaded_high[t0_idx];
    }
    if (t0_idx != -1 && t1_idx != -1) {
      if (t0_idx == t1_idx) { // ideally, this is not possible
        interpolate_lifetime.center = loaded_center[t0_idx];
        interpolate_lifetime.low = loaded_low[t0_idx];
        interpolate_lifetime.high = loaded_high[t0_idx];
      }
      else {
        // do linear interpolation using previous one and following one
        // refer: https://en.wikipedia.org/wiki/Linear_interpolation
        interpolate_lifetime.center = ( loaded_center[t0_idx]*(loaded_tv[t1_idx]-fCurrentTS) + loaded_center[t1_idx]*(fCurrentTS-loaded_tv[t0_idx]) ) / (loaded_tv[t1_idx]-loaded_tv[t0_idx]);
        interpolate_lifetime.low = ( loaded_low[t0_idx]*(loaded_tv[t1_idx]-fCurrentTS) + loaded_low[t1_idx]*(fCurrentTS-loaded_tv[t0_idx]) ) / (loaded_tv[t1_idx]-loaded_tv[t0_idx]);
        interpolate_lifetime.high = ( loaded_high[t0_idx]*(loaded_tv[t1_idx]-fCurrentTS) + loaded_high[t1_idx]*(fCurrentTS-loaded_tv[t0_idx]) ) / (loaded_tv[t1_idx]-loaded_tv[t0_idx]);
      }
    }

    fLifetimePurMon[chan] = interpolate_lifetime;
  } // end if fInterpolate

  fLifetimeLoaded = true;
  return true;
}



