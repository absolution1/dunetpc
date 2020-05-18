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
  fUseCondbLifetime = pset.get<bool>("UseCondbLifetime");
  fInterpolate     = pset.get<bool>("Interpolate");
  fLifetimeFileName   = pset.get<std::string>("LifetimeFileName");
  fLifetimeDBTag      = pset.get<std::string>("LifetimeDBTag");
  return true;
}

//------------------------------------------------
bool calib::LifetimeCalibProtoDUNE::Update(uint64_t ts) 
{

  if (fLifetimeLoaded && ts != fCurrentTS) {
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

  return 0;
}

//------------------------------------------------
double calib::LifetimeCalibProtoDUNE::GetLifetimeHigh() 
{
  if (!fLifetimeLoaded) this->LoadLifetime();

  return 0;
}

//------------------------------------------------
double calib::LifetimeCalibProtoDUNE::GetLifetimeLow() 
{
  if (!fLifetimeLoaded) this->LoadLifetime();

  return 0;
}

//------------------------------------------------
bool calib::LifetimeCalibProtoDUNE::LoadLifetime()
{
  if (!fUseCondbLifetime) return true;

  if (fLifetimeLoaded) return true;

  fLifetimeLoaded = true;
  return true;
}



