////////////////////////////////////////////////////////////////////////
// \file XYZCalibProtoDUNE.cxx
//
// \brief implementation of class for storing/accessing (x,y,z) corrections for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "dune/CalibServices/XYZCalibServiceProtoDUNE.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
calib::XYZCalibServiceProtoDUNE::XYZCalibServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new calib::XYZCalibProtoDUNE(pset));

  reg.sPreBeginRun.watch(this, &XYZCalibServiceProtoDUNE::preBeginRun);
}

//----------------------------------------------
void calib::XYZCalibServiceProtoDUNE::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
void calib::XYZCalibServiceProtoDUNE::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::XYZCalibServiceProtoDUNE, calib::XYZCalibService)
