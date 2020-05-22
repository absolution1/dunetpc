////////////////////////////////////////////////////////////////////////
// \file LifetimeCalibProtoDUNE_service.cc
//
// \brief implementation of class for storing/accessing lifetime corrections for ProtoDUNE
//
// \author wwu@fnal.gov
// \date May 19, 2020
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

#include "TTimeStamp.h"

// LArSoft includes
#include "dune/CalibServices/LifetimeCalibServiceProtoDUNE.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
calib::LifetimeCalibServiceProtoDUNE::LifetimeCalibServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new calib::LifetimeCalibProtoDUNE(pset));

  reg.sPreBeginRun.watch(this, &LifetimeCalibServiceProtoDUNE::preBeginRun);
}

//----------------------------------------------
void calib::LifetimeCalibServiceProtoDUNE::preBeginRun(const art::Run& run)
{
  art::Timestamp ts = run.beginTime();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  uint64_t  runtime = tts.AsDouble();
  
  // one can also consider using event time through "sPreProcessEvent" by define a "preProcessEvent" function, for exmaple.
  
  fProp->Update(runtime);

  //fProp->Update(run.id().run());
}

//------------------------------------------------
void calib::LifetimeCalibServiceProtoDUNE::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::LifetimeCalibServiceProtoDUNE, calib::LifetimeCalibService)
