////////////////////////////////////////////////////////////////////////
// \file TPCHVProtoDUNE.cxx
//
// \brief implementation of class for accessing TPC HV slow controls data for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "dune/SlowControlsServices/TPCHVServiceProtoDUNE.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
slowctrls::TPCHVServiceProtoDUNE::TPCHVServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new slowctrls::SlowControlsProtoDUNE(pset));

  //  reg.sPreBeginRun.watch(this, &TPCHVServiceProtoDUNE::preBeginRun);
}

//----------------------------------------------
//void slowctrls::TPCHVServiceProtoDUNE::preBeginRun(const art::Run& run)
//{
//  fProp->Update(run.id().run());
//}

//------------------------------------------------
void slowctrls::TPCHVServiceProtoDUNE::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(slowctrls::TPCHVServiceProtoDUNE, slowctrls::TPCHVService)
