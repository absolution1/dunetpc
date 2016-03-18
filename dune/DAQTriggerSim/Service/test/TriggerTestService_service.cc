// TriggerTestService_service.cc
//
// Michael Baird
// March 2016
//
// Test service class for the DAQ trigger service framework.
//

// DUNETPC specific includes
#include "dune/DAQTriggerSim/Service/test/TriggerTestService.h"

// Framework includes
#include "fhiclcpp/ParameterSet.h"

// C++ includes
#include <iostream>



//.....................................................................
TriggerTestService::TriggerTestService() :
  fIndex(0),
  fAString("no value...")
{ }



//.....................................................................
TriggerTestService::TriggerTestService(const fhicl::ParameterSet& pset, art::ActivityRegistry&) :
  fIndex(0),
  fAString("no value...")
{
  fAString = pset.get<std::string>("AString");
}



//.....................................................................
bool TriggerTestService::Trigger(art::Event const & evt)
{
  
  std::cout << "\n\nIn the Trigger function for event "
	    << evt.id().event()
	    << "\n\nfIndex = " << fIndex
	    << "\nString = " << fAString
	    << "\n\n\n";
  
  fIndex++;
  
  return true;
}



DEFINE_ART_SERVICE_INTERFACE_IMPL(TriggerTestService, BaseTriggerService)
