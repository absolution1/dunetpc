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
  fAString("")
{ }



//.....................................................................
TriggerTestService::TriggerTestService(const fhicl::ParameterSet& pset, art::ActivityRegistry&) :
  fIndex(0),
  fAString("")
{
  fAString = pset.get<std::string>("TriggerLabel");
}



//.....................................................................
bool TriggerTestService::TPCtrigger(art::Event const & evt)
{

  std::cout << "\n\nIn the TPCtrigger function for event "
	    << evt.id().event()
	    << "\n\nfIndex = " << fIndex
	    << "\nString = " << fAString
	    << "\n\n\n";

  fIndex++;

  return true;
}



//.....................................................................
bool TriggerTestService::PDtrigger(art::Event const & evt)
{

  std::cout << "\n\nIn the PDtrigger function for event "
	    << evt.id().event()
	    << "\n\n\n";

  return false;

}


  

DEFINE_ART_SERVICE_INTERFACE_IMPL(TriggerTestService, BaseTriggerService)
