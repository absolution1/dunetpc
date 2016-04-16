// TemplateTriggerService_service.cc
//
// Michael Baird
// March 2016
//
// Template trigger service class for the DAQ trigger service framework.
//

// DUNETPC specific includes
#include "dune/DAQTriggerSim/Service/TemplateTriggerService.h"

// Framework includes
#include "fhiclcpp/ParameterSet.h"

// C++ includes
#include <iostream>



//.....................................................................
TemplateTriggerService::TemplateTriggerService() :
  fIndex(0),
  fAString("no value...")
{ }



//.....................................................................
TemplateTriggerService::TemplateTriggerService(const fhicl::ParameterSet& pset, art::ActivityRegistry&) :
  fIndex(0),
  fAString("no value...")
{
  fAString = pset.get<std::string>("AString");
}



//.....................................................................
bool TemplateTriggerService::Trigger(art::Event const & evt)
{

  //
  // Into this function will go the guts of the trigger algorithm.
  //
  
  std::cout << "\n\nIn the Trigger function for event "
	    << evt.id().event()
	    << "\n\nfIndex = " << fIndex
	    << "\nString = " << fAString
	    << "\n\n\n";
  
  fIndex++;
  
  return true;
}



DEFINE_ART_SERVICE_INTERFACE_IMPL(TemplateTriggerService, BaseTriggerService)
