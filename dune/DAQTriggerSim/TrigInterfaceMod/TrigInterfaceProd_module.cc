////////////////////////////////////////////////////////////////////////
// Class:       TrigInterfaceProd
// Module Type: producer
// File:        TrigInterfaceProd_module.cc
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// Brief - Example of how to use BaseTriggerService as an "Interface" 
//       - This means the user can configure which trigger to use via
//         fcl configuration, and hence change the trigger at run time
//       - Changing triginterfaceprod_job.fcl:
//           services.user.BaseTriggerService: @local::activitytrigger
//           services.user.BaseTriggerService: @local::templatetrigger

//
////////////////////////////////////////////////////////////////////////

// C++ includes

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dune/DuneInterface/BaseTriggerService.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"



class TrigInterfaceProd;

class TrigInterfaceProd : public art::EDProducer {

public:

  explicit TrigInterfaceProd(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  TrigInterfaceProd(TrigInterfaceProd const &) = delete;
  TrigInterfaceProd(TrigInterfaceProd &&) = delete;
  TrigInterfaceProd & operator = (TrigInterfaceProd const &) = delete;
  TrigInterfaceProd & operator = (TrigInterfaceProd &&) = delete;

  // The main guts...
  void produce(art::Event& e) override;

  void reconfigure(fhicl::ParameterSet const & p);



private:

  // declare fcl input variables here
  std::string fAString;
  std::string fRawDigitLabel;
  
};



//......................................................
TrigInterfaceProd::TrigInterfaceProd(fhicl::ParameterSet const & p)
{
  produces< std::vector<triggersim::BasicTrigger> >();
  this->reconfigure(p);
}



//......................................................
void TrigInterfaceProd::reconfigure(fhicl::ParameterSet const & p)
{
  fAString       = p.get<std::string> ("AString");
  fRawDigitLabel = p.get<std::string> ("RawDigitLabel");
}



//......................................................
void TrigInterfaceProd::produce(art::Event& e)
{

  // As an example of how to read in fcl parameters, just print
  // the variable "AString."
  
  // make the vector of BasicTrigger objects:
  std::unique_ptr< std::vector<triggersim::BasicTrigger> >
    triggers(new std::vector<triggersim::BasicTrigger>);

  // make the trigger test service:
  art::ServiceHandle<BaseTriggerService> trigService;
  

  std::cout << "INFO  : " << "trigService->GetName(): " << trigService->GetName() << std::endl;
  std::cout << "INFO  : " << "fAString: " << fAString << std::endl;

  // make the BasicTrigger object:
  triggersim::BasicTrigger trig(trigService->Trigger(e),triggersim::kNucleon,triggersim::kProtonDecay);
  triggers->push_back(trig);

  // lift out the TPC raw digits:
  art::Handle<std::vector<raw::RawDigit>> digitsHandle;
  e.getByLabel(fRawDigitLabel, digitsHandle);

  // shuffle the digits into an easier container to use:
  std::vector<raw::RawDigit> rawdigits;
  for(unsigned int i = 0; i < digitsHandle->size(); ++i) {
    rawdigits.push_back((*digitsHandle)[i]);
  }



  // make a TPC only trigger object:
  triggersim::BasicTrigger tpctrig(trigService->TPCTrigger(rawdigits),triggersim::kNu,triggersim::kNuBeam);
  triggers->push_back(tpctrig);

  // put the triggers into the event:
  e.put(std::move(triggers));

}

DEFINE_ART_MODULE(TrigInterfaceProd)
