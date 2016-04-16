////////////////////////////////////////////////////////////////////////
// Class:       TemplateTriggerProd
// Module Type: producer
// File:        TemplateTriggerProd_module.cc
//
// Generated at Tuesday March 29 2016 by Michael Baird using the old
// copy and paste method...
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
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dune/DAQTriggerSim/Service/TemplateTriggerService.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"


class TemplateTriggerProd;

class TemplateTriggerProd : public art::EDProducer {

public:

  explicit TemplateTriggerProd(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  TemplateTriggerProd(TemplateTriggerProd const &) = delete;
  TemplateTriggerProd(TemplateTriggerProd &&) = delete;
  TemplateTriggerProd & operator = (TemplateTriggerProd const &) = delete;
  TemplateTriggerProd & operator = (TemplateTriggerProd &&) = delete;

  // The main guts...
  void produce(art::Event& e) override;

  void reconfigure(fhicl::ParameterSet const & p);



private:

  // declare fcl input variables here
  std::string fAString;
  
};



//......................................................
TemplateTriggerProd::TemplateTriggerProd(fhicl::ParameterSet const & p)
{
  produces< std::vector<triggersim::BasicTrigger> >();
  this->reconfigure(p);
}



//......................................................
void TemplateTriggerProd::reconfigure(fhicl::ParameterSet const & p)
{
  fAString = p.get<std::string> ("AString");
}



//......................................................
void TemplateTriggerProd::produce(art::Event& e)
{

  // As an example of how to read in fcl parameters, just print
  // the variable "AString."
  std::cout << "\n\nfAString = " << fAString << "\n\n";

  // make the vector of BasicTrigger objects:
  std::unique_ptr< std::vector<triggersim::BasicTrigger> >
    triggers(new std::vector<triggersim::BasicTrigger>);

  // make the trigger test service:
  art::ServiceHandle<TemplateTriggerService> TT;

  // make the BasicTrigger object:
  triggersim::BasicTrigger trig(triggersim::kProtonDecay,TT->Trigger(e));
  triggers->push_back(trig);

  // put the triggers into the event:
  e.put(std::move(triggers));

}

DEFINE_ART_MODULE(TemplateTriggerProd)
