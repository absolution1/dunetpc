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
#include "../Service/TemplateTriggerService.h"
#include "../TriggerDataProducts/TriggerTypes.h"
#include "../TriggerDataProducts/BasicTrigger.h"


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

  // label for module that made trigger objects
  std::string fTriggerLabel;
  
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
  fTriggerLabel = p.get<std::string> ("TriggerLabel");
}



//......................................................
void TemplateTriggerProd::produce(art::Event& e)
{

  // make the vector of BasicTrigger objects:
  std::unique_ptr< std::vector<triggersim::BasicTrigger> >
    triggers(new std::vector<triggersim::BasicTrigger>);

  // make the trigger test service:
  art::ServiceHandle<TemplateTriggerService> TT;

  // make the BasicTrigger object:
  triggersim::BasicTrigger trig(triggersim::kNullTrigger,TT->Trigger(e));
  triggers->push_back(trig);

  std::cout << "\n\n\nThe TemplateTriggerProd module says: \"Hello world!!!\"\n"
	    << fTriggerLabel << "\n\n\n";

  // put the triggers into the event:
  e.put(std::move(triggers));

}

DEFINE_ART_MODULE(TemplateTriggerProd)
