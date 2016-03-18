////////////////////////////////////////////////////////////////////////
// Class:       TemplateTriggerAna
// Module Type: analyzer
// File:        TemplateTriggerAna_module.cc
//
// Generated at Fri Feb 19 08:55:48 2016 by Michael Baird using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

// C++ includes

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
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



class TemplateTriggerAna;

class TemplateTriggerAna : public art::EDAnalyzer {

public:

  explicit TemplateTriggerAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  TemplateTriggerAna(TemplateTriggerAna const &) = delete;
  TemplateTriggerAna(TemplateTriggerAna &&) = delete;
  TemplateTriggerAna & operator = (TemplateTriggerAna const &) = delete;
  TemplateTriggerAna & operator = (TemplateTriggerAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p);



private:

  // label for module that made trigger objects
  std::string fTriggerLabel;
  
};



//......................................................
TemplateTriggerAna::TemplateTriggerAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}



//......................................................
void TemplateTriggerAna::reconfigure(fhicl::ParameterSet const & p)
{
  fTriggerLabel = p.get<std::string> ("TriggerLabel");
}



//......................................................
void TemplateTriggerAna::analyze(art::Event const & e)
{

  // make the trigger test service
  art::ServiceHandle<TemplateTriggerService> TT;
  TT->Trigger(e);

  std::cout << "\n\n\nThe TemplateTriggerAna module says: \"Hello world!!!\"\n"
	    << fTriggerLabel << "\n\n\n";

}

DEFINE_ART_MODULE(TemplateTriggerAna)
