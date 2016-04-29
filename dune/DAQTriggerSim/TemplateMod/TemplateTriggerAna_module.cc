////////////////////////////////////////////////////////////////////////
// Class:       TemplateTriggerAna
// Module Type: analyzer
// File:        TemplateTriggerAna_module.cc
//
// Generated at Fri Feb 19 08:55:48 2016 by Michael Baird using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

// C++ includes

// ROOT includes
#include "TH1F.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"



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

  void beginJob();



private:

  // label for module that made trigger objects
  std::string fTriggerLabel;

  // a simple histo to be filled
  TH1F *fTrigTypes;
  
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
void TemplateTriggerAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  fTrigTypes = tfs->make<TH1F>("fTrigTypes","Trigger Types;trig type;count",
			       101,-0.5,100.5);

}



//......................................................
void TemplateTriggerAna::analyze(art::Event const & e)
{

  //
  // Get trigger data products out of the event
  //
  art::Handle< std::vector< triggersim::BasicTrigger > > triggers;
  e.getByLabel(fTriggerLabel, triggers);

  

  //
  // Loop over trigger objects, print some info, and fill a histo...
  //
  for(unsigned int i = 0; i < triggers->size(); ++i) {

    std::cout << "\n\nInfo for trigger " << i << ":"
	      << "\nTrigger Type     = " << (*triggers)[i].TrigType()
	      << "\nTrigger Sub-Type = " << (*triggers)[i].TrigSubType()
	      << "\nTrigger Decision = " << (*triggers)[i].TrigDecision()
	      << "\n\n\n";

    fTrigTypes->Fill((*triggers)[i].TrigType());

  }

}

DEFINE_ART_MODULE(TemplateTriggerAna)
