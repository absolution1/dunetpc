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
#include "TTree.h"

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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"

namespace triggersim {

class TemplateTriggerAna : public art::EDAnalyzer {

public:

  explicit TemplateTriggerAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  TemplateTriggerAna(TemplateTriggerAna const &) = delete;
  TemplateTriggerAna(TemplateTriggerAna &&) = delete;
  TemplateTriggerAna & operator = (TemplateTriggerAna const &) = delete;
  TemplateTriggerAna & operator = (TemplateTriggerAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & evt) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;

  void ResetVars();

private:

  // label for module that made trigger objects
  std::string fTriggerLabel;

  // a simple histo to be filled
  TH1F *fTrigTypes;

  // The name of our TTree, and the variables which we want to put in it.
  TTree *TrigAnaTree;
  int  Run;           ///< The run which the trigger was taken from
  int  SubRun;        ///< The SubRun which the trigger was taken from
  int  Event;         ///< The event which the trigger was taken from
  bool TrigDecision;  ///< The decision of the trigger
  int  TrigType;      ///< The type of trigger decision which was made
  int  TrigSubType;   ///< The subtype of trigger decision which was made
};

//......................................................
void TemplateTriggerAna::ResetVars() {
  Run          = -1;
  SubRun       = -1;
  Event        = -1;
  TrigDecision = false;
  TrigType     = -1;
  TrigSubType  = -1;
}

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
  // --- Make a handle to the TFileService - the thing that handles all the histogram and TTree outputs...
  art::ServiceHandle<art::TFileService> tfs;
  // --- Now we can declare the histrogram which we want to make.
  fTrigTypes = tfs->make<TH1F>("fTrigTypes","Trigger Types;trig type;count",
			       101,-0.5,100.5);

  // --- Now we want to declare our TTree, and set the names of the variables in it...
  TrigAnaTree = tfs->make<TTree>("TrigAnaTree", "Analysis tree from triggers");
  TrigAnaTree -> Branch( "Run"         , &Run          );
  TrigAnaTree -> Branch( "SubRun"      , &SubRun       );
  TrigAnaTree -> Branch( "Event"       , &Event        );
  TrigAnaTree -> Branch( "TrigDecision", &TrigDecision );
  TrigAnaTree -> Branch( "TrigType"    , &TrigType     );
  TrigAnaTree -> Branch( "TrigSubType" , &TrigSubType  );  
}



//......................................................
void TemplateTriggerAna::analyze(art::Event const & evt) {

  // --- Reset all of our variables...
  ResetVars();
  
  // --- Get trigger data products out of the event...
  auto triggers = evt.getValidHandle<std::vector<triggersim::BasicTrigger> >(fTriggerLabel);

  // --- Loop over trigger objects...
  for(unsigned int i = 0; i < triggers->size(); ++i) {

    // --- Set our variables for the TTree...
    Run    = evt.run();
    SubRun = evt.subRun();
    Event  = evt.event();

    TrigDecision = (*triggers)[i].TrigDecision();
    TrigType     = (*triggers)[i].TrigType();
    TrigSubType  = (*triggers)[i].TrigSubType();
   // --- Print some info, and confirm to ourselves that our variables are correctly set...
    std::cout << "\n----------Info for trigger " << i << ":"
	      << "\nTrigger Decision   = " << (*triggers)[i].TrigDecision() << " " << TrigDecision
	      << "\nTrigger Type       = " << (*triggers)[i].TrigType()     << " " << TrigType
	      << "\nTrigger Sub-Type   = " << (*triggers)[i].TrigSubType()  << " " << TrigSubType
	      << "\nTrigger HardwareID = " << (*triggers)[i].TrigHardwareID()
	      << std::endl;


    std::cout << "\nAnother way of printing out a trigger" << std::endl;
    std::cout << (*triggers)[i] << std::endl;

    // --- Fill a histogram...
    fTrigTypes->Fill((*triggers)[i].TrigType());

    // --- Fill a TTree...
    TrigAnaTree -> Fill();
  }

}

DEFINE_ART_MODULE(TemplateTriggerAna)

}
