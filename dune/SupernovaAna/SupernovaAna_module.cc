////////////////////////////////////////////////////////////////////////
// Class:       SupernovaAna
// Module Type: analyzer
// File:        SupernovaAna_module.cc
//
// Generated at Mon Jul 11 21:36:48 2016 by Michael Baird using the old
// copy and paste...
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



class SupernovaAna;

class SupernovaAna : public art::EDAnalyzer {

public:

  explicit SupernovaAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  SupernovaAna(SupernovaAna const &) = delete;
  SupernovaAna(SupernovaAna &&) = delete;
  SupernovaAna & operator = (SupernovaAna const &) = delete;
  SupernovaAna & operator = (SupernovaAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob();



private:

  // label for module that made trigger objects
  std::string fTruthLabel;

  // a simple histo to be filled
  TH1F *fNMCTruths;
  
};



//......................................................
SupernovaAna::SupernovaAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}



//......................................................
void SupernovaAna::reconfigure(fhicl::ParameterSet const & p)
{
  fTruthLabel = p.get<std::string> ("TruthLabel");
}



//......................................................
void SupernovaAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  fNMCTruths = tfs->make<TH1F>("fNMCTruths","Number of MCTruths per event;N;count",
			       101,-0.5,100.5);

}



//......................................................
void SupernovaAna::analyze(art::Event const & e)
{

  //
  // Get the MCTruths out of the event
  //
  // art::Handle< std::vector< triggersim::BasicTrigger > > triggers;
  // e.getByLabel(fTriggerLabel, triggers);
  std::cout << "\n\n\n" << fTruthLabel << "\n\n\n";

  

  // do something...

}

DEFINE_ART_MODULE(SupernovaAna)
