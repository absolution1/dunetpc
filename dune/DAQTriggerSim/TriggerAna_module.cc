////////////////////////////////////////////////////////////////////////
// Class:       TriggerAna
// Module Type: analyzer
// File:        TriggerAna_module.cc
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
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes



class TriggerAna;

class TriggerAna : public art::EDAnalyzer {

public:

  explicit TriggerAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  TriggerAna(TriggerAna const &) = delete;
  TriggerAna(TriggerAna &&) = delete;
  TriggerAna & operator = (TriggerAna const &) = delete;
  TriggerAna & operator = (TriggerAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p);



private:

  // label for module that made trigger objects
  std::string fTriggerLabel;

};



//......................................................
TriggerAna::TriggerAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}



//......................................................
void TriggerAna::reconfigure(fhicl::ParameterSet const & p)
{
  fTriggerLabel = p.get<std::string> ("TriggerLabel");
}



//......................................................
void TriggerAna::analyze(art::Event const & e)
{

  std::cout << "\n\n\nHello world!!!\n\n\n"
	    << fTriggerLabel << "\n\n\n";

}

DEFINE_ART_MODULE(TriggerAna)
