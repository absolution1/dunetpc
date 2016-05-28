////////////////////////////////////////////////////////////////////////
// Class:       ActTrigProd
// Module Type: producer
// File:        ActTrigProd_module.cc
//
// Generated on Friday May 20 2016 by Michael Baird using the old
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
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dune/DAQTriggerSim/Service/ActivityTriggerService.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"

// ROOT includes
#include "TH1F.h"



class ActTrigProd;

class ActTrigProd : public art::EDProducer {

public:

  explicit ActTrigProd(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  ActTrigProd(ActTrigProd const &) = delete;
  ActTrigProd(ActTrigProd &&) = delete;
  ActTrigProd & operator = (ActTrigProd const &) = delete;
  ActTrigProd & operator = (ActTrigProd &&) = delete;

  void beginJob();

  // The main guts...
  void produce(art::Event& e) override;

  void reconfigure(fhicl::ParameterSet const & p);



private:

  // declare fcl input variables here
  std::string fRawDigitLabel;
  
  TH1F *fSpec0;
  TH1F *fSpec1;
  TH1F *fSpec2;
  TH1F *fSpec3;

};



//......................................................
ActTrigProd::ActTrigProd(fhicl::ParameterSet const & p)
{
  produces< std::vector<triggersim::BasicTrigger> >();
  this->reconfigure(p);
}



//......................................................
void ActTrigProd::reconfigure(fhicl::ParameterSet const & p)
{
  fRawDigitLabel = p.get<std::string> ("RawDigitLabel");
}



//......................................................
void ActTrigProd::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fSpec0 = tfs->make<TH1F>("fSpec0",
			   "ADC readout - wire 1600",
			   1025,-0.5,1024.5);

  fSpec1 = tfs->make<TH1F>("fSpec1",
			   "ADC readout - wire 1601",
			   1025,-0.5,1024.5);

  fSpec2 = tfs->make<TH1F>("fSpec2",
			   "ADC readout - wire 1602",
			   1025,-0.5,1024.5);

  fSpec3 = tfs->make<TH1F>("fSpec3",
			   "ADC readout - wire 1603",
			   1025,-0.5,1024.5);
}



//......................................................
void ActTrigProd::produce(art::Event& e)
{

  // make the vector of BasicTrigger objects:
  std::unique_ptr< std::vector<triggersim::BasicTrigger> >
    triggers(new std::vector<triggersim::BasicTrigger>);

  // make the trigger service:
  art::ServiceHandle<ActivityTriggerService> actTrig;





  // lift out the TPC raw digits:
  art::Handle<std::vector<raw::RawDigit>> digitsHandle;
  e.getByLabel(fRawDigitLabel, digitsHandle);

  // shuffle the digits into an easier container to use:
  std::vector<raw::RawDigit> rawdigits;
  for(unsigned int i = 0; i < digitsHandle->size(); ++i) {
    rawdigits.push_back((*digitsHandle)[i]);
  }



  // fill a few plots with ADC values...
  std::vector<short> ADCvec0 = rawdigits[1600].ADCs();
  std::vector<short> ADCvec1 = rawdigits[1601].ADCs();
  std::vector<short> ADCvec2 = rawdigits[1602].ADCs();
  std::vector<short> ADCvec3 = rawdigits[1603].ADCs();
  for(unsigned int i = 0; i < ADCvec0.size(); ++i) fSpec0->Fill(i,ADCvec0[i]);
  for(unsigned int i = 0; i < ADCvec1.size(); ++i) fSpec1->Fill(i,ADCvec1[i]);
  for(unsigned int i = 0; i < ADCvec2.size(); ++i) fSpec2->Fill(i,ADCvec2[i]);
  for(unsigned int i = 0; i < ADCvec3.size(); ++i) fSpec3->Fill(i,ADCvec3[i]);



  // make the trigger data product:
  // (just use kNullTrigger for the trigger subtype for now...)
  triggersim::BasicTrigger trigger(actTrig->TPCTrigger(rawdigits),triggersim::kActivity,triggersim::kNullTrigger);
  triggers->push_back(trigger);

  // put the triggers into the event:
  e.put(std::move(triggers));

}

DEFINE_ART_MODULE(ActTrigProd)
