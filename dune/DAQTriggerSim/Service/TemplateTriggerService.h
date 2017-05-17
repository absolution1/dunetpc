// TemplateTriggerService.h
//
// Michael Baird
// March 2016
//
// Template trigger service class for the DAQ trigger service framework.
//

#ifndef TEMPLATETRIGGERSERVICE_H
#define TEMPLATETRIGGERSERVICE_H

#include "dune/DuneInterface/BaseTriggerService.h"

#include <memory>
#include <string>
#include <iostream>

namespace fhicl {
  class ParameterSet;
}
namespace art {
  class ActivityRegistry;
}

class TemplateTriggerService : public BaseTriggerService {
  
 public:
  
  // Simple constructor...
  TemplateTriggerService();
  
  // Constructor from fcl parameters...
  TemplateTriggerService(fhicl::ParameterSet const & pset, art::ActivityRegistry&);
  
  // Make the trigger decision.
  //
  //   params:
  //   evt - The ART event, to be sliced and diced as the trigger algorithm
  //         sees fit...
  bool Trigger(art::Event const & evt);

  // Make a trigger decision based on TPC info only.
  //
  //   params:
  //   raw - A vector of raw digits, to be sliced and diced as the trigger algorithm
  //         sees fit...
  virtual bool TPCTrigger(std::vector<raw::RawDigit> const & raw);

  virtual std::string GetName(){ return "TemplateTriggerService";}

 private:
  
  // Parameters.
  unsigned int fIndex;
  std::string  fAString;
  
};

DECLARE_ART_SERVICE_INTERFACE_IMPL(TemplateTriggerService, BaseTriggerService, LEGACY)

#endif
