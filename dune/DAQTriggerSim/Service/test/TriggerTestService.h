// TriggerTestService.h
//
// Michael Baird
// March 2016
//
// Test service class for the DAQ trigger service framework.
//

#ifndef TRIGGERTESTSERVICE_H
#define TRIGGERTESTSERVICE_H

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

class TriggerTestService : public BaseTriggerService {

public:

  // Simple constructor...
  TriggerTestService();

  // Constructor from fcl parameters...
  TriggerTestService(fhicl::ParameterSet const & pset, art::ActivityRegistry&);

  // Make a trigger decision based on TPC info only.
  //
  //   params:
  //   evt - The ART event, to be sliced and diced as the trigger algorithm
  //         sees fit...
  bool TPCtrigger(art::Event const & evt);

  // Make a trigger decision based on PD info only.
  //
  //   params:
  //   evt - The ART event, to be sliced and diced as the trigger algorithm
  //         sees fit...
  bool PDtrigger(art::Event const & evt);

private:

  // Parameters.
  unsigned int fIndex;
  std::string  fAString;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(TriggerTestService, BaseTriggerService, LEGACY)

#endif
