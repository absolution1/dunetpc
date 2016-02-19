// BaseTriggerService.h
//
// Michael Baird
// Feb. 2016
//
// Base trigger service class for other trigger services to inherit from.

#ifndef BaseTriggerService_H
#define BaseTriggerService_H

#include <memory>
#include "larcore/SimpleTypesAndConstants/RawTypes.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Event.h"

namespace fhicl {
  class ParameterSet;
}

namespace art {
  class ActivityRegistry;
}

class BaseTriggerService {

public:

  // Destructor:
  virtual ~BaseTriggerService() = default;



  // Make a trigger decision based on TPC info only.
  //
  //   params:
  //   evt - The ART event, to be sliced and diced as the trigger algorithm
  //         sees fit...
  virtual bool TPCtrigger(art::Event const & evt)
  { return false; }



  // Make a trigger decision based on PD info only.
  //
  //   params:
  //   evt - The ART event, to be sliced and diced as the trigger algorithm
  //         sees fit...
  virtual bool PDtrigger(art::Event const & evt)
  { return false; }



  // Make a trigger decision from both TPC and PD info.
  //
  //   params:
  //   evt - The ART event, to be sliced and diced as the trigger algorithm
  //         sees fit...
  virtual bool TPCPDtrigger(art::Event const & evt)
  { return false; }

};

DECLARE_ART_SERVICE_INTERFACE(BaseTriggerService, LEGACY)

#endif
