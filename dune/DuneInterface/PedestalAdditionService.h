// PedestalAdditionService.h

// David Adams
// December 2015
//
// Interface for a service that add pedestals to an array
// ADC signals.

#ifndef PedestalAdditionService_H
#define PedestalAdditionService_H

#include <string>
#include <iostream>
#include <vector>
#include "dune/DuneInterface/AdcTypes.h"

namespace sim {
class SimChannel;
}

class PedestalAdditionService {

public:

  typedef unsigned int Channel;

  // Dtor.
  virtual ~PedestalAdditionService() =default;

  // Add pedestal to each tick in sigs.
  // The pedestal and its RMS are returned in ped and pedrms.
  virtual int addPedestal(Channel chan, AdcSignalVector& sigs, float& ped, float& pedrms) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(PedestalAdditionService, LEGACY)
#endif

#endif

