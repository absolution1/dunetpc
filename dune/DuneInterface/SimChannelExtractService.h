// SimChannelExtractService.h

// David Adams
// December 2015
//
// Interface for a service that extracts charge for one channel from
// a SimChannel object and assigns it to ticks.

#ifndef SimChannelExtractService_H
#define SimChannelExtractService_H

#include <vector>
#include <iostream>
#include "dune/DuneInterface/AdcTypes.h"

namespace sim {
class SimChannel;
}

class SimChannelExtractService {

public:

  virtual ~SimChannelExtractService() =default;

  virtual int extract(const sim::SimChannel* psc, AdcSignalVector& sig) const =0;

  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(SimChannelExtractService, LEGACY)
#endif

#endif

