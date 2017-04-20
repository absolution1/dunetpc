// AdcChannelDataCopyService.h

#ifndef AdcChannelDataCopyService_H
#define AdcChannelDataCopyService_H

// David Adams
// September 2016
//
// Interface for a service that copies AdcChannelData.

#include <iostream>
#include "dune/DuneInterface/AdcChannelData.h"

class AdcChannelDataCopyService {

public:

  // Copy the data from one AdcChannelData object to another.
  // Returns nonzero for error.
  virtual int copy(const AdcChannelData& oldacd, AdcChannelData& newacd) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcChannelDataCopyService, LEGACY)
#endif

#endif
