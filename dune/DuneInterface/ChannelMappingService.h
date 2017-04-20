// ChannelMappingService.h

// David Adams
// December 2015
//
// Interface for a service that maps online to offline channels.

#ifndef ChannelMappingService_H
#define ChannelMappingService_H

#include <vector>
#include <iostream>
#include "dune/DuneInterface/AdcTypes.h"

class ChannelMappingService {

public:

  typedef unsigned int Channel;

  // Map online to offline.
  virtual Channel offline(Channel onlineChannel) const =0;

  // Map offline to online.
  virtual Channel online(Channel offlineChannel) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(ChannelMappingService, LEGACY)
#endif

#endif

