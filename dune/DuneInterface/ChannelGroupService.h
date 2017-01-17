// ChannelGroupService.h

// David Adams
// October 2016
//
// Interface for a service that groups channels.

#ifndef ChannelGroupService_H
#define ChannelGroupService_H

#include <vector>
#include <iostream>

class ChannelGroupService {

public:

  typedef unsigned int Index;
  typedef unsigned int Channel;
  typedef std::vector<Channel> ChannelVector;
  typedef std::string Name;

  // Return the # groups.
  virtual Index size() const =0;

  // Return a group name.
  virtual Name name(Index igrp) const =0;

  // Return a group.
  virtual const ChannelVector& channels(Index igrp) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(ChannelGroupService, LEGACY)
#endif

#endif

