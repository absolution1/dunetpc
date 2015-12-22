// ChannelNoiseService.h

// David Adams
// December 2015
//
// Interface for a service that adds noise to a TPC signal vector.

#ifndef ChannelNoiseService_H
#define ChannelNoiseService_H

#include <vector>
#include <iostream>
#include "DZSInterface/AdcTypes.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

class ChannelNoiseService {

public:

  typedef unsigned int Channel;

  virtual ~ChannelNoiseService() =default;

  // Add noise to a signal vector sigs appropriate for channel chan.
  // Noise is added for all entries in the input vector.
  virtual int addNoise(Channel chan, AdcSignalVector& sigs) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

DECLARE_ART_SERVICE_INTERFACE(ChannelNoiseService, LEGACY)

#endif

