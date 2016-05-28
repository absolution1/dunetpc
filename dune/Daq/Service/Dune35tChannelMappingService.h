// Dune35tChannelMappingService

// David Adams
// May 2016
//
// Implementation of online-offline channel mapping that use the service in lbne_raw_data.
//
// Parameters:
//   LogLevel [1] - Message logging level: 0=none, 1=init, 2+=every call to interface
//
// The parameter LogLevel is optional and has the indicated default values.

#ifndef Dune35tChannelMappingService_H
#define Dune35tChannelMappingService_H

#include "dune/DuneInterface/ChannelMappingService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class Dune35tChannelMappingService : public ChannelMappingService {

  typedef std::vector<Channel> ChannelMap;

public:

  // Ctor.
  Dune35tChannelMappingService(fhicl::ParameterSet const& pset);

  // Ctor.
  Dune35tChannelMappingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Map online to offline.
  Channel offline(Channel onlineChannel) const;

  // Map offline to online.
  Channel online(Channel offlineChannel) const;

  // Print the parameters.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  int m_LogLevel;

  art::ServiceHandle<lbne::ChannelMapService> hlbnesvc;
 
};

DECLARE_ART_SERVICE_INTERFACE_IMPL(Dune35tChannelMappingService, ChannelMappingService, LEGACY)

#endif
