// Dune35tChannelGroupService.h

// David Adams
// January 2017
//
// Channel groups for DUNE 35t. Grouping may be by APA, Regulator or ASIC.
// Channels may be online or offline. Group numbers are online.
//
// Parameters:
//   UseOffline - true for offline channel numbers, false for online
//   Grouping - Grouping option: APA, ASIC or Regulator
//   SplitByPlane - split groups by plane (u, v, z1, z2)
//   LogLevel - 0=quiet, 1=init only, 2=every call

#ifndef Dune35tChannelGroupService_H
#define Dune35tChannelGroupService_H

#include "dune/DuneInterface/ChannelGroupService.h"
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

class Dune35tChannelGroupService : public ChannelGroupService {

public:

  // Ctor.
  Dune35tChannelGroupService(fhicl::ParameterSet const& pset);

  // Ctor.
  Dune35tChannelGroupService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~Dune35tChannelGroupService() =default;

  // Return the # groups.
  Index size() const;

  // Return a group name.
  Name name(Index igrp) const;

  // Return a group.
  const ChannelVector& channels(Index igrp) const;

  // Print parameters.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  bool m_UseOffline;
  std::string m_Grouping;      // APA, ASIC, Regulator.
  bool m_SplitByPlane;
  int m_LogLevel;

  typedef std::vector<Name> NameVector;
  typedef std::vector<ChannelVector> ChannelVectorVector;
 
  art::ServiceHandle<lbne::ChannelMapService> m_pcms;

  NameVector m_names;
  ChannelVectorVector m_chanvecs;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(Dune35tChannelGroupService, ChannelGroupService, LEGACY)

#endif
