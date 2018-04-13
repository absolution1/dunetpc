// FixedChannelGroupService.h

// David Adams
// November 2016
//
// APA channel groups with fixed mapping defined in FCL.

#ifndef FixedChannelGroupService_H
#define FixedChannelGroupService_H

#include "dune/DuneInterface/ChannelGroupService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

class FixedChannelGroupService : public ChannelGroupService {

public:

  // Ctor.
  FixedChannelGroupService(fhicl::ParameterSet const& pset);

  // Ctor.
  FixedChannelGroupService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~FixedChannelGroupService() =default;

  // Return the # groups.
  Index size() const;

  // Return a group name.
  Name name(Index igrp) const;

  // Return a group.
  const ChannelVector& channels(Index igrp) const;

  // Print parameters.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  typedef std::vector<Name> NameVector;
  typedef std::vector<ChannelVector> ChannelVectorVector;
 
  //Index m_size; // unused
  NameVector m_names;
  ChannelVectorVector m_chanvecs;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(FixedChannelGroupService, ChannelGroupService, LEGACY)

#endif
