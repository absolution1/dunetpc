// GeoApaChannelGroupService.h

// David Adams
// October 2016
//
// APA channel groups extracted from larsoft geometry.

#ifndef GeoApaChannelGroupService_H
#define GeoApaChannelGroupService_H

#include "dune/DuneInterface/ChannelGroupService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

class GeoApaChannelGroupService : public ChannelGroupService {

public:

  // Ctor.
  GeoApaChannelGroupService(fhicl::ParameterSet const& pset);

  // Ctor.
  GeoApaChannelGroupService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~GeoApaChannelGroupService() =default;

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
 
  art::ServiceHandle<geo::Geometry> m_pgeo;

  Index m_size;
  NameVector m_names;
  ChannelVectorVector m_chanvecs;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(GeoApaChannelGroupService, ChannelGroupService, LEGACY)

#endif
