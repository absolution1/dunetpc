// GeoApaChannelGroupService.h

// David Adams
// October 2016
//
// APA channel groups extracted from larsoft geometry.
// Group igrp holds the channels from APA igrp.
//
// The array ApaRops specifies which ROPs are kept for each group.
// If it is empty, all ROPS are retained.
// Otherwise the bit pattern at ApaRops[igrp%nent] (nent = ApaRops.size())
// specifies which ROPS to to keep. E.g. 15 = 0xa = b_1111 to keep all four or
// 11 = 0xb = 1011 to keep ROPs 0, 1 and 3.
// Note the DUNE plane ordering from highest to lowest is x2-x1-v-u so that
// 7 keeps planes from the first TPC and 11 keeps the planes for the 2nd TPC.
//
// Parameters:
//   LogLevel - 0=silent, 1=init only, 2=every call
//   ApaRops - Bit pattern indicating which ROPs to keep (see above).

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

  // Parameters.
  int m_LogLevel;
  std::vector<int> m_ApaRops;

  Index m_size;
  NameVector m_names;
  ChannelVectorVector m_chanvecs;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(GeoApaChannelGroupService, ChannelGroupService, LEGACY)

#endif
