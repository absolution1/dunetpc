// DuneDPhase3x1x1NoiseRemovalService.h
//
// Robert Sulej, Aug 2017
//
// Remove coherent noise from 3x1x1 data.
//

#ifndef DuneDPhase3x1x1NoiseRemovalService_H
#define DuneDPhase3x1x1NoiseRemovalService_H

#include "dune/DuneInterface/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/AdcTypes.h"

namespace geo {
  class Geometry;
}

class DuneDPhase3x1x1NoiseRemovalService : public AdcNoiseRemovalService {

public:

  DuneDPhase3x1x1NoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int update(AdcChannelDataMap& datamap) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix = "") const;

private:

  // Configuration parameters.

  // Services.
  const geo::Geometry* fGeometry;

  // List of channels for each orientation and group.
  //std::vector<std::vector<AdcChannelVector>> m_GroupChannels;
};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DuneDPhase3x1x1NoiseRemovalService, AdcNoiseRemovalService, LEGACY)

#endif
