// ThresholdNoiseRemovalService.h
//
// David Adams
// August 2016
//
// Implementation of service that removes samples with magnitude below a threshold.
// Intended for demonstration and test.
//
// Configuration:
//          LogLevel - message logging level: 0=none, 1=initialization, 2+=every event
//         Threshold - samples with |value| < Threshold are set to zero

#ifndef ThresholdNoiseRemovalService_H
#define ThresholdNoiseRemovalService_H

#include "dune/DuneInterface/AdcChannelNoiseRemovalService.h"
#include "dune/DuneInterface/AdcTypes.h"

namespace geo {
  class Geometry;
}

namespace lbne {
  class ChannelMapService;
}

class ThresholdNoiseRemovalService : public AdcChannelNoiseRemovalService {

public:

  ThresholdNoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int update(AdcChannelData& data) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int  m_LogLevel;
  AdcSignal m_Threshold;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ThresholdNoiseRemovalService, AdcChannelNoiseRemovalService, LEGACY)

#endif
