// MedianPedestalService.h
//
// David Adams
// June 2016
// Updated Jan 2017 to add UseMean.
//
// Implementation of service that evluates a pedestal as the median of
// an input ADC signal vector.
//
// Configuration:
//             LogLevel - Message logging level: 0=none, 1=initialization, 2+=every event
//              UseMean - Evaluate pedestal with mean instead of median.
//   SkipFlaggedSamples - Exclude ticks for which flag != AdcGood
//          SkipSignals - Exclude ticks identified as signal

#ifndef MedianPedestalService_H
#define MedianPedestalService_H

#include "dune/DuneInterface/PedestalEvaluationService.h"
#include "dune/DuneInterface/AdcTypes.h"

namespace geo {
  class Geometry;
}

namespace lbne {
  class ChannelMapService;
}

class MedianPedestalService : public PedestalEvaluationService {

public:

  MedianPedestalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int evaluate(const AdcChannelData& data,
               AdcSignal* pped =nullptr, AdcSignal* prms =nullptr,
               AdcSignal* ppederr =nullptr, AdcSignal* prmserr =nullptr) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int m_LogLevel;
  bool m_UseMean;
  bool m_SkipFlaggedSamples;
  bool m_SkipSignals;

  // Services.
  //const geo::Geometry* m_pGeometry;	// unused
  //const lbne::ChannelMapService* m_pChannelMap; // unused

  // List of channels for each orientation and group.
  std::vector<std::vector<AdcChannelVector>> m_GroupChannels;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(MedianPedestalService, PedestalEvaluationService, LEGACY)

#endif
