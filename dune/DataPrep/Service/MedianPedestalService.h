// MedianPedestalService.h
//
// David Adams
// June 2016
//
// Implementation of service that evluates a pedestal as the median of
// an input ADC signal vector.
//
// Configuration:
//          LogLevel - message logging level: 0=none, 1=initialization, 2+=every event

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

  int evaluate(const AdcSignalVector& sigs,
               AdcSignal* pped =nullptr, AdcSignal* prms =nullptr,
               AdcSignal* ppederr =nullptr, AdcSignal* prmserr =nullptr) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int  m_LogLevel;
  int  m_GroupingFlag;
  bool m_SkipStuckCodes;
  bool m_CorrectStuckCodes;
  int  m_ShowGroups;

  // Services.
  const geo::Geometry* m_pGeometry;
  const lbne::ChannelMapService* m_pChannelMap;

  // List of channels for each orientation and group.
  std::vector<std::vector<AdcChannelVector>> m_GroupChannels;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(MedianPedestalService, PedestalEvaluationService, LEGACY)

#endif
