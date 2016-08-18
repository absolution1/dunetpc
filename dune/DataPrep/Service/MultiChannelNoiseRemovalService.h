// MultiChannelNoiseRemovalService.h
//
// David Adams
// June 2016
//
// Implementation of service that removesfrom and ADC channel map by calling
// an AdcChannnelNoiseRemovalService for each channel.
//
// Configuration:
//          LogLevel - message logging level: 0=none, 1=initialization, 2+=every event

#ifndef MultiChannelNoiseRemovalService_H
#define MultiChannelNoiseRemovalService_H

#include "dune/DuneInterface/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/AdcTypes.h"

class AdcChannelNoiseRemovalService;

class MultiChannelNoiseRemovalService : public AdcNoiseRemovalService {

public:

  MultiChannelNoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int update(AdcChannelDataMap& datamap) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int  m_LogLevel;

  // The single channel service.
  const AdcChannelNoiseRemovalService* m_pAdcChannelNoiseRemovalService;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(MultiChannelNoiseRemovalService, AdcNoiseRemovalService, LEGACY)

#endif
