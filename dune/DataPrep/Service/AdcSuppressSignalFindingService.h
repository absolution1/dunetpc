// AdcSuppressSignalFindingService.h
//
// David Adams
// June 2016
//
// Implementation of service to find signals in AdcChannelData using
// an AdcSuppressService.
//
#ifndef AdcSuppressSignalFindingService_H
#define AdcSuppressSignalFindingService_H

#include "dune/DuneInterface/AdcSignalFindingService.h"

class AdcSuppressService;

class AdcSuppressSignalFindingService : public AdcSignalFindingService {

public:

  AdcSuppressSignalFindingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int find(AdcChannelData& data) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Services.
  const AdcSuppressService* m_psup;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(AdcSuppressSignalFindingService, AdcSignalFindingService, LEGACY)

#endif
