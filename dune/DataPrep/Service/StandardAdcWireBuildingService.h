// StandardAdcWireBuildingService.h
//
// David Adams
// June 2016
//
// Implementation of service to build the recob::Wire in AdcChannelData.
//
// Configuration:
//   LogLevel    - usual log level
//
#ifndef StandardAdcWireBuildingService_H
#define StandardAdcWireBuildingService_H

#include "dune/DuneInterface/AdcWireBuildingService.h"

class AdcSuppressService;

class StandardAdcWireBuildingService : public AdcWireBuildingService {

public:

  StandardAdcWireBuildingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int build(AdcChannelData& data, WireVector* pwires) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Parameters.
  int m_LogLevel;
  bool m_SaveChanPedRMS;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(StandardAdcWireBuildingService, AdcWireBuildingService, LEGACY)

#endif
