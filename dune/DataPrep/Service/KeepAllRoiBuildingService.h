// KeepAllRoiBuildingService.h
//
// David Adams
// September 2016
//
// Implementation of service to build ROIs in AdcChannelData that simply
// creates on ROI that includes all samples. Both signal and rois are updated.
//
// Configuration:
//   LogLevel - Usual log level
//
#ifndef KeepAllRoiBuildingService_H
#define KeepAllRoiBuildingService_H

#include "dune/DuneInterface/AdcRoiBuildingService.h"

class AdcSuppressService;

class KeepAllRoiBuildingService : public AdcRoiBuildingService {

public:

  KeepAllRoiBuildingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int build(AdcChannelData& data) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Parameters.
  int m_LogLevel;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(KeepAllRoiBuildingService, AdcRoiBuildingService, LEGACY)

#endif
