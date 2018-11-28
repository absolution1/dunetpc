////////////////////////////////////////////////////////////////////////
// TPCHVService.h
//
// Pure virtual service interface for slow control TPC HV data
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef TPCHVSERVICE_H
#define TPCHVSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "dune/SlowControls/SlowControls.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace slowctrls{
    class TPCHVService {
      public:
      typedef slowctrls::SlowControls provider_type;

      public:
      virtual ~TPCHVService() = default;

      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual provider_type* provider() = 0;

      }; // class TPCHVService
    } //namespace slowctrl
DECLARE_ART_SERVICE_INTERFACE(slowctrls::TPCHVService, LEGACY)
#endif // TPCHVSERVICE_H
