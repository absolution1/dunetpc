////////////////////////////////////////////////////////////////////////
// LifetimeCalibService.h
//
// Pure virtual service interface for lifetime calibration functions
//
//  wwu@fnal.gov
//  date: May 19, 2020
//
////////////////////////////////////////////////////////////////////////
#ifndef LIFETIMECALIBSERVICE_H
#define LIFETIMECALIBSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "dune/Calib/LifetimeCalib.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace calib{
    class LifetimeCalibService {
      public:
      typedef calib::LifetimeCalib provider_type;

      public:
      virtual ~LifetimeCalibService() = default;

      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual calib::LifetimeCalib* provider() const = 0;

      }; // class LifetimeCalibService
    } //namespace detinfo
DECLARE_ART_SERVICE_INTERFACE(calib::LifetimeCalibService, LEGACY)
#endif // LIFETIMECALIBSERVICE_H
