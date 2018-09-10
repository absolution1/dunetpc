////////////////////////////////////////////////////////////////////////
// XYZCalibService.h
//
// Pure virtual service interface for XYZ calibration functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef XYZCALIBSERVICE_H
#define XYZCALIBSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "dune/Calib/XYZCalib.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace calib{
    class XYZCalibService {
      public:
      typedef calib::XYZCalib provider_type;

      public:
      virtual ~XYZCalibService() = default;

      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual calib::XYZCalib* provider() const = 0;

      }; // class XYZCalibService
    } //namespace detinfo
DECLARE_ART_SERVICE_INTERFACE(calib::XYZCalibService, LEGACY)
#endif // XYZCALIBSERVICE_H
