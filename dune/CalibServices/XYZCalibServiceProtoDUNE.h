////////////////////////////////////////////////////////////////////////
// \file XYZCalibServiceProtoDUNE.h
//
// \brief header of service for storing/accessing (x,y,z) calibration corrections for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef XYZCALIBSERVICEPROTODUNE_H
#define XYZCALIBSERVICEPROTODUNE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "dune/Calib/XYZCalibProtoDUNE.h"
#include "dune/CalibServices/XYZCalibService.h"

namespace calib{
  class XYZCalibServiceProtoDUNE : public XYZCalibService {
    public:
      
      XYZCalibServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset)  override;
      void   preBeginRun(const art::Run& run);


      virtual provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<calib::XYZCalibProtoDUNE> fProp;

    }; // class XYZCalibServiceProtoDUNE
} //namespace calib
DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::XYZCalibServiceProtoDUNE, calib::XYZCalibService, LEGACY)
#endif // XYZCALIBSERVICEPROTODUNE_H
