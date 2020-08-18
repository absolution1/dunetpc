////////////////////////////////////////////////////////////////////////
// \file LifetimeCalibServiceProtoDUNE.h
//
// \brief header of service for storing/accessing lifetime calibration corrections for ProtoDUNE
//
// \author wwu@fnal.gov
// \date May 19, 2020
// 
////////////////////////////////////////////////////////////////////////
#ifndef LIFETIMECALIBSERVICEPROTODUNE_H
#define LIFETIMECALIBSERVICEPROTODUNE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "dune/Calib/LifetimeCalibProtoDUNE.h"
#include "dune/CalibServices/LifetimeCalibService.h"

namespace calib{
  class LifetimeCalibServiceProtoDUNE : public LifetimeCalibService {
    public:
      
      LifetimeCalibServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset)  override;
      void   preBeginRun(const art::Run& run); 
      //void   preBeginRun(const art::Event& evt);


      virtual provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<calib::LifetimeCalibProtoDUNE> fProp;

    }; // class LifetimeCalibServiceProtoDUNE
} //namespace calib
DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::LifetimeCalibServiceProtoDUNE, calib::LifetimeCalibService, LEGACY)
#endif // LIFETIMECALIBSERVICEPROTODUNE_H
