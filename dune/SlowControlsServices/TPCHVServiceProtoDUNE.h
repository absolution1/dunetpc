////////////////////////////////////////////////////////////////////////
// \file TPCHVServiceProtoDUNE.h
//
// \brief header of service for storing/accessing (x,y,z) calibration corrections for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef TPCHVSERVICEPROTODUNE_H
#define TPCHVSERVICEPROTODUNE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "dune/SlowControls/SlowControlsProtoDUNE.h"
#include "dune/SlowControlsServices/TPCHVService.h"

namespace slowctrls{

  class TPCHVServiceProtoDUNE : public TPCHVService {
    public:
      
      TPCHVServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset)  override;
      //      void   preBeginRun(const art::Run& run);

      virtual provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<slowctrls::SlowControlsProtoDUNE> fProp;

    }; // class TPCHVServiceProtoDUNE
} //namespace slowctrl
DECLARE_ART_SERVICE_INTERFACE_IMPL(slowctrls::TPCHVServiceProtoDUNE, slowctrls::TPCHVService, LEGACY)
#endif // TPCHVSERVICEPROTODUNE_H
