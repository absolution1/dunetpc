#ifndef DETPEDESTALSERVICEDUNE_CC
#define DETPEDESTALSERVICEDUNE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/ParameterSet.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Providers/DetPedestalRetrievalAlg.h"
#include "dune/RunHistory/DetPedestalDUNE.h"
#include "Geometry/Geometry.h"

namespace lariov{

  class DetPedestalServiceDUNE : public IDetPedestalService {
  
    public:
    
    DetPedestalServiceDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~DetPedestalServiceDUNE(){}
    
    void reconfigure(fhicl::ParameterSet const& pset) {
      fProvider.Configure(pset);
    }

    void preBeginRun(const art::Run& run) {
      fProvider.Update(run.id().run());
    }
    
    private:
    
      const IDetPedestalProvider& DoGetPedestalProvider() const override {
        return fProvider;
      }    
    
    dune::DetPedestalDUNE fProvider;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::DetPedestalServiceDUNE, lariov::IDetPedestalService, LEGACY)
      

namespace lariov{

  DetPedestalServiceDUNE::DetPedestalServiceDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("DetPedestalDUNE"))
  {

    //register callback to update local database cache before each event is processed
    //reg.sPreProcessEvent.watch(&SIOVDetPedestalService::PreProcessEvent, *this);
    reg.sPreBeginRun.watch(this, &DetPedestalServiceDUNE::preBeginRun);

    art::ServiceHandle<geo::Geometry> geom;
    std::string detName = geom->DetectorName();
    std::string d35tName = "dune35t";
    
    if (detName.find(std::string(d35tName)) != std::string::npos)
      fProvider.SetDetName(d35tName);
    
  }

}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::DetPedestalServiceDUNE, lariov::IDetPedestalService)

#endif
