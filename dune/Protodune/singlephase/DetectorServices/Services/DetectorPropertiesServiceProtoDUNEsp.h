////////////////////////////////////////////////////////////////////////
// DetectorPropertiesServiceProtoDUNEsp.h
//
// Service interface for DetectorProperties functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORPROPERTIESSERVICEPROTODUNESP_H
#define DETECTORPROPERTIESSERVICEPROTODUNESP_H
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"

#include "art/Framework/Principal/Event.h"
//#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "dune/Protodune/singlephase/DetectorServices/Providers/DetectorPropertiesProtoDUNEsp.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


#include "art/Persistency/Provenance/ScheduleContext.h"


///General LArSoft Utilities
namespace spdp{
  
  /**
   * @brief "Standard" implementation of DetectorProperties service
   * 
   * This class wraps DetectorPropertiesStandard provider into a art service.
   * It delivers the provider via the standard interface:
   *     
   *     detinfo::DetectorProperties const* detprop
   *       = art::ServiceHandle<detinfo::DetectorPropertiesStandard>()
   *       ->provider();
   *     
   * or, using the standard interface in "CoreUtils/ServiceUtil.h":
   *     
   *     auto const* detprop
   *       = lar::providerFrom<detinfo::DetectorPropertiesStandard>();
   *     
   * In addition to the functionality of the provider, this service allows
   * to read the configuration from the input file, inherited from a previous
   * run.
   * 
   * Configuration parameters
   * -------------------------
   * 
   * This service passes the whole configuration down to its service provider,
   * but it also reacts to:
   * - *InheritNumberTimeSamples* (boolean; default: false): if true, the
   *   configuration database in the ROOT input file is queried and if a
   *   configuration for this service is found, it's used instead of the
   *   one from the current FHiCL configuration
   * 
   */
  
  class DetectorPropertiesServiceProtoDUNEsp : public detinfo::DetectorPropertiesService {
    public:
      
      // the following is currently not used for validation,
      // but only for documentation
      struct ServiceConfiguration_t {
        
        // service-specific configuration
        fhicl::Atom<bool> InheritNumberTimeSamples {
          fhicl::Name("InheritNumberTimeSamples"),
          fhicl::Comment(""),
          false /* default value */
        };
        
        // provider configuration
        spdp::DetectorPropertiesProtoDUNEsp::Configuration_t ProviderConfiguration;
        
      }; // ServiceConfiguration_t
      
      
      // this enables art to print the configuration help:
      using Parameters = art::ServiceTable<ServiceConfiguration_t>;
      
      DetectorPropertiesServiceProtoDUNEsp(fhicl::ParameterSet const& pset,
                                art::ActivityRegistry& reg);
      void   reconfigure(fhicl::ParameterSet const& pset);
      void   preProcessEvent(const art::Event& evt, art::ScheduleContext);
      void   postOpenFile(const std::string& filename);
      void   preOpenFile(const std::string& filename);
      void preBeginRun(const art::Run& run);
      
      virtual const provider_type* provider() const override { return fProp.get();}
      
    private:
      std::unique_ptr<spdp::DetectorPropertiesProtoDUNEsp> fProp;
      fhicl::ParameterSet   fPS;       ///< Original parameter set.
      bool isNewRun;
      bool fInheritNumberTimeSamples; ///< Flag saying whether to inherit NumberTimeSamples
      
      bool isDetectorPropertiesServiceProtoDUNEsp(const fhicl::ParameterSet& ps) const;
      
    }; // class DetectorPropertiesService
} //namespace detinfo
DECLARE_ART_SERVICE_INTERFACE_IMPL(spdp::DetectorPropertiesServiceProtoDUNEsp, detinfo::DetectorPropertiesService, LEGACY)
#endif // DETECTORPROPERTIESSERVICESTANDARD_H