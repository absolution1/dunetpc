/**
 * @file   StandardGeometryHelper.h
 * @brief  Geometry helper service for DUNE geometries
 * @author rs@fnal.gov
 * 
 * Handles DUNE-specific information for the generic Geometry service
 * within LArSoft. Derived from the ExptGeoHelperInterface class.
 */

#ifndef LBNE_ExptGeoHelperInterface_h
#define LBNE_ExptGeoHelperInterface_h

#include "Geometry/ExptGeoHelperInterface.h"

#include <memory>

// Forward declarations
namespace geo
{
  class ChannelMapAlg;
  class GeometryCore;
}

// Declaration
//
namespace lbne
{
  class LBNEGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
    
    LBNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry &reg );
    
    /*
      Public interface for ExptGeoHelperInterface (for reference purposes)
      
      Configure, initialize and return the channel map:
      
      void ConfigureChannelMapAlg
        (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom);
      
      Returns null pointer if the initialization failed:
      
      ChannelMapAlgPtr_t GetChannelMapAlg() const;
    */
  
  private:
    
    virtual void doConfigureChannelMapAlg
      (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
      override;
    virtual ChannelMapAlgPtr_t doGetChannelMapAlg() const override;
    
    fhicl::ParameterSet const & fPset;
    std::shared_ptr<geo::ChannelMapAlg> fChannelMap;
    
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(lbne::LBNEGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // LBNE_ExptGeoHelperInterface_h
