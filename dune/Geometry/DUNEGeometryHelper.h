/**
 * @file   StandardGeometryHelper.h
 * @brief  Geometry helper service for DUNE geometries
 * @author rs@fnal.gov
 * 
 * Handles DUNE-specific information for the generic Geometry service
 * within LArSoft. Derived from the ExptGeoHelperInterface class.
 */

// Modified Oct 2016 by David Adams.
// Add param ChannelMapClass to explicitly specify the mapping class.

#ifndef DUNE_ExptGeoHelperInterface_h
#define DUNE_ExptGeoHelperInterface_h

#include "larcore/Geometry/ExptGeoHelperInterface.h"

#include <memory>

// Forward declarations
namespace geo
{
  class ChannelMapAlg;
  class GeometryCore;
}

#if defined __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-private-field"
#endif

// Declaration
//
namespace dune {

class DUNEGeometryHelper : public geo::ExptGeoHelperInterface {

public:
    
  DUNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry &reg );
    
  /*
    Public interface for ExptGeoHelperInterface (for reference purposes)
    
    Configure, initialize and return the channel map:
      
    void ConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom);
      
    Returns null pointer if the initialization failed:
      
      ChannelMapAlgPtr_t GetChannelMapAlg() const;
  */
  
private:
    
  virtual void
  doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom) override;
  virtual ChannelMapAlgPtr_t doGetChannelMapAlg() const override;
    
  fhicl::ParameterSet const & fPset; 
  std::shared_ptr<geo::ChannelMapAlg> fChannelMap;

  // FCL params.
  std::string fChannelMapClass;   // Name of the class that does the channel mapping.
  std::string fGeoSorterClass;    // Name of the class that does the sorting.
    
};

}  // end dune namespace
#if defined __clang__
  #pragma clang diagnostic pop
#endif

DECLARE_ART_SERVICE_INTERFACE_IMPL(dune::DUNEGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // DUNE_ExptGeoHelperInterface_h
