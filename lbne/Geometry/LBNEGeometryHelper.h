////////////////////////////////////////////////////////////////////////////////
/// \file LBNEGeometryHelper.h
/// \brief Geometry helper service for LBNE geometries. 
/// 
/// Handles LBNE-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class
///
/// \verion $Id
/// \author rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef LBNE_ExptGeoHelperInterface_h
#define LBNE_ExptGeoHelperInterface_h

#include "Geometry/ExptGeoHelperInterface.h"

#include <memory>
#include <vector>

// Forward declarations
//
class TString;

namespace geo
{
  class ChannelMapAlg;
  class CryostaGeo;
  class ExptGeoHelperInterface;
}

namespace geo
{
  class ChannelMapAlg;
}

// Declaration
//
namespace lbne
{
  class LBNEGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
  
    LBNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry &reg );
    ~LBNEGeometryHelper() throw();

    // Public interface for ExptGeoHelperInterface (for reference purposes)
    //
    // Configure and initialize the channel map.
    //
    // void  ConfigureChannelMapAlg( const TString & detectorName, 
    //                               fhicl::ParameterSet const & sortingParam,
    //                               std::vector<geo::CryostatGeo*> & c );
    //
    // Returns null pointer if the initialization failed
    // NOTE:  the sub-class owns the ChannelMapAlg object
    //
    // std::shared_ptr<const geo::ChannelMapAlg> & GetChannelMapAlg() const;
  
  private:
    
    void  doConfigureChannelMapAlg( const TString & detectorName,
                                    fhicl::ParameterSet const & sortingParam,
                                    std::vector<geo::CryostatGeo*> & c ) override;
    std::shared_ptr<const geo::ChannelMapAlg> doGetChannelMapAlg() const override;
    
    fhicl::ParameterSet const & fPset;
    art::ActivityRegistry & fReg;
    std::shared_ptr<geo::ChannelMapAlg> fChannelMap;
  
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(lbne::LBNEGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // LBNE_ExptGeoHelperInterface_h
