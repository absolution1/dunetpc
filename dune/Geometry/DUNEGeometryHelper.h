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
}

// Declaration
//
namespace dune {

class DUNEGeometryHelper : public geo::ExptGeoHelperInterface {
public:

  explicit DUNEGeometryHelper(fhicl::ParameterSet const& pset);

private:

  std::unique_ptr<geo::ChannelMapAlg>
  doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                           std::string const& detectorName) const override;

  // FCL params.
  std::string fChannelMapClass;   // Name of the class that does the channel mapping.
  std::string fGeoSorterClass;    // Name of the class that does the sorting.
};

}  // end dune namespace

DECLARE_ART_SERVICE_INTERFACE_IMPL(dune::DUNEGeometryHelper,
                                   geo::ExptGeoHelperInterface,
                                   SHARED)

#endif // DUNE_ExptGeoHelperInterface_h
