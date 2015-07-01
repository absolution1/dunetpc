////////////////////////////////////////////////////////////////////////////////
/// \file LBNEGeometryHelper_service.cc
///
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// class header
#include "lbne/Geometry/LBNEGeometryHelper.h"

// LArSoft libraries
#include "Geometry/GeometryCore.h"
#include "Geometry/ChannelMapAlg.h"
#include "lbne/Geometry/ChannelMap35Alg.h"
#include "lbne/Geometry/ChannelMap35OptAlg.h"
#include "lbne/Geometry/ChannelMapAPAAlg.h"

// C/C++ standard libraries
#include <string>


namespace lbne
{

  LBNEGeometryHelper::LBNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & )
  :  fPset( pset ),
     fChannelMap()
  {}

  void LBNEGeometryHelper::doConfigureChannelMapAlg
    (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
  {
    fChannelMap.reset();
    
    std::string const detectorName = geom->DetectorName();
    
    //
    // DUNE 35t prototype
    //
    if (( detectorName.find("lbne35t") != std::string::npos )
      || ( detectorName.find("dune35t") != std::string::npos )
      )
    {
      std::string const detectorVersion
        = sortingParameters.get< std::string >("DetectorVersion");
      
      if (( detectorVersion.find("v3") != std::string::npos )
        || ( detectorVersion.find("v4") != std::string::npos )
        || ( detectorVersion.find("v5") != std::string::npos ))
        fChannelMap = std::make_shared<geo::ChannelMap35OptAlg>(sortingParameters);
      else
        fChannelMap = std::make_shared<geo::ChannelMap35Alg>(sortingParameters);
      
    }
    //
    // DUNE 10kt
    //
    else if (( detectorName.find("lbne10kt") != std::string::npos ) 
      || ( detectorName.find("dune10kt") != std::string::npos ))
    {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);
    }
    //
    // DUNE 34kt
    //
    else if (( detectorName.find("lbne34kt") != std::string::npos )
      || ( detectorName.find("dune34kt") != std::string::npos ))
    {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);
    }
    
    if ( fChannelMap )
    {
      geom->ApplyChannelMap(fChannelMap);
    }
  }
  
  
  LBNEGeometryHelper::ChannelMapAlgPtr_t
  LBNEGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(lbne::LBNEGeometryHelper, geo::ExptGeoHelperInterface)
