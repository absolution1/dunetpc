////////////////////////////////////////////////////////////////////////////////
/// \file DUNEGeometryHelper_service.cc
///
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// class header
#include "dune/Geometry/DUNEGeometryHelper.h"

// LArSoft libraries
#include "Geometry/GeometryCore.h"
#include "Geometry/ChannelMapAlg.h"
#include "dune/Geometry/ChannelMap35Alg.h"
#include "dune/Geometry/ChannelMap35OptAlg.h"
#include "dune/Geometry/ChannelMapAPAAlg.h"

// C/C++ standard libraries
#include <string>


namespace dune
{

  DUNEGeometryHelper::DUNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & )
  :  fPset( pset ),
     fChannelMap()
  {}

  void DUNEGeometryHelper::doConfigureChannelMapAlg
    (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
  {
    fChannelMap.reset();
    
    std::string const detectorName = geom->DetectorName();
    
    //
    // DUNE 35t prototype
    //
    if (( detectorName.find("dune35t") != std::string::npos )
      || ( detectorName.find("lbne35t") != std::string::npos )
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
    else if (( detectorName.find("dune10kt") != std::string::npos ) 
      || ( detectorName.find("lbne10kt") != std::string::npos ))
    {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);
    }
    //
    // DUNE 34kt
    //
    else if (( detectorName.find("dune34kt") != std::string::npos )
      || ( detectorName.find("lbne34kt") != std::string::npos ))
    {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);
    }
    //
    // protoDUNE
    //
    else if (( detectorName.find("protodune") != std::string::npos )
      || ( detectorName.find("protolbne") != std::string::npos ))
    {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);
    }
    else {
      throw cet::exception("DUNEGeometryHelper")
        << "Unsupported geometry from input file: '" << detectorName << "'\n";
    }
    
    if ( fChannelMap )
    {
      geom->ApplyChannelMap(fChannelMap);
    }
  }
  
  
  DUNEGeometryHelper::ChannelMapAlgPtr_t
  DUNEGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(dune::DUNEGeometryHelper, geo::ExptGeoHelperInterface)
