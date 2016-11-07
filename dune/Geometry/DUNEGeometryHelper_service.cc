////////////////////////////////////////////////////////////////////////////////
/// \file DUNEGeometryHelper_service.cc
///
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// class header
#include "dune/Geometry/DUNEGeometryHelper.h"

// LArSoft libraries
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/ChannelMapAlg.h"
#include "larcore/Geometry/ChannelMapStandardAlg.h"
#include "dune/Geometry/Dune35tChannelMapAlg.h"
#include "dune/Geometry/DuneApaChannelMapAlg.h"
#include "dune/Geometry/ChannelMap35Alg.h"
#include "dune/Geometry/ChannelMap35OptAlg.h"
#include "dune/Geometry/ChannelMapAPAAlg.h"
#include "dune/Geometry/ChannelMapCRMAlg.h"

// C/C++ standard libraries
#include <string>

using std::string;
using dune::DUNEGeometryHelper;

//**********************************************************************

DUNEGeometryHelper::DUNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & )
: fPset( pset ),
  fChannelMap() {
  pset.get_if_present<string>("ChannelMapClass", fChannelMapClass);
}

//**********************************************************************

void DUNEGeometryHelper::
doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom) {
  fChannelMap.reset();
 
  // If class name is given use it.
  if ( fChannelMapClass.size() ) {
    if        ( fChannelMapClass == "Dune35tChannelMapAlg" ) {
      fChannelMap = std::make_shared<geo::Dune35tChannelMapAlg>(sortingParameters);
    } else if ( fChannelMapClass == "DuneApaChannelMapAlg" ) {
      fChannelMap = std::make_shared<geo::DuneApaChannelMapAlg>(sortingParameters);
    } else if ( fChannelMapClass == "ChannelMap35Alg" ) {
      fChannelMap = std::make_shared<geo::ChannelMap35Alg>(sortingParameters);
    } else if ( fChannelMapClass == "ChannelMap35OptAlg" ) {
      fChannelMap = std::make_shared<geo::ChannelMap35OptAlg>(sortingParameters);
    } else if ( fChannelMapClass == "ChannelMapAPAAlg" ) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);
    } else {
      throw cet::exception("DUNEGeometryHelper") << "Invalid channel map class name:" << fChannelMapClass << "\n";
    }

  // Otherwise derive class name from detector name.
  } else {
    std::string const detectorName = geom->DetectorName();

    // DUNE 35t prototype
    if ( ( detectorName.find("dune35t") != std::string::npos ) ||
         ( detectorName.find("lbne35t") != std::string::npos ) ) {
      std::string const detectorVersion = sortingParameters.get< std::string >("DetectorVersion");
      if (( detectorVersion.find("v3") != std::string::npos )
       || ( detectorVersion.find("v4") != std::string::npos )
       || ( detectorVersion.find("v5") != std::string::npos )
       || ( detectorVersion.find("v6") != std::string::npos )) {
        fChannelMap = std::make_shared<geo::ChannelMap35OptAlg>(sortingParameters);
      } else {
        fChannelMap = std::make_shared<geo::ChannelMap35Alg>(sortingParameters);
      }

    // DUNE 10kt
    } else if (( detectorName.find("dune10kt") != std::string::npos ) 
          || ( detectorName.find("lbne10kt") != std::string::npos )) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);

    // DUNE 10kt dual-phase
    } else if ( detectorName.find("dunedphase10kt") != std::string::npos ) {
      fChannelMap = std::make_shared<geo::ChannelMapCRMAlg>(sortingParameters);

    // DUNE 34kt
    } else if ( ( detectorName.find("dune34kt") != std::string::npos )
             || ( detectorName.find("lbne34kt") != std::string::npos ) ) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);

    // protoDUNE
    } else if ( ( detectorName.find("protodune") != std::string::npos )
             || ( detectorName.find("protolbne") != std::string::npos ) ) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(sortingParameters);

    // LArND
    } else if ( detectorName.find("larnd") != std::string::npos ) {
      fChannelMap = std::make_shared<geo::ChannelMapStandardAlg>(sortingParameters);
    }
    else {
      throw cet::exception("DUNEGeometryHelper") << "Unsupported detector: '" << detectorName << "'\n";
    }    
  }
  
  if ( fChannelMap ) {
    geom->ApplyChannelMap(fChannelMap);
  }
}  
//**********************************************************************
  
DUNEGeometryHelper::ChannelMapAlgPtr_t DUNEGeometryHelper::doGetChannelMapAlg() const {
  return fChannelMap;
}

//**********************************************************************
  
DEFINE_ART_SERVICE_INTERFACE_IMPL(dune::DUNEGeometryHelper, geo::ExptGeoHelperInterface)

//**********************************************************************
