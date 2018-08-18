////////////////////////////////////////////////////////////////////////////////
/// \file DUNEGeometryHelper_service.cc
///
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// class header
#include "dune/Geometry/DUNEGeometryHelper.h"

// LArSoft libraries
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"
#include "dune/Geometry/DuneApaChannelMapAlg.h"
#include "dune/Geometry/ChannelMap35Alg.h"
#include "dune/Geometry/ChannelMap35OptAlg.h"
#include "dune/Geometry/ChannelMapAPAAlg.h"
#include "dune/Geometry/ChannelMapCRMAlg.h"
#include "dune/Geometry/ProtoDUNEChannelMapAlg.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "dune/Geometry/GeoObjectSorterAPA.h"
#include "dune/Geometry/GeoObjectSorter35.h"

// C/C++ standard libraries
#include <string>

using std::string;
using geo::GeoObjectSorter;
using geo::GeoObjectSorterAPA;
using geo::GeoObjectSorter35;
using dune::DUNEGeometryHelper;

//**********************************************************************

DUNEGeometryHelper::DUNEGeometryHelper(const fhicl::ParameterSet& pset, art::ActivityRegistry&)
: fPset( pset ),
  fChannelMap() {
  pset.get_if_present<string>("ChannelMapClass", fChannelMapClass);
  pset.get_if_present<string>("GeoSorterClass", fGeoSorterClass);
}

//**********************************************************************

void DUNEGeometryHelper::
doConfigureChannelMapAlg(fhicl::ParameterSet const& pset, geo::GeometryCore* geom) {
  fChannelMap.reset();
 
  // The APA sorting algorithm requires special handling.
  // This flag is set if that map is used.
  bool useApaMap = false;
  bool is35t = false;

  // If class name is given use it.
  if ( fChannelMapClass.size() ) {
    if ( fChannelMapClass == "DuneApaChannelMapAlg" ) {
      useApaMap = true;
    } else if ( fChannelMapClass == "ChannelMap35Alg" ) {
      fChannelMap = std::make_shared<geo::ChannelMap35Alg>(pset);
    } else if ( fChannelMapClass == "ChannelMap35OptAlg" ) {
      fChannelMap = std::make_shared<geo::ChannelMap35OptAlg>(pset);
    } else if ( fChannelMapClass == "ChannelMapAPAAlg" ) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(pset);
    } else if ( fChannelMapClass == "ProtoDUNEChannelMapAlg" ) {
      fChannelMap = std::make_shared<geo::ProtoDUNEChannelMapAlg>(pset);
    } else {
      throw cet::exception("DUNEGeometryHelper") << "Invalid channel map class name:" << fChannelMapClass << "\n";
    }

  // Otherwise derive class name from detector name.
  } else {
    std::string const detectorName = geom->DetectorName();

    // DUNE 35t prototype
    if ( ( detectorName.find("dune35t") != std::string::npos ) ||
         ( detectorName.find("lbne35t") != std::string::npos ) ) {
      std::string const detectorVersion = pset.get< std::string >("DetectorVersion");
      if (( detectorVersion.find("v3") != std::string::npos )
       || ( detectorVersion.find("v4") != std::string::npos )
       || ( detectorVersion.find("v5") != std::string::npos )
       || ( detectorVersion.find("v6") != std::string::npos )) {
        useApaMap = true;
        is35t = true;
      } else {
        fChannelMap = std::make_shared<geo::ChannelMap35Alg>(pset);
      }

    // DUNE 10kt
    } else if (( detectorName.find("dune10kt") != std::string::npos ) 
          || ( detectorName.find("lbne10kt") != std::string::npos )) {
      useApaMap = true;

    // DUNE 10kt dual phase
    } else if ( detectorName.find("dunedphase10kt") != std::string::npos ) {
      fChannelMap = std::make_shared<geo::ChannelMapCRMAlg>(pset);

    // protoDUNE 6x6x6 dual phase
    } else if ( detectorName.find("protodunedphase") != std::string::npos ) {
      fChannelMap = std::make_shared<geo::ChannelMapCRMAlg>(pset);

    // 3x1x1 dual phase
    } else if ( detectorName.find("3x1x1dphase") != std::string::npos ) {
      fChannelMap = std::make_shared<geo::ChannelMapCRMAlg>(pset);

    // DUNE 34kt
    } else if ( ( detectorName.find("dune34kt") != std::string::npos )
             || ( detectorName.find("lbne34kt") != std::string::npos ) ) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(pset);

    // protoDUNE
    } else if ( ( detectorName.find("protodune") != std::string::npos )
             || ( detectorName.find("protolbne") != std::string::npos ) ) {
      fChannelMap = std::make_shared<geo::ProtoDUNEChannelMapAlg>(pset);

    // LArND
    } else if ( detectorName.find("larnd") != std::string::npos ) {
      fChannelMap = std::make_shared<geo::ChannelMapStandardAlg>(pset);
    }
    else {
      throw cet::exception("DUNEGeometryHelper") << "Unsupported detector: '" << detectorName << "'\n";
    }    
  }
  
  if ( useApaMap ) {
    // Find the sorter.
    bool useApaSort = false;
    bool use35tSort = false;
    GeoObjectSorter* psort = nullptr;
    if ( fGeoSorterClass.size() ) {
      if ( fGeoSorterClass == "GeoObjectSorterAPA" ) {
        useApaSort = true;
      } else if ( fGeoSorterClass == "GeoObjectSorter35" ) {
        use35tSort = true;
      }
    } else {
      if ( is35t ) use35tSort = true;
      else useApaSort = true;
    }
    // Construct the sorter.
    if ( useApaSort ) {
      psort = new GeoObjectSorterAPA(pset);
    } else if ( use35tSort ) {
      psort = new GeoObjectSorter35(pset);
    }
    // Create channel map and set sorter.
    geo::DuneApaChannelMapAlg* pmap = new geo::DuneApaChannelMapAlg(pset);
    pmap->setSorter(*psort);
    fChannelMap.reset(pmap);
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
