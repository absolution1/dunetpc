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
#include "larcore/Geometry/GeoObjectSorter.h"
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
 
  // Find the sorter.
  GeoObjectSorter* psort = nullptr;
  if ( fGeoSorterClass.size() ) {
    if ( fGeoSorterClass == "GeoObjectSorterAPA" ) {
      psort = new GeoObjectSorterAPA(pset);
    } else if ( fGeoSorterClass == "GeoObjectSorter35" ) {
      psort = new GeoObjectSorter35(pset);
    } else {
      throw cet::exception("DUNEGeometryHelper") << "Invalid sorter class name:" << fGeoSorterClass << "\n";
    }
  }

  // If class name is given use it.
  if ( fChannelMapClass.size() ) {
    if        ( fChannelMapClass == "Dune35tChannelMapAlg" ) {
      fChannelMap = std::make_shared<geo::Dune35tChannelMapAlg>(pset);
    } else if ( fChannelMapClass == "DuneApaChannelMapAlg" ) {
      geo::DuneApaChannelMapAlg* pmap = new geo::DuneApaChannelMapAlg(pset);
      pmap->setSorter(*psort);
      fChannelMap.reset(pmap);
    } else if ( fChannelMapClass == "ChannelMap35Alg" ) {
      fChannelMap = std::make_shared<geo::ChannelMap35Alg>(pset);
    } else if ( fChannelMapClass == "ChannelMap35OptAlg" ) {
      fChannelMap = std::make_shared<geo::ChannelMap35OptAlg>(pset);
    } else if ( fChannelMapClass == "ChannelMapAPAAlg" ) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(pset);
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
        fChannelMap = std::make_shared<geo::ChannelMap35OptAlg>(pset);
      } else {
        fChannelMap = std::make_shared<geo::ChannelMap35Alg>(pset);
      }

    // DUNE 10kt
    } else if (( detectorName.find("dune10kt") != std::string::npos ) 
          || ( detectorName.find("lbne10kt") != std::string::npos )) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(pset);

    // DUNE 10kt dual-phase
    } else if ( detectorName.find("dunedphase10kt") != std::string::npos ) {
      fChannelMap = std::make_shared<geo::ChannelMapCRMAlg>(pset);

    // DUNE 34kt
    } else if ( ( detectorName.find("dune34kt") != std::string::npos )
             || ( detectorName.find("lbne34kt") != std::string::npos ) ) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(pset);

    // protoDUNE
    } else if ( ( detectorName.find("protodune") != std::string::npos )
             || ( detectorName.find("protolbne") != std::string::npos ) ) {
      fChannelMap = std::make_shared<geo::ChannelMapAPAAlg>(pset);

    // LArND
    } else if ( detectorName.find("larnd") != std::string::npos ) {
      fChannelMap = std::make_shared<geo::ChannelMapStandardAlg>(pset);
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
