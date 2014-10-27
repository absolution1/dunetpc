////////////////////////////////////////////////////////////////////////
// $Id: GeometryTest35_module.cc,v 1.1 2011/02/17 01:45:48 brebel Exp $
//
//
// geometry unit tests
//
// tylerdalion@gmail.com
//
////////////////////////////////////////////////////////////////////////

#ifndef GEO_GEOMETRYTEST35_H
#define GEO_GEOMETRYTEST35_H
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

// ROOT includes
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TNtuple.h"
#include "TGeoManager.h"
#include "TStopwatch.h"
#include "TMath.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/OpDetGeo.h"
#include "Geometry/geo.h"
#include "SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDAnalyzer.h"

namespace geo { class Geometry; }

///tracking algorithms
namespace geo {
  class GeometryTest35 : public art::EDAnalyzer {
  public:
    explicit GeometryTest35(fhicl::ParameterSet const& pset);
    virtual ~GeometryTest35();

    virtual void analyze(art::Event const&) {}
    virtual void beginJob();

  private:


  };
}

namespace geo{

  //......................................................................
  GeometryTest35::GeometryTest35(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset)
  {
  }

  //......................................................................
  GeometryTest35::~GeometryTest35()
  {
  }

  //......................................................................
  void GeometryTest35::beginJob()
  {
    //art::ServiceHandle<geo::Geometry> geom;

    std::cout << "35t specific testing...\n";
    mf::LogVerbatim("GeometryTest35") << "35t specific testing...\n";

    try{

      LOG_DEBUG("GeometryTest35") << "test ...";

      LOG_DEBUG("GeometryTest35") << "complete.";

    }
    catch (cet::exception &e) {
      mf::LogWarning("GeometryTest35") << "exception caught: \n" << e;
    }
    
    return;
  }


}//end namespace


namespace geo{

  DEFINE_ART_MODULE(GeometryTest35)

}

#endif
