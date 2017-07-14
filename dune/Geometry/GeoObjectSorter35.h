////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter35.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTER35_H
#define GEO_GEOOBJECTSORTER35_H

#include <vector>
#include <string>
#include "larcorealg/Geometry/GeoObjectSorter.h"

namespace geo{

  class GeoObjectSorter35 : public GeoObjectSorter {

  public:

    GeoObjectSorter35(fhicl::ParameterSet const& p);
    ~GeoObjectSorter35();

    void SortAuxDets        (std::vector<geo::AuxDetGeo*>            & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*>   & adgeo)    const;
    void SortCryostats      (std::vector<geo::CryostatGeo*>          & cgeo)     const;
    void SortTPCs     	    (std::vector<geo::TPCGeo*>      	     & tgeo)     const;
    void SortPlanes   	    (std::vector<geo::PlaneGeo*>    	     & pgeo,	      
		      	     geo::DriftDirection_t     	        const& driftDir) const;
    void SortWires    	    (std::vector<geo::WireGeo*>     	     & wgeo)     const;
    
  private:
    
    std::string fDetVersion;  ///< string of the detector version
  };

}

#endif // GEO_GEOOBJECTSORTER35_H
