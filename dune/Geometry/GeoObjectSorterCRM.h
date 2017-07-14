////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterCRM.h
/// \brief Interface to algorithm class for sorting dual-phase CRM of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov vgalymov@ipnl.in2p3.fr
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERCRM_H
#define GEO_GEOOBJECTSORTERCRM_H

#include <vector>

#include "larcorealg/Geometry/GeoObjectSorter.h"

namespace geo{

  class GeoObjectSorterCRM : public GeoObjectSorter {

  public:

    GeoObjectSorterCRM(fhicl::ParameterSet const& p);
    ~GeoObjectSorterCRM();

    void SortAuxDets        (std::vector<geo::AuxDetGeo*>          & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo)   const;
    void SortCryostats      (std::vector<geo::CryostatGeo*>        & cgeo)     const;
    void SortTPCs     	    (std::vector<geo::TPCGeo*>      	   & tgeo)     const;
    void SortPlanes   	    (std::vector<geo::PlaneGeo*>    	   & pgeo,	      
		      	     geo::DriftDirection_t     	     const & driftDir) const;
    void SortWires    	    (std::vector<geo::WireGeo*>     	   & wgeo)     const;
    
  private:
    
  };

}

#endif // GEO_GEOOBJECTSORTERCRM_H
