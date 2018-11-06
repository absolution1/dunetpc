////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterProtoDUNESP.h
/// \brief A GeoObjectSorterAPA that also knows how to sort AuxDets 
///        for conveniently working with the ProtoDUNE-SP Cosmic Ray 
///        Tagger (CRT).
///
/// \version $Id:  $
/// \author  aolivier@ur.rochester.edu
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERPROTODUNESP_H
#define GEO_GEOOBJECTSORTERPROTODUNESP_H

#include <vector>

#include "dune/Geometry/GeoObjectSorterAPA.h"

namespace geo{

  class GeoObjectSorterProtoDUNESP : public GeoObjectSorterAPA {

  public:

    GeoObjectSorterProtoDUNESP(fhicl::ParameterSet const& p);
    ~GeoObjectSorterProtoDUNESP();

    void SortAuxDets        (std::vector<geo::AuxDetGeo*>            & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*>   & adgeo)    const;
        
  private:
    
  };

}

#endif // GEO_GEOOBJECTSORTERPROTODUNESP_H
