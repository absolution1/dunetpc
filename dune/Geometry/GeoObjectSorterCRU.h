////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterCRU.h
/// \brief Interface to algorithm class for sorting VD CRUs of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author vgalymov@ipnl.in2p3.fr
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERCRU_H
#define GEO_GEOOBJECTSORTERCRU_H

#include <vector>

#include "fhiclcpp/fwd.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"

namespace geo{

  class GeoObjectSorterCRU : public GeoObjectSorter {
  public:

    GeoObjectSorterCRU(fhicl::ParameterSet const& p);

    void SortAuxDets        (std::vector<geo::AuxDetGeo>         & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo>& adsgeo)   const;
    void SortCryostats      (std::vector<geo::CryostatGeo>       & cgeo)     const;
    void SortTPCs           (std::vector<geo::TPCGeo>            & tgeo)     const;
    void SortPlanes         (std::vector<geo::PlaneGeo>          & pgeo,
                             geo::DriftDirection_t                 driftDir) const;
    void SortWires          (std::vector<geo::WireGeo>           & wgeo)     const;
    void SortOpDets         (std::vector<geo::OpDetGeo> & opdet) const;

  };

}

#endif // GEO_GEOOBJECTSORTERCRM_H
