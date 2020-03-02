////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterAPA.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERAPA_H
#define GEO_GEOOBJECTSORTERAPA_H

#include <vector>

#include "fhiclcpp/fwd.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"

namespace geo{

  class GeoObjectSorterAPA : public GeoObjectSorter {
  public:

    GeoObjectSorterAPA(fhicl::ParameterSet const& p);

    void SortAuxDets        (std::vector<geo::AuxDetGeo>            & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo>   & adgeo)    const;
    void SortCryostats      (std::vector<geo::CryostatGeo>          & cgeo)     const;
    void SortTPCs           (std::vector<geo::TPCGeo>               & tgeo)     const;
    void SortPlanes         (std::vector<geo::PlaneGeo>             & pgeo,
                             geo::DriftDirection_t                    driftDir) const;
    void SortWires          (std::vector<geo::WireGeo>              & wgeo)     const;
    void SortOpDets         (std::vector<geo::OpDetGeo> & opdet) const;

  };

}

#endif // GEO_GEOOBJECTSORTERAPA_H
