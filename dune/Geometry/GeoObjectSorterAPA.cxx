////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterAPA.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "dune/Geometry/GeoObjectSorterAPA.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

namespace geo{

  //----------------------------------------------------------------------------
  static bool sortAuxDetAPA(const AuxDetGeo* ad1, const AuxDetGeo* ad2)
  {
  
    // NOTE: This initial sorting is arbitrary, there
    //       are no APAAlg users using AuxDet yet

    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    ad1->LocalToWorld(local, xyz1);
    ad2->LocalToWorld(local, xyz2);

    // First sort all AuxDets into same-y groups
    if(xyz1[1] < xyz2[1]) return true;
 
    // Within a same-y group, sort TPCs into same-z groups
    if(xyz1[1] == xyz2[1] && xyz1[2] < xyz2[2]) return true;
 
    // Within a same-z, same-y group, sort TPCs according to x
    if(xyz1[2] == xyz2[2] && xyz1[1] == xyz2[1] && xyz1[0] < xyz2[0]) return true;      

    // none of those are true, so return false
    return false;

  }

  //----------------------------------------------------------------------------
  static bool sortAuxDetSensitiveAPA(const AuxDetSensitiveGeo* ad1, const AuxDetSensitiveGeo* ad2)
  {
  
    // NOTE: This initial sorting is arbitrary, there
    //       are no APAAlg users using AuxDet yet

    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    ad1->LocalToWorld(local, xyz1);
    ad2->LocalToWorld(local, xyz2);

    // First sort all AuxDets into same-y groups
    if(xyz1[1] < xyz2[1]) return true;
 
    // Within a same-y group, sort TPCs into same-z groups
    if(xyz1[1] == xyz2[1] && xyz1[2] < xyz2[2]) return true;
 
    // Within a same-z, same-y group, sort TPCs according to x
    if(xyz1[2] == xyz2[2] && xyz1[1] == xyz2[1] && xyz1[0] < xyz2[0]) return true;      

    // none of those are true, so return false
    return false;

  }

  //----------------------------------------------------------------------------
  // Define sort order for cryostats in APA configuration
  //   same as standard
  static bool sortCryoAPA(const CryostatGeo* c1, const CryostatGeo* c2)
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.}; 
    c1->LocalToWorld(local, xyz1);
    c2->LocalToWorld(local, xyz2);

    return xyz1[0] < xyz2[0];   
  }


  //----------------------------------------------------------------------------
  // Define sort order for tpcs in APA configuration.
  static bool sortTPCAPA(const TPCGeo* t1, const TPCGeo* t2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    // The goal is to number TPCs first in the x direction so that,
    // in the case of APA configuration, TPCs 2c and 2c+1 make up APA c.
    // then numbering will go in y then in z direction.

    // First sort all TPCs into same-z groups
    if(xyz1[2] < xyz2[2]) return true;
 
    // Within a same-z group, sort TPCs into same-y groups
    if(xyz1[2] == xyz2[2] && xyz1[1] < xyz2[1]) return true;
 
    // Within a same-z, same-y group, sort TPCs according to x
    if(xyz1[2] == xyz2[2] && xyz1[1] == xyz2[1] && xyz1[0] < xyz2[0]) return true;
 
    // none of those are true, so return false
    return false;
  }


 //----------------------------------------------------------------------------
  // Define sort order for planes in APA configuration
  //   same as standard, but implemented differently
  static bool sortPlaneAPA(const PlaneGeo* p1, const PlaneGeo* p2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    p1->LocalToWorld(local, xyz1);
    p2->LocalToWorld(local, xyz2);

    return xyz1[0] > xyz2[0];
  }


  //----------------------------------------------------------------------------
  bool sortWireAPA(WireGeo* w1, WireGeo* w2){
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};

    w1->GetCenter(xyz1); w2->GetCenter(xyz2);

    // we want the wires to be sorted such that the smallest corner wire
    // on the readout end of a plane is wire zero, with wire number
    // increasing away from that wire. 

    // if a top TPC, count from top down
    if( xyz1[1]>0 && xyz1[1] > xyz2[1] ) return true;

    // if a bottom TPC, count from bottom up
    if( xyz1[1]<0 && xyz1[1] < xyz2[1] ) return true;

    // sort the all vertical wires to increase in +z direction
    if( xyz1[1] == xyz2[1] && xyz1[2] < xyz2[2] ) return true;

    return false;
  }


  //----------------------------------------------------------------------------
  GeoObjectSorterAPA::GeoObjectSorterAPA(fhicl::ParameterSet const&)
  {
  }

  //----------------------------------------------------------------------------
  GeoObjectSorterAPA::~GeoObjectSorterAPA()
  {
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterAPA::SortAuxDets(std::vector<geo::AuxDetGeo*> & adgeo) const
  {
    std::sort(adgeo.begin(), adgeo.end(), sortAuxDetAPA);
    
    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterAPA::SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adgeo) const
  {
    std::sort(adgeo.begin(), adgeo.end(), sortAuxDetSensitiveAPA);
    
    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterAPA::SortCryostats(std::vector<geo::CryostatGeo*> & cgeo) const
  {
    std::sort(cgeo.begin(), cgeo.end(), sortCryoAPA);
    
    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterAPA::SortTPCs(std::vector<geo::TPCGeo*>  & tgeo) const
  {
    
    std::sort(tgeo.begin(), tgeo.end(), sortTPCAPA);

    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterAPA::SortPlanes(std::vector<geo::PlaneGeo*> & pgeo,
				      geo::DriftDirection_t  const& driftDir) const
  {
    // sort the planes to increase in drift direction
    // The drift direction has to be set before this method is called.  It is set when
    // the CryostatGeo objects are sorted by the CryostatGeo::SortSubVolumes method
    if     (driftDir == geo::kPosX) std::sort(pgeo.rbegin(), pgeo.rend(), sortPlaneAPA);
    else if(driftDir == geo::kNegX) std::sort(pgeo.begin(),  pgeo.end(),  sortPlaneAPA);
    else if(driftDir == geo::kUnknownDrift)
      throw cet::exception("TPCGeo") << "Drift direction is unknown, can't sort the planes\n";

    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterAPA::SortWires(std::vector<geo::WireGeo*> & wgeo) const
  {
    std::sort(wgeo.begin(), wgeo.end(), sortWireAPA);

    return;
  }

}
