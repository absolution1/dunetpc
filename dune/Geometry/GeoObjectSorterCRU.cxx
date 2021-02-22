////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterCRU.cxx
/// \brief Interface to algorithm class for sorting VD CRUs of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author vgalymov@ipnl.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "dune/Geometry/GeoObjectSorterCRU.h"
#include "dune/Geometry/OpDetSorter.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include <array>
#include <string>
#include <cmath> // std::abs()


// comparison functions for sorting various geo objects
namespace geo { 
  
  namespace CRU {
    //----------------------------------------------------------------------------
    // Define sort order for auxdet in VD configuration
    static bool sortAuxDet(const AuxDetGeo& ad1, const AuxDetGeo& ad2)
    {

      // sort based off of GDML name, assuming ordering is encoded
      std::string ad1name = (ad1.TotalVolume())->GetName();
      std::string ad2name = (ad2.TotalVolume())->GetName();

      // assume volume name is "volAuxDet##"
      int ad1Num = atoi( ad1name.substr( 9, ad1name.size()).c_str() );
      int ad2Num = atoi( ad2name.substr( 9, ad2name.size()).c_str() );

      return ad1Num < ad2Num;

    }

    //----------------------------------------------------------------------------
    // Define sort order for auxsensitve det in VD configuration
    static bool sortAuxDetSensitive(const AuxDetSensitiveGeo& ad1, const AuxDetSensitiveGeo& ad2)
    {

      // sort based off of GDML name, assuming ordering is encoded
      std::string ad1name = (ad1.TotalVolume())->GetName();
      std::string ad2name = (ad2.TotalVolume())->GetName();

      // assume volume name is "volAuxDetSensitive##"
      int ad1Num = atoi( ad1name.substr( 9, ad1name.size()).c_str() );
      int ad2Num = atoi( ad2name.substr( 9, ad2name.size()).c_str() );

      return ad1Num < ad2Num;

    }

    //----------------------------------------------------------------------------
    // Define sort order for cryostats in VD configuration
    static bool sortCryo(const CryostatGeo& c1, const CryostatGeo& c2)
    {
      double xyz1[3] = {0.}, xyz2[3] = {0.};
      double local[3] = {0.};
      c1.LocalToWorld(local, xyz1);
      c2.LocalToWorld(local, xyz2);

      return xyz1[0] < xyz2[0];
    }


    //----------------------------------------------------------------------------
    // Define sort order for tpcs (CRUs) in VD configuration
    static bool sortTPC(const TPCGeo& t1, const TPCGeo& t2)
    {
      double xyz1[3] = {0.};
      double xyz2[3] = {0.};
      double local[3] = {0.};
      t1.LocalToWorld(local, xyz1);
      t2.LocalToWorld(local, xyz2);

      // First sort all TPCs into same-z groups
      if(xyz1[2]<xyz2[2]) return true;

      // Within a same-z group, sort TPCs into same-y groups
      if(xyz1[2] == xyz2[2] && xyz1[1] < xyz2[1]) return true;

      return false;
    }


    //----------------------------------------------------------------------------
    // Define sort order for planes in VD configuration
    static bool sortPlane(const PlaneGeo& p1, const PlaneGeo& p2)
    {
      double xyz1[3] = {0.};
      double xyz2[3] = {0.};
      double local[3] = {0.};
      p1.LocalToWorld(local, xyz1);
      p2.LocalToWorld(local, xyz2);

      // drift direction is negative, plane number increases in drift direction
      return xyz1[0] > xyz2[0];
    }


    //----------------------------------------------------------------------------
    bool sortWire(WireGeo const& w1, WireGeo const& w2){
      // wire sorting algorithm
      // z_low -> z_high
      // for same-z group 
      //  the sorting either from bottom (y) left (z) corner up for strip angles > pi/2
      //  and horizontal wires
      //  or from top corner to bottom otherwise
      

      std::array<double, 3> xyz1, xyz2;
      w1.GetCenter(xyz1.data()); w2.GetCenter(xyz2.data());
      
      // always sort in the increasing z
      if( xyz1[2] < xyz2[2] ){
	return true;
      }
      
      if( w1.isHorizontal() && xyz1[1] < xyz2[1] ){
	return true;
      }
      
      if( w1.ThetaZ(true) > 90.0 && xyz1[1] < xyz2[1] ){
	return true;
      }
      
      if( w1.ThetaZ(true) < 90.0 && xyz1[1] > xyz2[1] ){
	return true;
      }
      
      return false;
      
      // sort 
      
      /*
      // same z
      if ( std::abs(xyz1[2] - xyz2[2]) < 1.0E-4 ){
	if (      std::abs(xyz1[1] - xyz2[1]) < 1.0E-4 ) // same y
	  return xyz1[0] < xyz2[0];
	else if ( std::abs(xyz1[0] - xyz2[0]) < 1.0E-4 ) // same x
	  return xyz1[1] < xyz2[1];
      }
      
      return xyz1[2] < xyz2[2];
      */
    }

  } // namespace geo::CRU

  //----------------------------------------------------------------------------
  GeoObjectSorterCRU::GeoObjectSorterCRU(fhicl::ParameterSet const&)
  {}

  //----------------------------------------------------------------------------
  void GeoObjectSorterCRU::SortAuxDets(std::vector<geo::AuxDetGeo> & adgeo) const
  {
    std::sort(adgeo.begin(), adgeo.end(), CRU::sortAuxDet);
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterCRU::SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo> & adsgeo) const
  {
    std::sort(adsgeo.begin(), adsgeo.end(), CRU::sortAuxDetSensitive);
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterCRU::SortCryostats(std::vector<geo::CryostatGeo> & cgeo) const
  {
    std::sort(cgeo.begin(), cgeo.end(), CRU::sortCryo);
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterCRU::SortTPCs(std::vector<geo::TPCGeo>  & tgeo) const
  {
    std::sort(tgeo.begin(), tgeo.end(), CRU::sortTPC);
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterCRU::SortPlanes(std::vector<geo::PlaneGeo> & pgeo,
                                           geo::DriftDirection_t  const driftDir) const
  {
    // sort the planes to increase in drift direction
    // The drift direction has to be set before this method is called.  It is set when
    // the CryostatGeo objects are sorted by the CryostatGeo::SortSubVolumes method
    if     (driftDir == geo::kPosX) std::sort(pgeo.rbegin(), pgeo.rend(), CRU::sortPlane);
    else if(driftDir == geo::kNegX) std::sort(pgeo.begin(),  pgeo.end(),  CRU::sortPlane);
    else if(driftDir == geo::kUnknownDrift)
      throw cet::exception("TPCGeo") << "Drift direction is unknown, can't sort the planes\n";
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterCRU::SortWires(std::vector<geo::WireGeo> & wgeo) const
  {
    std::sort(wgeo.begin(), wgeo.end(), CRU::sortWire);
  }

  //----------------------------------------------------------------------------
 void GeoObjectSorterCRU::SortOpDets         (std::vector<geo::OpDetGeo> & opdet) const {
   std::sort(opdet.begin(), opdet.end(), sortorderOpDet);
 }

}
