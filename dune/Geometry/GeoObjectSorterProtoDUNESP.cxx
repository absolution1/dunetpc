////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterProtoDUNESP.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

//LArSoft includes
#include "dune/Geometry/GeoObjectSorterProtoDUNESP.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
//#include "larcoreobj/SimpleTypesAndConstants/geo_vectors_utils.h" //Would provide geo::vect::middlePoint()

//c++ includes
#include <numeric> //std::accumulate was moved here in c++14 

namespace
{
  template <class ITERATOR>
  void SortInSpiral(const ITERATOR begin, const ITERATOR end)
  {
    const auto center = std::accumulate(begin, end, geo::Vector_t(0., 0., 0.), [](geo::Vector_t sum, const auto geo) 
                                                                               { return sum+geo->GetCenter(); });

    const geo::Vector_t vertical(0., 1., 0.);

    //Sort by angle between vector from center of all AuxDetGeos to center of each AuxDetGeo 
    //and vertical unit vector.  If this is really a bottleneck, one could probably write some 
    //multi-return-statement version that's more efficient.  
    std::sort(begin, end, [&center, &vertical](const auto first, const auto second)
                          {
                            const auto firstVec = first->GetCenter()-center;
                            const auto secondVec = second->GetCenter()-center;
                             
                            return firstVec.Dot(vertical) < secondVec.Dot(vertical);
                          });
  }

  //Return the first ITERATOR to pointer-to-AuxDetGeo that is more than zTolerance distant from the one before 
  //it and end if no element between begin and end is far enough from its neighbor.  
  template <class ITERATOR>
  ITERATOR FindFirstDownstream(const ITERATOR begin, const ITERATOR end, const double zTolerance)
  {
    if(begin == end) return end; //In a container of size 1, there is no way to use this algorithm to 
                                 //find out whether this range contains an upstream or a downstream volume

    ITERATOR prev = begin;
    for(auto geo = std::next(begin); geo != end; ++geo)
    {
      if(fabs(((*geo)->GetCenter() - (*prev)->GetCenter()).Z()) > zTolerance) return geo;
      prev = geo;
    }

    return end;
  }
}

namespace geo{

  //----------------------------------------------------------------------------
  GeoObjectSorterProtoDUNESP::GeoObjectSorterProtoDUNESP(fhicl::ParameterSet const& p): GeoObjectSorterAPA(p)
  {
  }

  //----------------------------------------------------------------------------
  GeoObjectSorterProtoDUNESP::~GeoObjectSorterProtoDUNESP()
  {
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterProtoDUNESP::SortAuxDets(std::vector<geo::AuxDetGeo*> & adgeo) const
  {
    //Since this only happens once and I expect to have to sort only 32 AuxDetGeos, I've written 
    //an algorithm that prioritizes readability over efficiency.  The result is that it loops over 
    //each volume at least 3 times when it might only need to loop once.  If this is a performance bottleneck, 
    //consider a "return early" sorting strategy or something based on caching information.  

    //First, find the thickest AuxDetGeo in z
    const double thickest = 800.; //0; //Some upstream CRTs are up to 7m apart, but upstream is at least 10m from downstream.  
    //If anyone wants to use this algorithm in a more general setting where my assumptions about rough positions of 
    //the CRTs might not be true, the loop below might provide a good starting point.  
    /*for(const auto& geo: adgeo)
    {
      const auto& center = geo->GetCenter();

      //Assume that all elements of adgeo are represented well by box geometries.
      const LocalPoint_t maxZ(0., 0., geo.GetLength()), maxY(0., geo.GetHalfHeight(), 0.), 
                         maxX(std::max(geo.GetHalfWidth1(), geo.GetHalfWidth2()), 0., 0.);
      //TODO: Translate each of these to global coordinates and set deltaZ to largest z distance from center

      const auto deltaZ = ;
      if(deltaZ > thickest) thickest = deltaZ;
    }*/

    //Next, separate AuxDetGeos into upstream and downstream.  
    std::sort(adgeo.begin(), adgeo.end(), [&thickest](const auto& first, const auto& second) 
                                          { return fabs(first->GetCenter().Z()-second->GetCenter().Z()) < thickest; });

    //Then, find an iterator to the first downstream AuxDetGeo
    const auto firstDownstream = ::FindFirstDownstream(adgeo.begin(), adgeo.end(), thickest);

    //Last, sort upstream and downstream containers in a spiral about their center.  
    ::SortInSpiral(adgeo.begin(), firstDownstream); //upstream
    ::SortInSpiral(firstDownstream, adgeo.end()); //downstream
    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterProtoDUNESP::SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adgeo) const
  {
    std::sort(adgeo.begin(), adgeo.end(), [](const auto& first, const auto& second)
                                          {
                                            static constexpr double tolerance = 0.5;
                                            //Since I'm working with the CRT, I know that AuxDetSensitiveGeos will be about 
                                            //1cm "thick" in z and 5cm in x/y once they are rotated to the global 
                                            //frame.  Since I'm working with double-precision numbers, define tolerance of 
                                            //comparisons so that it is about half of the smallest dimension of a CRT strip.


                                            //Calculate all of the differences between centers I need at once.  
                                            const auto diff = second->GetCenter()-first->GetCenter();
                                            if(fabs(diff.Z()) > tolerance) return diff.Z() > 0; //Sort into layers in z
                                            
                                            //Each CRT module has only strips in one direction.  I don't know how this module is 
                                            //oriented, so arbitrarily choose to try to sort in x first.  If you suspect that this 
                                            //is a problem, you could try converting one detector's center to the other detector's 
                                            //local coordinate system.  
                                            if(fabs(diff.X()) > tolerance) return diff.X() > 0; //Sort from beam-left to beam-right
                                            return diff.Y() > 0; //Sort from top to bottom 
                                          });
    
    return;
  }
}
