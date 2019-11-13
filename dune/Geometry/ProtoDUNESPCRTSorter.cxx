///////////////////////////////////////////////////////////////////////                                                       
/// \file  CRTSorter.cxx    
/// \brief A function to sort the CRTs for ProtoDUNESP   
/// \version $Id:  $         
/// \author  e.tyley@sheffield.ac.uk  
////////////////////////////////////////////////////////////////////////   
//LArSoft includes                                                                                                            

#include "larcorealg/Geometry/AuxDetGeo.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "larcoreobj/SimpleTypesAndConstants/geo_vectors_utils.h" 

#include "dune/Geometry/ProtoDUNESPCRTSorter.h"
                        

//c++ includes                                                                                                                 
#include <numeric> //std::accumulate was moved here in c++14                                                                   

namespace
{
  // function to work out the edge where the CRT reads out, the outside edge of the CRT                                 
  geo::Vector_t FindEdge(const geo::AuxDetGeo& det,const geo::Vector_t Center)
  {
    // Work out whether the CRT is vertical or horizontal                                                                      
    const bool yPos = (det.toWorldCoords(geo::AuxDetGeo::LocalPoint_t(0, 0, det.Length()/2.)).Y() - det.GetCenter().Y()
                       < det.toWorldCoords(geo::AuxDetGeo::LocalPoint_t(0, 0, det.Length()/2.)).X() - det.GetCenter().X());

    double xCo = det.GetCenter().X();
    double yCo = det.GetCenter().Y();
    double len = det.Length();

    // Work out whoch way to push the CRT position to make it represent the correct edge                                       
    // if yPos(Reads out in X) push in x to the edge of the detector                                                           
 
    if (yPos)
      {
        if (det.GetCenter().X()>Center.X()) xCo+=len/2.;
        else xCo-=len/2.;
      }
    else
      {
        if (det.GetCenter().Y()>Center.Y()) yCo+=len/2.;
        else yCo-=len/2.;
      }
    //Return the vector of the edge of the CRT                                                                                 

    return geo::Vector_t(xCo,yCo,det.GetCenter().Z());
  }
  
 
  template <class ITERATOR>
  void SortInSpiral(const ITERATOR begin, const ITERATOR end)
  {

    int counter = end-begin; //Calculate number of CRT panels                                             
    //Calculate the center point by averaging the position of the panels                                
 
    const geo::Vector_t center = std::accumulate(begin, end, geo::Vector_t(0., 0., 0.), [&counter](geo::Vector_t sum, const auto geo)
                                                 {return sum+geo->GetCenter()/counter;});

    //const geo::Vector_t vertical(0., 1., 0.);                                                                               
    //std::cout<<"Center: "<<center.X()<<" X, "<<center.Y()<<"Y and "<<center.Z()<<"Z With counter: "<<counter<<std::endl; 
    //Sort by angle between vector from center of all AuxDetGeos to edge of each AuxDetGeo                   
    //If this is really a bottleneck, one could probably write some                                         
    //multi-return-statement version that's more efficient.                                                                    
  
    std::sort(begin, end, [&center](const geo::AuxDetGeo* first,const geo::AuxDetGeo* second)
              {
                //Find the readout edge of each CRT                                                                            
		geo::Vector_t firstEdge = FindEdge(*first,center);
		geo::Vector_t secondEdge = FindEdge(*second,center);

                //Calculate the vector from the center to each edge                                                            
		geo::Vector_t firstVec = firstEdge-center;
		geo::Vector_t secondVec = secondEdge-center;

                //Calculate the angle Phi (Angle in X-Y Plane)                                                                 
                double firstPhi = firstVec.Phi();
                double secondPhi = secondVec.Phi();


                //We want to start at verticle (pi) and go in a spiral from there                                                              //so this means pushing the angles of >pi to the end of the array                                              
                //by subtracting 2*pi                                                                                          
                                                                                                              
                if (firstPhi>1.57) firstPhi-=6.28;
                if (secondPhi>1.57) secondPhi-=6.28;

		if ((*first).GetCenter().Z() <0.0){
		  return firstPhi > secondPhi;
		} else {
                  return firstPhi < secondPhi;
		}

              });

    // loop through the sorted vector and print out to check it is working as intended                                                                                                                                    
    /*
    for(ITERATOR geo = (begin); geo != end; ++geo)                                                                                                                                                                    
    {
      geo::Vector_t firstEdge = FindEdge(*(*geo),center);                                                                       
      geo::Vector_t firstVec = firstEdge-center;                                                                               
      double firstPhi = firstVec.Phi();                                                                                        
      
      if (firstPhi>1.57) {                                                                                                     
	firstPhi-=3.14;                                                                                                        
	std::cout<<"Angle Phi -2*pi: "<<firstPhi<<std::endl;                                                                  
	}                                                                                                                       
                                                                                                             
      else{std::cout<<"Angle Phi: "<<firstPhi<<std::endl;}                                                                     
      //std::cout<<"Angle Phi: "<<firstPhi<<std::endl;                      
      //std::cout<<(*geo)->GetCenter().X()<<"X, "<<(*geo)->GetCenter().Y()<<"Y and "<<(*geo)->GetCenter().Z()<<"Z"<<std::endl;
      //std::cout<<"Edge: "<<firstEdge.X()<<"X, "<<firstEdge.Y()<<"Y and "<<firstEdge.Z()<<"Z"<<std::endl;                 
      //std::cout<<"Vec:  "<<firstVec.X()<<"X, "<<firstVec.Y()<<"Y and "<<firstVec.Z()<<"Z"<<std::endl;                 
      }
    */
  }


  //Return the first ITERATOR to pointer-to-AuxDetGeo that is more than zTolerance distant from the one before               
  //it and end if no element between begin and end is far enough from its neighbor.                                           

  template <class ITERATOR>
  ITERATOR FindFirstDownstream(const ITERATOR begin, const ITERATOR end, const double zTolerance)
  {
    if(begin == end) return end; //In a container of size 1, there is no way to use this algorithm to                          
                                 //find out whether this range contains an upstream or a downstream volume                     

    //std::cout<<"test2"<<std::endl;
    ITERATOR prev = begin;
    for(ITERATOR geo = std::next(begin); geo != end; ++geo)
      {
        //std::cout<<fabs(((*geo)->GetCenter() - (*prev)->GetCenter()).Z())<<std::endl;                                                                                                                                                     
        //std::cout<<((*prev)->GetCenter()).Z()<<std::endl;                                                                                                                                                                               
  
        if(fabs(((*geo)->GetCenter() - (*prev)->GetCenter()).Z()) > zTolerance)
	  {
	    return geo;
	  }
	prev = geo;
      }
    return end;
  }
  
  template <class CONTAINER>
  std::vector<size_t> map(const CONTAINER& sorted, const CONTAINER& unsorted)
  {
    std::vector<size_t> mymap;
    for(size_t i = 0; i < sorted.size(); ++i)
      {
	auto ThisAdgeo = sorted[i];
        for(size_t j = 0; j < unsorted.size(); ++j)
          {
	    auto ThisUnsorted  = unsorted[j];

            if (ThisAdgeo->Name() == ThisUnsorted->Name())
              {
                mymap.push_back(j);
		//std::cout<<ThisAdgeo->Name()<<" Sorted CRT "<< i <<" Unsorted CRT "<< j <<std::endl;
              }
          }
      }
    return mymap;
  }
  
  //Wrapper so that I don't have to care whether this is a std:vector<const geo::AuxDetGeo*> or std::vector<geo::AuxDetGeo*>.  
  //I'm not going to change the AuxDetGeos, so I just require that iterators passed dereference to a const AuxDetGeo*.  
  template <class CONTAINER>
  std::vector<size_t> doMapping(CONTAINER& unsorted)
  {
    std::cout<<"called"<<std::endl;
    auto adgeo = unsorted;
    //auto unsorted = adgeo;
    const double thickest = 800.;

    //std::cout<<"test"<<std::endl;
    std::sort(adgeo.begin(), adgeo.end(), [](const auto& first, const auto& second)
              { return (first->GetCenter()).Z() < (second->GetCenter()).Z() ;} );

    //Then, find an iterator to the first downstream AuxDetGeo                                                          
    //std::cout<<"test"<<std::endl;
    const auto firstDownstream = ::FindFirstDownstream(adgeo.begin(), adgeo.end(), thickest);
    std::cout<<"test"<<std::endl;

    //Last, sort upstream and downstream containers in a spiral about their center.    
    ::SortInSpiral(adgeo.begin(), firstDownstream); //upstream
    ::SortInSpiral(firstDownstream, adgeo.end()); //downstream
    
    std::vector<size_t> mapping = ::map(adgeo,unsorted); //TODO: Rewrite this to use iterators?
    std::cout<<"end"<<std::endl;
    return mapping;
  }
}

namespace CRTSorter
{
  std::vector<size_t> Mapping(std::vector<geo::AuxDetGeo*>& adgeo) //Overload for validation via geometry sorting
  {
    return doMapping(adgeo);
  }

  std::vector<size_t> Mapping(std::vector<const geo::AuxDetGeo*>& adgeo) //Overload for kludge after geometry sorting where only 
                                                                        //const geo::AuxDetGeo& is available
  {
    return doMapping(adgeo);
  }
}
