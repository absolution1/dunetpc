////////////////////////////////////////////////////////////////////////                                                                                                                                                                    
/// \file  CRTSorter.h                                                                                                                                                                                                     
/// \brief A function to sort the CRTs for ProtoDUNESP
/// \version $Id:  $                                                                                                                                                                                                                        
/// \author  e.tyley@sheffield.ac.uk                                                                                                                                                                                                      
////////////////////////////////////////////////////////////////////////                                                                                                                                                                     
#ifndef CRT_CRTSORTER_H
#define CRT_CRTSORTER_H

#include <vector>
#include "dune/Geometry/GeoObjectSorterAPA.h"

namespace CRTSorter{

  std::vector<size_t> Mapping(std::vector<geo::AuxDetGeo*>& adgeo); //Overload for validation via geometry sorting

  std::vector<size_t> Mapping(std::vector<const geo::AuxDetGeo*>& adgeo); //Overload for kludge after geometry sorting where only 
                                                                          //const geo::AuxDetGeo& is available.  You're probably using 
                                                                          //this one.  
}

#endif // CRT_CRTSORTER_H











