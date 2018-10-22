//File: Geometry.cpp
//Brief: Interface that converts CRT electronics channels to 
//       sorting criteria for finding overlaps.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#include "Geometry.h" //Header

namespace CRT
{
  ModuleID Geometry::ModuleID(const size_t module) const
  {
    return doModuleID(module);
  }  

  StripID Geometry::StripID(const size_t module, const size_t channel) const
  {
    //I can do any common post-processing of results here
    return doStripID(ModuleID(module), channel);  
  }
} 
