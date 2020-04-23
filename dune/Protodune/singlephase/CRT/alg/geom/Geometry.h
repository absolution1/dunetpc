//File: Geometry.h
//Brief: A CRT::Geometry is a mapping from some channel numbering system 
//       to a labelling system for CRT detector components and a mapping 
//       from that labelling system to strip geometries.  It serves as a 
//       channel-map-like abstraction that doesn't need to know about the 
//       offline framework.  
//
//       A CRT::Geometry is just an interface to such a mapping, so it needs 
//       concrete implmentation(s) to be used.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRT_GEOMETRY_H
#define CRT_GEOMETRY_H

#include "dune/Protodune/singlephase/CRT/alg/geom/CRTID.h"

namespace CRT
{
  class Geometry
  {
    public:
      //TODO: Decide on a construction process for Geometry
      Geometry() = default; 
      virtual ~Geometry() = default;

      //Mapping from (module, channel) pair, like in CRT::Trigger, to strip identifier in 
      //"geometry space".  StripID is designed to be more convenient for finding overlaps 
      //between CRT channels.  
      CRT::ModuleID ModuleID(const size_t module) const;
      CRT::StripID StripID(const size_t module, const size_t channel) const;
      //TODO: Some struct that encapsulates a (module, channel) pair to make it harder to 
      //      get this wrong in the future.  Probably entails an update to CRT::Trigger 
      //      to do this right.  

      //Mapping from StripID to StripGeo.  StripGeos should give all of 
      //the information needed to turn a set of 4 overlapping StripIDs into 
      //a 3D position in the offline coordinate system.  
      //CRT::StripGeo& StripGeo(const CRT::StripID& id) const; //TODO: Write StripGeo

    private:
      //Private implementation.  Interface implementations must implement these pure 
      //virtual functions.  
      virtual CRT::StripID doStripID(const CRT::ModuleID module, const size_t channel) const = 0;
      virtual CRT::ModuleID doModuleID(const size_t module) const = 0;

      //virtual CRT::StripGeo& StripGeo(const CRT::StripID& id) const = 0; //TODO: Write StripGeo
  };
}

#endif //CRT_GEOMETRY_H  
