//File: HardCodedGeometry.h
//Brief: CRT geometry from hard-coded values.  Converts (module, channel) pairs into 
//       sorting criteria for identifying overlapping hits in the detector.  Prefer 
//       values from the offline framework's geometry service when possible.
//Author: Andrew Olivier aolivier@ur.rochester.edu 

#ifndef CRT_HARDCODEDGEOMETRY_H
#define CRT_HARDCODEDGEOMETRY_H
#include "dunetpc/dune/Protodune/singlephase/CRT/alg/geom/Geometry.h" //Interface that this header describes an implementation for.  

namespace CRT
{
  class HardCodedGeometry: public Geometry
  {
    public: 
      HardCodedGeometry();
      virtual ~HardCodedGeometry() = default;

      //Public interface provided by Geometry
    private:
      virtual CRT::ModuleID doModuleID(const size_t module) const override;
      virtual CRT::StripID doStripID(const CRT::ModuleID module, const size_t channel) const override;

      //Data used for mapping calculation
      const size_t fModulesPerPlane;
      const size_t fStripsPerLayer;
  };
}

#endif //CRT_HARDCODEDGEOMETRY_H
