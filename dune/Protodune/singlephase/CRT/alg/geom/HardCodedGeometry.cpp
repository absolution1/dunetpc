//File: HardCodedGeometry.cpp
//Brief: Implementation of a mapping from CRT electronics identifiers to 
//       sorting criteria for adjacency.  Uses hard-coded values from commissioning, 
//       so you should prefer something like an offline-geometry-based 
//       implementation. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

#include "HardCodedGeometry.h" //Header

namespace CRT
{
  HardCodedGeometry::HardCodedGeometry(): fModulesPerPlane(2), fStripsPerLayer(32)
  {
  }

  ModuleID HardCodedGeometry::doModuleID(const size_t module) const
  {
    const size_t nOrientations = 2; //2 orientations: horizontal and vertical
    const size_t modulesPerFrame = nOrientations*fModulesPerPlane;
    const size_t frameNum = module/modulesPerFrame;
    const FrameID frame(frameNum);
    const bool frameOrientation = frameNum % nOrientations; //0 is horizontal, 1 is vertical
                                                            //TODO: Maybe an Orientation type by itself would be useful
    const size_t frameLocal = module % modulesPerFrame;
    const PlaneID plane(frame, frameOrientation?frameLocal/fModulesPerPlane:!(frameLocal/fModulesPerPlane)); //TODO: Is this just an and?
    const CRT::ModuleID modID(plane, frameLocal % nOrientations);
    return modID;
  }

  StripID HardCodedGeometry::doStripID(const CRT::ModuleID module, const size_t channel) const
  {
    const LayerID layer(module, (channel > fStripsPerLayer));
    const CRT::StripID strip(layer, channel % fStripsPerLayer); //Implicitly defining the 0-31 layer of strips as the one 
                                                                //closest to the HV connect in data.  
    return strip;
  }
} 
