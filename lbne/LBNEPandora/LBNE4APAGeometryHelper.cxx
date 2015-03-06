/**
 *  @file  lbnecode/LBNEPandora/LBNE4APAGeometryHelper.cxx
 *
 *  @brief helper function for LBNE 4APA geometry
 *
 */

#include "LBNE4APAGeometryHelper.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"

#include "cetlib/exception.h"

#include <limits>
#include <iostream>

namespace lar_pandora {

LBNE4APAGeometryHelper::LBNE4APAVolume LBNE4APAGeometryHelper::GetVolumeID(const unsigned int cstat, const unsigned int tpc)
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc,cstat));

    // Left drift volume: negative drift direction
    if (theTpcGeo.DriftDirection() == geo::kNegX)
    {
        return LBNE4APAGeometryHelper::kLeftVolume;
    }

    // Right drift volume: positive drift direction
    if (theTpcGeo.DriftDirection() == geo::kPosX)
    {
        return LBNE4APAGeometryHelper::kRightVolume;
    }

    throw cet::exception("LArPandora") << " LBNE4APAGeometryHelper::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
