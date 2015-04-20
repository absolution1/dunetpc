/**
 *  @file  lbnecode/LBNEPandora/LBNE35tGeometryHelper.cxx
 *
 *  @brief helper function for LBNE 35t geometry
 *
 */

#include "lbne/LBNEPandora/LBNE35tGeometryHelper.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"

#include "cetlib/exception.h"

#include <limits>
#include <iostream>

namespace lar_pandora {

LBNE35tGeometryHelper::LBNE35tVolume LBNE35tGeometryHelper::GetVolumeID(const unsigned int cstat, const unsigned int tpc)
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc,cstat));

    // Long drift volume: negative drift direction, odd TPC numbers (1 == tpc%2)
    if (theTpcGeo.DriftDirection() == geo::kNegX)
    {
        return LBNE35tGeometryHelper::kLongVolume;
    }

    // Short drift volume: positive drift direction, even TPC numbers (0 == tpc%2)
    if (theTpcGeo.DriftDirection() == geo::kPosX)
    {
        return LBNE35tGeometryHelper::kShortVolume;
    }

    throw cet::exception("LArPandora") << " LBNE35tGeometryHelper::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
