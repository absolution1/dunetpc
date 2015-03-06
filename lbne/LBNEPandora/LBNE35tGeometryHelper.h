/**
 *  @file  lbnecode/LBNEPandora/LBNE35tGeometryHelper.h
 *
 *  @brief helper function for LBNE 35t geometry
 *
 */
#ifndef LBNE_35T_PANDORA_HELPER_H
#define LBNE_35T_PANDORA_HELPER_H

namespace lar_pandora 
{

class LBNE35tGeometryHelper 
{
public:

    enum LBNE35tVolume
    {
        kShortVolume = 0,
        kLongVolume = 1,
        kUnknownVolume = 2
    };

    /**
     *  @brief Assign a drift volume ID based on cryostate and TPC
     *
     *  @param cstat the cryostat
     *  @param tpc the tpc 
     */
     static LBNE35tGeometryHelper::LBNE35tVolume GetVolumeID(const unsigned int cstat, const unsigned int tpc);
};

} // namespace lar_pandora

#endif //  LBNE_35T_PANDORA_HELPER_H
