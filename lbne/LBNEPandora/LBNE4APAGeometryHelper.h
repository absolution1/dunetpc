/**
 *  @file  lbnecode/LBNEPandora/LBNE4APAGeometryHelper.h
 *
 *  @brief helper function for LBNE 4APA geometry
 *
 */
#ifndef LBNE_4APA_PANDORA_HELPER_H
#define LBNE_4APA_PANDORA_HELPER_H

namespace lar_pandora 
{

class LBNE4APAGeometryHelper 
{
public:

    enum LBNE4APAVolume
    {
        kLeftVolume = 0,
        kRightVolume = 1,
        kUnknownVolume = 2
    };

    /**
     *  @brief Assign a drift volume ID based on cryostate and TPC
     *
     *  @param cstat the cryostat
     *  @param tpc the tpc 
     */
     static LBNE4APAGeometryHelper::LBNE4APAVolume GetVolumeID(const unsigned int cstat, const unsigned int tpc);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_LBNE_4APA_PANDORA_HELPER_H
