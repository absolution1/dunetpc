/**
 *  @file   lbnecode/LBNEPandora/LBNE35tParticleStitcher_module.cc
 *
 *  @brief  Stitching module for LBNE 35t detector.
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local includes
#include "PFParticleStitcher.h"

#include "cetlib/exception.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  LBNE35tParticleStitcher class
 */
class LBNE35tParticleStitcher : public PFParticleStitcher
{
public: 

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    LBNE35tParticleStitcher(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~LBNE35tParticleStitcher();

private:
    unsigned int GetVolumeID(const unsigned int cstat, const unsigned int tpc) const;

};

DEFINE_ART_MODULE(LBNE35tParticleStitcher)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "LBNE35tGeometryHelper.h"

namespace lar_pandora {

LBNE35tParticleStitcher::LBNE35tParticleStitcher(fhicl::ParameterSet const &pset) : PFParticleStitcher(pset)
{
  
}

//------------------------------------------------------------------------------------------------------------------------------------------

LBNE35tParticleStitcher::~LBNE35tParticleStitcher()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LBNE35tParticleStitcher::GetVolumeID(const unsigned int cstat, const unsigned int tpc) const
{
    LBNE35tGeometryHelper::LBNE35tVolume volumeID(LBNE35tGeometryHelper::GetVolumeID(cstat, tpc));

    if (LBNE35tGeometryHelper::kShortVolume == volumeID) 
        return 0;

    if (LBNE35tGeometryHelper::kLongVolume == volumeID) 
        return 1;

    throw cet::exception("LArPandora") << " LBNE35tParticleStitcher::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
