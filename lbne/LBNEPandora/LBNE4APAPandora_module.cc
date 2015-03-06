/**
 *  @file   lbnecode/LBNEPandora/LBNE4APAPandora_module.cc
 *
 *  @brief  Producer module for LBNE 4APA detector.
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local includes
#include "LArPandoraParticleCreator.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  LBNE4APAPandora class
 */
class LBNE4APAPandora : public LArPandoraParticleCreator
{
public: 

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    LBNE4APAPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~LBNE4APAPandora();

private:

    unsigned int GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const;
    void ConfigurePandoraGeometry() const;

    bool            m_useLeftVolume;      ///<
    bool            m_useRightVolume;       ///<
};

DEFINE_ART_MODULE(LBNE4APAPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// LArSoft includes
#include "Geometry/Geometry.h"

// Local includes (LArContent) 
#include "LArContent.h"

// Local includes (LArPandora)
#include "LBNE4APAPseudoLayerPlugin.h"
#include "LBNE4APATransformationPlugin.h"
#include "LBNE4APAGeometryHelper.h"

namespace lar_pandora {

LBNE4APAPandora::LBNE4APAPandora(fhicl::ParameterSet const &pset) : LArPandoraParticleCreator(pset)
{
    m_useLeftVolume = pset.get<bool>("UseLeftVolume",true);
    m_useRightVolume = pset.get<bool>("UseRightVolume",true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LBNE4APAPandora::~LBNE4APAPandora()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LBNE4APAPandora::ConfigurePandoraGeometry() const
{
    mf::LogDebug("LArPandora") << " *** LBNE4APAPandora::ConfigurePandoraGeometry(...) *** " << std::endl;

    // Identify the Geometry and load the plugins
    art::ServiceHandle<geo::Geometry> theGeometry;

    if (std::string::npos == theGeometry->DetectorName().find("lbne10kt"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " LBNE4APAPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end(); 
        pIter != pIterEnd; ++pIter)
    {
        const unsigned int      volumeID = pIter->first;
        const pandora::Pandora *pPandora = pIter->second;

        const bool isForward((0 == volumeID) ? true : false); // ATTN: Sign of rotation matrix is taken from Volume ID
    
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, 
            new LBNE4APAPseudoLayerPlugin));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, 
            new LBNE4APATransformationPlugin(isForward)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LBNE4APAPandora::GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const
{    
    const LBNE4APAGeometryHelper::LBNE4APAVolume volumeID(LBNE4APAGeometryHelper::GetVolumeID(cstat, tpc));

    if (LBNE4APAGeometryHelper::kLeftVolume == volumeID) 
    {
        if (m_useLeftVolume) 
            return 0;
    }

    if (LBNE4APAGeometryHelper::kRightVolume == volumeID) 
    {
        if (m_useRightVolume) 
            return 1;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lar_pandora
