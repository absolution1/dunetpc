/**
 *  @file   lbnecode/LBNEPandora/LBNE35tPandora_module.cc
 *
 *  @brief  Producer module for LBNE 35t detector.
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local includes
#include "LArPandoraInterface/LArPandoraParticleCreator.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  LBNE35tPandora class
 */
class LBNE35tPandora : public LArPandoraParticleCreator
{
public: 

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    LBNE35tPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~LBNE35tPandora();

private:

    unsigned int GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const;
    void ConfigurePandoraGeometry() const;

    bool            m_useShortVolume;      ///<
    bool            m_useLongVolume;       ///<
};

DEFINE_ART_MODULE(LBNE35tPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// LArSoft includes
#include "Geometry/Geometry.h"

// Local includes (LArContent) 
#include "LArPandoraAlgorithms/LArContent.h"

// Local includes (LArPandora)
#include "LBNE35tPseudoLayerPlugin.h"
#include "LBNE35tTransformationPlugin.h"
#include "LBNE35tGeometryHelper.h"

namespace lar_pandora {

LBNE35tPandora::LBNE35tPandora(fhicl::ParameterSet const &pset) : LArPandoraParticleCreator(pset)
{
    m_useShortVolume = pset.get<bool>("UseShortVolume",true);
    m_useLongVolume = pset.get<bool>("UseLongVolume",true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LBNE35tPandora::~LBNE35tPandora()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LBNE35tPandora::ConfigurePandoraGeometry() const
{
    mf::LogDebug("LArPandora") << " *** LBNE35tPandora::ConfigurePandoraGeometry(...) *** " << std::endl;

    // Identify the Geometry and load the plugins
    art::ServiceHandle<geo::Geometry> theGeometry;

    if (std::string::npos == theGeometry->DetectorName().find("lbne35t"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " LBNE35tPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end(); 
        pIter != pIterEnd; ++pIter)
    {
        const unsigned int      volumeID = pIter->first;
        const pandora::Pandora *pPandora = pIter->second;

        const bool isForward((0 == volumeID) ? true : false); // ATTN: Sign of rotation matrix is taken from Volume ID
    
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, 
            new LBNE35tPseudoLayerPlugin));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, 
            new LBNE35tTransformationPlugin(isForward)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LBNE35tPandora::GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const
{    
    const LBNE35tGeometryHelper::LBNE35tVolume volumeID(LBNE35tGeometryHelper::GetVolumeID(cstat, tpc));

    if (LBNE35tGeometryHelper::kShortVolume == volumeID) 
    {
        if (m_useShortVolume) 
            return 0;
    }

    if (LBNE35tGeometryHelper::kLongVolume == volumeID) 
    {
        if (m_useLongVolume) 
            return 1;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lar_pandora
