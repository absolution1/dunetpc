////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapAPAAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELAPAMAPALG_H
#define GEO_CHANNELAPAMAPALG_H

#include <vector>
#include <stdint.h>

#include "Geometry/ChannelMapAlg.h"

namespace geo{

  class ChannelMapAPAAlg : public ChannelMapAlg{

  public:

    ChannelMapAPAAlg();
    ~ChannelMapAPAAlg();
    
    void                     Initialize(std::vector<geo::CryostatGeo*> & cgeo);
    void                     Uninitialize();
    std::vector<WireID>      ChannelToWire(uint32_t channel)           const;
    uint32_t                 Nchannels()                               const;
    WireID                   NearestWireID(const TVector3& worldPos,
					   unsigned int    PlaneNo,
					   unsigned int    TPCNo,
					   unsigned int    cstat)        const;
    uint32_t                 PlaneWireToChannel(unsigned int plane,
						unsigned int wire,
						unsigned int tpc,
						unsigned int cstat)    const;
   const View_t              View( uint32_t const channel )            const;
   const SigType_t           SignalType( uint32_t const channel )      const;
    
  private:
    
    unsigned int                                         fNcryostat;      ///< number of cryostats in the detector
    unsigned int                                         fNchannels;      ///< number of channels in the detector
    uint32_t                                             fTopChannel;     ///< book keeping highest channel #
    std::vector<unsigned int>                            fNTPC;           ///< number of TPCs in each cryostat

    // Assuming all APA's are identical
    std::vector< unsigned int >				 fWiresInPlane;
    unsigned int					 fPlanesPerAPA;   
    uint32_t					         fChannelsPerAPA;
    std::vector< unsigned int >				 nAnchoredWires;

    std::vector<std::vector<std::vector<unsigned int>>>  fWiresPerPlane;  ///< The number of wires in this plane 
                                                                          ///< in the heirachy

    std::vector<std::vector<std::vector<double>>> fFirstWireCenterY;
    std::vector<std::vector<std::vector<double>>> fFirstWireCenterZ;
    std::vector< double > fWirePitch;
    std::vector< double > fOrientation;
    std::vector< double > fTanOrientation; // to explore improving speed
    std::vector< double > fCosOrientation; // to explore improving speed



  };

}
#endif // GEO_CHANNELMAPAPAALG_H

