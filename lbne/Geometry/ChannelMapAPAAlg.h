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

#include "Geometry/ChannelMapAlg.h"

namespace geo{

  class ChannelMapAPAAlg : public ChannelMapAlg{

  public:

    ChannelMapAPAAlg();
    ~ChannelMapAPAAlg();
    
    void                     Initialize(std::vector<geo::CryostatGeo*> const& cgeo);
    void                     Uninitialize();
    std::vector<WireID>      ChannelToWire(unsigned int channel)    const;
    unsigned int             Nchannels()                            const;
    unsigned int             NearestWire(const TVector3& worldPos,
                                         unsigned int    PlaneNo,
                                         unsigned int    TPCNo,
                                         unsigned int    cstat)     const;
    unsigned int             PlaneWireToChannel(unsigned int plane,
						unsigned int wire,
						unsigned int tpc,
						unsigned int cstat) const;
    
  private:
    
    unsigned int                                         fNcryostat;      ///< number of cryostats in the detector
    unsigned int                                         fNchannels;      ///< number of channels in the detector
    unsigned int                                         fTopChannel;     ///< book keeping highest channel #
    std::vector<unsigned int>                            fNTPC;           ///< number of TPCs in each cryostat

	// Assuming all APA's are identical
    std::vector< unsigned int >				 fWiresInPlane;
    unsigned int					 fPlanesPerAPA;   
    unsigned int					 fChannelsPerAPA;
    std::vector< unsigned int >				 nAnchoredWires;
    std::vector<std::vector<std::vector<unsigned int>>>	 fAPAs;		  ///< list of TPC pairs that make APAs


    std::vector<std::vector<std::vector<unsigned int>>>  fWiresPerPlane;  ///< The number of wires in this plane 
                                                                          ///< in the heirachy

    std::vector<std::vector<std::vector<double>>> fFirstWireCenterY;
    std::vector<std::vector<std::vector<double>>> fFirstWireCenterZ;
    std::vector< double > fWirePitch;
    std::vector< double > fOrientation;
    std::vector< double > fTanOrientation; // to explore improving speed



  };

}
#endif // GEO_CHANNELMAPAPAALG_H

