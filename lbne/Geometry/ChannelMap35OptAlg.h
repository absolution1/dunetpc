////////////////////////////////////////////////////////////////////////
/// \file  ChannelMap35OptAlg.h
/// \brief The class of 35t specific algorithms, optimized
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////
///
/// This class is starting as a copy of ChannelMap35Alg, plus one bug fix
/// in the loop that counts the number of anchored wires in an APA, or 
/// rather the number of channels per APA.
///
/// NOTE: Actual optimization still needs to be done. Much more generality
/// than actually needed is carried over from older ChannelMaps.
///
/// Any gdml before v3 should stay configured to use ChannelMap35Alg, and 
/// any gdml v3 or later should be configured to use ChannelMap35OptAlg.
/// This is done in LBNEGeometryHelper using the fcl parameter DetectorVersion
/// in the SortingParameters pset.
///
///
#ifndef GEO_CHANNEL35OPTMAPALG_H
#define GEO_CHANNEL35OPTMAPALG_H

#include <vector>
#include <set>
#include <stdint.h>

#include "Geometry/ChannelMapAlg.h"
#include "lbne/Geometry/GeoObjectSorter35.h"
#include "fhiclcpp/ParameterSet.h"

namespace geo{

  class ChannelMap35OptAlg : public ChannelMapAlg{

  public:

    ChannelMap35OptAlg(fhicl::ParameterSet const& p);
    ~ChannelMap35OptAlg();
    
    void                     Initialize( std::vector<geo::CryostatGeo*> & cgeo,
					 std::vector<geo::AuxDetGeo*>   & adgeo );
    void                     Uninitialize();
    std::vector<WireID>      ChannelToWire(uint32_t channel)        const;
    uint32_t                 Nchannels()                            const;
    virtual double WireCoordinate(double YPos, double ZPos,
                                  unsigned int PlaneNo,
                                  unsigned int TPCNo,
                                  unsigned int cstat) const override;
    WireID                   NearestWireID(const TVector3& worldPos,
					   unsigned int    PlaneNo,
					   unsigned int    TPCNo,
					   unsigned int    cstat)   const;
    uint32_t                 PlaneWireToChannel(unsigned int plane,
						unsigned int wire,
						unsigned int tpc,
						unsigned int cstat) const;
    View_t                   View( uint32_t const channel )         const;
    SigType_t                SignalType( uint32_t const channel)    const;
    std::set<View_t>  const& Views()                                const;
    std::set<PlaneID> const& PlaneIDs()                             const;

    unsigned int NOpChannels(int NOpDets) const;

    unsigned int OpChannel(int detNum, int channel = 0) const;
    unsigned int OpDetFromOpChannel(int opChannel) const;
    unsigned int HardwareChannelFromOpChannel(int opChannel) const;
    
  private:
    
    unsigned int                                         fNcryostat;      ///< number of cryostats in the detector
    uint32_t                                             fNchannels;      ///< number of channels in the detector
    uint32_t                                             fTopChannel;     ///< book keeping highest channel #
    std::vector<unsigned int>                            fNTPC;           ///< number of TPCs in each cryostat
    std::set<View_t>                                     fViews;          ///< vector of the views present in the detector
    std::set<PlaneID>                                    fPlaneIDs;       ///< vector of the PlaneIDs present in the detector

    std::vector< unsigned int >				 fWiresInPlane;
    unsigned int					 fPlanesPerAPA;   
    uint32_t					         fChannelsPerAPA;
    std::vector<std::vector<std::vector<unsigned int>>>	 nAnchoredWires;

    std::vector<std::vector<std::vector<unsigned int>>>  fWiresPerPlane;  ///< The number of wires in this plane 
                                                                          ///< in the heirachy
    geo::GeoObjectSorter35                               fSorter;         ///< sorts geo::XXXGeo objects
    
    /// all data we need for each APA
    typedef struct {
      double fFirstWireCenterY;
      double fFirstWireCenterZ;
      /// +1 if the wire ID order follow z (larger z, or smaller intercept => larger wire ID); -1 otherwise
      float fWireSortingInZ;
    } PlaneData_t;

    ///< collects all data we need for each plane (indices: c t p)
    std::vector<std::vector<std::vector<PlaneData_t>>> fPlaneData;
    
    std::vector< double > fWirePitch;
    std::vector< double > fOrientation;
    std::vector< double > fSinOrientation; // to explore improving speed
    std::vector< double > fCosOrientation; // to explore improving speed



  };

}
#endif // GEO_CHANNELMAP35OPTALG_H

