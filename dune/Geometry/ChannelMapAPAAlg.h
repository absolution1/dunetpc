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
#include <set>

#include "cetlib/exception.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/ChannelMapAlg.h"
#include "dune/Geometry/GeoObjectSorterAPA.h"

#include "fhiclcpp/ParameterSet.h"

namespace geo{

  class ChannelMapAPAAlg : public ChannelMapAlg{

  public:

    ChannelMapAPAAlg(fhicl::ParameterSet const& p);
    
    void                     Initialize( GeometryData_t& geodata ) override;
    void                     Uninitialize();
    std::vector<WireID>      ChannelToWire(raw::ChannelID_t channel) const;
    unsigned int             Nchannels()                            const;
    //@{
    virtual double WireCoordinate(double YPos,
                                  double ZPos,
                                  unsigned int PlaneNo,
                                  unsigned int TPCNo,
                                  unsigned int cstat) const override
      { return WireCoordinate(YPos, ZPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    virtual double WireCoordinate(double YPos,
                                  double ZPos,
                                  geo::PlaneID const& planeID) const override;
    //@}
    
    //@{
    virtual WireID NearestWireID(const TVector3& worldPos,
                                 unsigned int    PlaneNo,
                                 unsigned int    TPCNo,
                                 unsigned int    cstat) const override
      { return NearestWireID(worldPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    virtual WireID NearestWireID
      (const TVector3& worldPos, geo::PlaneID const& planeID) const override;
    //@}
    //@{
    virtual raw::ChannelID_t PlaneWireToChannel(unsigned int plane,
                                                unsigned int wire,
                                                unsigned int tpc,
                                                unsigned int cstat) const override
      { return PlaneWireToChannel(geo::WireID(cstat, tpc, plane, wire)); }
    virtual raw::ChannelID_t PlaneWireToChannel(geo::WireID const& wireID) const override;
    //@}
    View_t                   View( raw::ChannelID_t const channel )      const;
    SigType_t                SignalType( raw::ChannelID_t const channel) const;
    std::set<View_t>  const& Views()                                     const;
    std::set<PlaneID> const& PlaneIDs()                                  const;

    unsigned int NOpChannels(unsigned int NOpDets)                        const;
    unsigned int NOpHardwareChannels(unsigned int opDet)                  const;
    unsigned int OpChannel(unsigned int detNum, unsigned int channel = 0) const;
    unsigned int OpDetFromOpChannel(unsigned int opChannel)               const;
    unsigned int HardwareChannelFromOpChannel(unsigned int opChannel)     const;
    
  private:
    
    unsigned int                                         fChannelsPerOpDet;
    unsigned int                                         fNcryostat;      ///< number of cryostats in the detector
    unsigned int                                         fNchannels;      ///< number of channels in the detector
    raw::ChannelID_t                                     fTopChannel;     ///< book keeping highest channel #
    std::vector<unsigned int>                            fNTPC;           ///< number of TPCs in each cryostat
    std::set<View_t>                                     fViews;          ///< vector of the views present in the detector
    std::set<PlaneID>                                    fPlaneIDs;       ///< vector of the PlaneIDs present in the detector
    // Assuming all APA's are identical
    std::vector< unsigned int >                                 fWiresInPlane;
    unsigned int                                         fPlanesPerAPA;   
    unsigned int                                             fChannelsPerAPA;
    std::vector< unsigned int >                                 nAnchoredWires;

    PlaneInfoMap_t<unsigned int>                         fWiresPerPlane;  ///< The number of wires in this plane 
                                                                          ///< in the heirachy

    geo::GeoObjectSorterAPA                              fSorter;         ///< sorts geo::XXXGeo objects

    /// all data we need for each APA
    typedef struct {
      double fFirstWireCenterY;
      double fFirstWireCenterZ;
      /// +1 if the wire ID order follow z (larger z, or smaller intercept => larger wire ID); -1 otherwise
      float fWireSortingInZ;
    } PlaneData_t;
    
    ///< collects all data we need for each plane (indices: c t p)
    PlaneInfoMap_t<PlaneData_t>                          fPlaneData;
    
    std::vector< double > fWirePitch;
    std::vector< double > fOrientation;
    std::vector< double > fSinOrientation; // to explore improving speed
    std::vector< double > fCosOrientation; // to explore improving speed

  };

}
#endif // GEO_CHANNELMAPAPAALG_H

