////////////////////////////////////////////////////////////////////////
/// \file  ChannelMap35Alg.h
/// \brief The class of 35t specific algorithms
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////
///
/// Any gdml before v3 should stay configured to use ChannelMap35Alg, and 
/// any gdml v3 or later should be configured to use ChannelMap35OptAlg.
/// This is done in DUNEGeometryHelper using the fcl parameter DetectorVersion
/// in the SortingParameters pset.
///
#ifndef GEO_CHANNEL35MAPALG_H
#define GEO_CHANNEL35MAPALG_H

#include <vector>
#include <set>

#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/ChannelMapAlg.h"
#include "dune/Geometry/GeoObjectSorter35.h"
#include "fhiclcpp/ParameterSet.h"

namespace geo{

  class ChannelMap35Alg : public ChannelMapAlg{

  public:

    ChannelMap35Alg(fhicl::ParameterSet const& p);
    
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
    
    unsigned int                                         fNcryostat;      ///< number of cryostats in the detector
    unsigned int                                         fNchannels;      ///< number of channels in the detector
    raw::ChannelID_t                                     fTopChannel;     ///< book keeping highest channel #
    std::vector<unsigned int>                            fNTPC;           ///< number of TPCs in each cryostat
    std::set<View_t>                                     fViews;          ///< vector of the views present in the detector
    std::set<PlaneID>                                    fPlaneIDs;       ///< vector of the PlaneIDs present in the detector

    unsigned int                                         fPlanesPerAPA;
    raw::ChannelID_t                                     fChannelsPerAPA;
    PlaneInfoMap_t<unsigned int>                         nAnchoredWires;

    PlaneInfoMap_t<unsigned int>                         fWiresPerPlane;  ///< The number of wires in this plane 
                                                                          ///< in the heirachy
    geo::GeoObjectSorter35                               fSorter;         ///< sorts geo::XXXGeo objects
    
    /// all data we need for each APA
    typedef struct {
      double fFirstWireCenterY;
      double fFirstWireCenterZ;
      /// +1 if the wire ID order follow z (larger z, or smaller intercept => larger wire ID); -1 otherwise
      float fWireSortingInZ;
      double fYmax;
      double fYmin;
      double fZmax;
      double fZmin;
    } PlaneData_t;
    
    ///< collects all data we need for each plane (indices: c t p)
    PlaneInfoMap_t<PlaneData_t>                          fPlaneData;
    
    std::vector< double > fWirePitch;
    std::vector< double > fOrientation;
    std::vector< double > fSinOrientation; // to explore improving speed
    std::vector< double > fCosOrientation; // to explore improving speed
    
    template <typename T>
    T const& AccessAPAelement
      (PlaneInfoMap_t<T> const& data, geo::PlaneID planeid) const
      { planeid.TPC /= 2; return AccessElement(data, planeid); }
    unsigned int WiresPerPlane(geo::PlaneID const& planeid) const
      { return AccessAPAelement(fWiresPerPlane, planeid); }
    unsigned int AnchoredWires(geo::PlaneID const& planeid) const
      { return AccessAPAelement(nAnchoredWires, planeid); }
    
    
  };

}
#endif // GEO_CHANNELMAP35ALG_H

