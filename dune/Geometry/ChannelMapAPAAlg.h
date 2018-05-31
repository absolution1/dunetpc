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

#include "cetlib_except/exception.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h" // readout::ROPID, ...
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "dune/Geometry/GeoObjectSorterAPA.h"

#include "fhiclcpp/ParameterSet.h"

namespace geo{

  class ChannelMapAPAAlg : public ChannelMapAlg{

  public:

    ChannelMapAPAAlg(fhicl::ParameterSet const& p);
    
    void                     Initialize( GeometryData_t const& geodata ) override;
    void                     Uninitialize() override;
    
    /// Returns a list of TPC wires connected to the specified readout channel ID
    /// @throws cet::exception (category: "Geometry") if non-existent channel
    std::vector<WireID>      ChannelToWire(raw::ChannelID_t channel) const override;
    
    unsigned int             Nchannels()                            const override;
    
    /// @brief Returns the number of channels in the specified ROP
    /// @return number of channels in the specified ROP, 0 if non-existent
    /// @todo to be completed
    virtual unsigned int     Nchannels(readout::ROPID const& ropid) const override;
    
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
    SigType_t                SignalTypeForChannelImpl( raw::ChannelID_t const channel) const override;
    std::set<View_t>  const& Views()                                     const;
    std::set<PlaneID> const& PlaneIDs()                                  const override;

    //
    // TPC set interface
    //
    /// @name TPC set mapping
    /// @{
    /**
     * @brief Returns the total number of TPC sets in the specified cryostat
     * @param cryoid cryostat ID
     * @return number of TPC sets in the cryostat, or 0 if no cryostat found
     *
     * @todo to be completed
     */
    virtual unsigned int NTPCsets
      (readout::CryostatID const& cryoid) const override;
    
    /// Returns the largest number of TPC sets any cryostat in the detector has
    /// @todo to be completed
    virtual unsigned int MaxTPCsets() const override;
    
    /// Returns whether we have the specified TPC set
    /// @return whether the TPC set is valid and exists
    /// @todo to be completed
    virtual bool HasTPCset(readout::TPCsetID const& tpcsetid) const override;
    
    /**
     * @brief Returns the ID of the TPC set the specified TPC belongs to
     * @param tpcid ID of the TPC
     * @return the ID of the corresponding TPC set
     *
     * Note that the check is performed on the validity of the TPC set ID, that
     * does not necessarily imply that the TPC set specified by the ID actually
     * exists. Check the existence of the TPC set first (HasTPCset()).
     * Behaviour on valid, non-existent TPC set IDs is undefined.
     * 
     * @todo to be completed
     */
    virtual readout::TPCsetID TPCtoTPCset
      (geo::TPCID const& tpcid) const override;
    
    /**
     * @brief Returns a list of ID of TPCs belonging to the specified TPC set
     * @param tpcsetid ID of the TPC set to convert into TPC IDs
     * @return the list of TPCs, empty if TPC set is invalid
     *
     * Note that the check is performed on the validity of the TPC set ID, that
     * does not necessarily imply that the TPC set specified by the ID actually
     * exists. Check the existence of the TPC set first (HasTPCset()).
     * Behaviour on valid, non-existent TPC set IDs is undefined.
     *
     * @todo to be completed
     */
    virtual std::vector<geo::TPCID> TPCsetToTPCs
      (readout::TPCsetID const& tpcsetid) const override;
    
    /// Returns the ID of the first TPC belonging to the specified TPC set
    /// @todo to be completed
    virtual geo::TPCID FirstTPCinTPCset
      (readout::TPCsetID const& tpcsetid) const override;
    
    /// @} TPC set mapping
    
    
    
    //
    // Readout plane interface
    //
    /// @name Readout plane mapping
    /// @{

    /**
     * @brief Returns the total number of ROP in the specified TPC set
     * @param tpcsetid TPC set ID
     * @return number of readout planes in the TPC set, or 0 if no TPC set found
     * 
     * Note that this methods explicitly check the existence of the TPC set.
     *
     * @todo to be completed
     */
    virtual unsigned int NROPs
      (readout::TPCsetID const& tpcsetid) const override;
    
    /// Returns the largest number of ROPs a TPC set in the detector has
    /// @todo to be completed
    virtual unsigned int MaxROPs() const override;
    
    /// Returns whether we have the specified ROP
    /// @return whether the readout plane is valid and exists
    virtual bool HasROP(readout::ROPID const& ropid) const override;
    
    /**
     * @brief Returns the ID of the ROP planeid belongs to, or invalid if none
     * @param planeid ID of the wire plane
     * @return the ID of the corresponding ROP, or invalid ID when planeid is
     *
     * Note that this check is performed on the validity of the wire plane
     * ID, that does not necessarily imply that the wire plane specified by
     * the ID actually exists. Check if the wire plane exists with HasPlane().
     * The behaviour on non-existing wire planes is undefined.
     * 
     * @todo to be completed
     */
    virtual readout::ROPID WirePlaneToROP
      (geo::PlaneID const& planeid) const override;
    
    /**
     * @brief Returns a list of ID of wire planes belonging to the specified ROP
     * @param ropid ID of the readout plane to convert into wire planes
     * @return the list of wire plane IDs, empty if readout plane ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour on non-existing readout planes is undefined.
     * 
     * @todo to be completed
     */
    virtual std::vector<geo::PlaneID> ROPtoWirePlanes
      (readout::ROPID const& ropid) const override;
    
    /**
     * @brief Returns a list of ID of TPCs the specified ROP spans
     * @param ropid ID of the readout plane
     * @return the list of TPC IDs, empty if readout plane ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour on non-existing readout planes is undefined.
     *
     * @todo to be completed
     */
    virtual std::vector<geo::TPCID> ROPtoTPCs
      (readout::ROPID const& ropid) const override;
    
    /// Returns the ID of the ROP the channel belongs to (invalid if none)
    /// @throws cet::exception (category: "Geometry") if non-existent channel
    virtual readout::ROPID ChannelToROP
      (raw::ChannelID_t channel) const override;
    
    /**
     * @brief Returns the ID of the first channel in the specified readout plane
     * @param ropid ID of the readout plane
     * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
     * 
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour for non-existing readout planes is undefined.
     *
     * @todo to be completed
     */
    virtual raw::ChannelID_t FirstChannelInROP
      (readout::ROPID const& ropid) const override;
    
    /**
     * @brief Returns the ID of the first plane belonging to the specified ROP
     * @param ropid ID of the readout plane
     * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
     * 
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour for non-existing readout planes is undefined.
     *
     * @todo to be completed
     */
    virtual geo::PlaneID FirstWirePlaneInROP
      (readout::ROPID const& ropid) const override;
    
    /// @} readout plane mapping
    
    /// Returns the object to sort geometry with
    virtual geo::GeoObjectSorter const& Sorter() const override
      { return fSorter; }
    
    
    
    unsigned int NOpChannels(unsigned int NOpDets)                        const override;
    unsigned int NOpHardwareChannels(unsigned int opDet)                  const override;
    unsigned int OpChannel(unsigned int detNum, unsigned int channel = 0) const override;
    unsigned int OpDetFromOpChannel(unsigned int opChannel)               const override;
    unsigned int HardwareChannelFromOpChannel(unsigned int opChannel)     const override;
    
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

  };

}
#endif // GEO_CHANNELMAPAPAALG_H

