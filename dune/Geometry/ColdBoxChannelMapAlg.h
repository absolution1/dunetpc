/*!
 * \file    ColdBoxChannelMapAlg.h
 * \brief   Channel mapping algorithms for VD ColdBox CRP.
 * \details Adapted from ICARUSChannelMapAlg by G. Petrillo
 * \date    October 25, 2021
 * \warning Only one CRP is currently supported with 48 deg ind1
 * 
 *  vgalymov@ipnl.in2p3.fr
 */

#ifndef DUNE_COLDBOXCHANNELMAPALG_H
#define DUNE_COLDBOXCHANNELMAPALG_H

// LArSoft libraries
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryData.h"
#include "larcorealg/Geometry/GeometryDataContainers.h"
#include "larcorealg/Geometry/ReadoutDataContainers.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// framework libraries
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"

// dune specific
#include "dune/Geometry/GeoObjectSorterCRU.h"

// C/C++ standard libraries
#include <vector>
#include <cassert>
#include <utility>


// -----------------------------------------------------------------------------
// forward declarations
namespace geo {
  
  class ColdBoxChannelMapAlg;
  
} // namespace dune

// -----------------------------------------------------------------------------
namespace geo { namespace dune{ namespace vd {

  class ChannelToWireMap
  {
    
  public:

    struct ChannelsInROPStruct
    {
      raw::ChannelID_t firstChannel = raw::InvalidChannelID;
      unsigned int nChannels        = 0U;
      readout::ROPID ropid; ///< ID of the ROP we cover.
      ChannelsInROPStruct() = default;
      ChannelsInROPStruct(
			  raw::ChannelID_t firstChannel,
			  unsigned int nChannels,
			  readout::ROPID const& ropid )
	: firstChannel(firstChannel), nChannels(nChannels), ropid(ropid)
      {}
    };
    

    /// Returns data of the ROP including `channel`, `nullptr` if none.
    ChannelsInROPStruct const* find(raw::ChannelID_t channel) const
    {
      auto const dend = fROPfirstChannel.end();
      auto const iNextData = fROPfirstChannel.upper_bound( channel );
      assert(iNextData != dbegin);
      return ((iNextData == dend) && (channel >= endChannel()))
	? nullptr: &(std::prev(iNextData)->second);
    }

    
    /// Returns data of the ROP `ropid`, `nullptr` if none.
    ChannelsInROPStruct const* find(readout::ROPID const& ropid) const
    {
      auto const dbegin = fROPfirstChannel.begin(), dend = fROPfirstChannel.end();
      for( auto it = dbegin; it != dend; ++it ){
	if (it->second.ropid == ropid) return &(it->second);
      }
      return nullptr;
    }
    
    /// Returns the ID of the first invalid channel (the last channel, plus 1).
    raw::ChannelID_t endChannel() const { return fEndChannel; }
    
    /// Returns the number of mapped channels.
    unsigned int nChannels() const
    { return endChannel(); }
    
    /// Sets the ID of the channels after the last valid one.
    void setEndChannel(raw::ChannelID_t channel) { fEndChannel = channel; }
    
    /// Resets the data of the map to like just constructed.
    void clear(){
      fROPfirstChannel.clear();
      fEndChannel = raw::ChannelID_t{ 0 };
    }

    void addROP( readout::ROPID const& rid,
		 raw::ChannelID_t firstROPchannel, unsigned int nChannels )
    {
      assert( fROPfirstChannel.find( firstROPchannel ) != fROPfirstChannel.end() );
      fROPfirstChannel[ firstROPchannel ] = {firstROPchannel, nChannels, rid};
    }

  private:
				       
    std::map<raw::ChannelID_t, ChannelsInROPStruct> fROPfirstChannel;
    raw::ChannelID_t fEndChannel = 0;
  };
    }}} //namespaces

// -----------------------------------------------------------------------------

class geo::ColdBoxChannelMapAlg: public geo::ChannelMapAlg {
  
  // import definitions
  using TPCColl_t   = std::vector<geo::TPCID>;
  using PlaneColl_t = std::vector<geo::PlaneID>;
  
 public:
  
  ColdBoxChannelMapAlg(fhicl::ParameterSet const& p);
  
  /// Prepares the algorithm extracting information from the geometry.
  virtual void Initialize(geo::GeometryData_t const& geodata) override;
  
  
  /// Frees the resources of this algorithm.
  virtual void Uninitialize() override;
  
  
  // --- BEGIN -- Channel mapping ----------------------------------------------
  /**
   * @brief Returns a collection of ID of wires connected to the `channel`.
   * @param channel TPC readout channel number
   * @return collection of the wire IDs associated with `channel`
   * @throws cet::exception (category: "Geometry") if non-existent channel
   * 
   * If the TPC readout `channel` is invalid or non-existing, an exception is
   * thrown.
   * Channels are expected to be associated with at least one wire.
   */
  virtual std::vector<geo::WireID> ChannelToWire(raw::ChannelID_t channel) const
    override;
  
  /// Returns the number of readout channels (ID's go `0` to `Nchannels()`).
  virtual unsigned int Nchannels() const override;
  
  /// @brief Returns the number of channels in the specified ROP
  /// @return number of channels in the specified ROP, 0 if non-existent
  virtual unsigned int Nchannels(readout::ROPID const& ropid) const override;
  
  //@{
  virtual raw::ChannelID_t PlaneWireToChannel
    (geo::WireID const& wireID) const override;
  virtual raw::ChannelID_t PlaneWireToChannel(unsigned int plane,
                                              unsigned int wire,
                                              unsigned int tpc,
                                              unsigned int cstat) const override
    { return PlaneWireToChannel(geo::WireID(cstat, tpc, plane, wire)); }
  //@}

  /// @}
  // --- END -- Channel mapping ------------------------------------------------
  
  
  /**
   * @name Deprecated functions.
   *
   * These methods are legacy and might be replaced by `geo::GeometryCore`
   * calls.
   */
  /// @{
  //@{
  virtual double WireCoordinate
    (double YPos, double ZPos, geo::PlaneID const& planeID) const override;
  virtual double WireCoordinate(double YPos, double ZPos,
                               unsigned int PlaneNo,
                               unsigned int TPCNo,
                               unsigned int cstat) const override
    { return WireCoordinate(YPos, ZPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
  //@}
  
  //@{
  virtual geo::WireID NearestWireID
    (const TVector3& worldPos, geo::PlaneID const& planeID) const override;
  virtual geo::WireID NearestWireID(const TVector3& worldPos,
                               unsigned int    PlaneNo,
                               unsigned int    TPCNo,
                               unsigned int    cstat) const override
    { return NearestWireID(worldPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
  //@}
  
  virtual std::set<geo::PlaneID> const& PlaneIDs() const override;
  
  /// @}
  
  
  //
  // TPC set interface
  //
  /// @name TPC set mapping
  /// @{
  /**
   * @brief Returns the total number of TPC sets in the specified cryostat.
   * @param cryoid cryostat ID
   * @return number of TPC sets in the cryostat, or 0 if no cryostat found
   */
  virtual unsigned int NTPCsets
    (readout::CryostatID const& cryoid) const override;

  /// Returns the largest number of TPC sets any cryostat in the detector has.
  virtual unsigned int MaxTPCsets() const override;

  /// Returns whether we have the specified TPC set.
  /// @return whether the TPC set is valid and exists
  virtual bool HasTPCset(readout::TPCsetID const& tpcsetid) const override;

  /**
   * @brief Returns the ID of the TPC set the specified TPC belongs to
   * @param tpcid ID of the TPC
   * @return the ID of the corresponding TPC set, or invalid ID when tpcid is
   *
   * If the TPC ID is not valid, an invalid TPC set ID is returned.
   * Note that this check is performed on the validity of the TPC ID, that
   * does not necessarily imply that the TPC specified by the ID actually
   * exists.
   */
  virtual readout::TPCsetID TPCtoTPCset
    (geo::TPCID const& tpcid) const override;

  /**
   * @brief Returns a list of ID of TPCs belonging to the specified TPC set
   * @param tpcsetid ID of the TPC set to convert into TPC IDs
   * @return the list of TPCs, empty if TPC set is invalid
   *
   * Each TPC set contains one TPC if on first induction plane, where wires are
   * split, or two in the second induction and collection planes, which have
   * shared channels.
   */
  virtual std::vector<geo::TPCID> TPCsetToTPCs
    (readout::TPCsetID const& tpcsetid) const override;

  /// Returns the ID of the first TPC belonging to the specified TPC set
  virtual geo::TPCID FirstTPCinTPCset
    (readout::TPCsetID const& tpcsetid) const override;

  /// @} TPC set mapping


  
  // --- BEGIN -- Readout plane interface --------------------------------------
  /// @name Readout plane mapping
  /// @{

  /**
   * @brief Returns the total number of readout planes in the specified TPC set.
   * @param tpcsetid TPC set ID
   * @return number of readout planes in the TPC sets, or `0` if ID is invalid
   *
   * The validity of the requested TPC set is performed like in `HasTPCset()`.
   */
  virtual unsigned int NROPs
    (readout::TPCsetID const& tpcsetid) const override;

  /// Returns the largest number of ROPs a TPC set in the detector has.
  virtual unsigned int MaxROPs() const override;

  /// Returns whether we have the specified ROP
  /// @return whether the readout plane is valid and exists
  virtual bool HasROP(readout::ROPID const& ropid) const override;

  /**
   * @brief Returns the ID of the ROP planeid belongs to, or invalid if none.
   * @param planeid ID of the plane
   * @return the ID of the corresponding ROP, or invalid ID when `planeid` is
   *
   * If the plane ID is not valid, an invalid readout plane ID is returned.
   * Note that this check is performed on the validity of the plane ID, that
   * does not necessarily imply that the plane specified by the ID actually
   * exists.
   */
  virtual readout::ROPID WirePlaneToROP
    (geo::PlaneID const& planeid) const override;

  /**
   * @brief Returns a list of ID of wire planes belonging to the specified ROP.
   * @param ropid ID of the readout plane to convert into wire planes
   * @return the list of wire plane IDs, empty if readout plane ID is invalid
   *
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists.
   * 
   * In this mapping, readout planes contain one or two wire planes each,
   * depending on the view.
   */
  virtual std::vector<geo::PlaneID> ROPtoWirePlanes
    (readout::ROPID const& ropid) const override;

  /**
   * @brief Returns a list of ID of TPCs the specified ROP spans
   * @param ropid ID of the readout plane
   * @return the list of TPC IDs, empty if readout plane ID is invalid
   *
   * In this mapping, readout planes and wire planes are mapped one-to-one.
   * The returned list contains always one entry, unless the specified readout
   * plane ID is invalid, in which case the list is empty.
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists. Check if the ROP exists with HasROP().
   * The behaviour on non-existing readout planes is undefined.
   */
  virtual std::vector<geo::TPCID> ROPtoTPCs
    (readout::ROPID const& ropid) const override;

  /// Returns the ID of the ROP the channel belongs to (invalid if none).
  virtual readout::ROPID ChannelToROP
    (raw::ChannelID_t channel) const override;

  /**
   * @brief Returns the ID of the first channel in the specified readout plane.
   * @param ropid ID of the readout plane
   * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
   *
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists. Check if the ROP exists with HasROP().
   * The behaviour for non-existing readout planes is undefined.
   */
  virtual raw::ChannelID_t FirstChannelInROP
    (readout::ROPID const& ropid) const override;

  /**
   * @brief Returns the ID of the first plane belonging to the specified ROP.
   * 
   * The wire planes within a readout plane are supposed to be sorted by beam
   * (_z_) coordinate, so that the first one should be the closest to the
   * beam entrance point.
   */
  virtual geo::PlaneID FirstWirePlaneInROP
    (readout::ROPID const& ropid) const override;

  /// @}
  // --- END -- Readout plane interface ----------------------------------------

  
  /// Return the sorter.
  virtual geo::GeoObjectSorter const& Sorter() const override
    { return fSorter; }

    private:
  
  /// Type for counts of wireless channels: per TPC set (even/odd), then per
  /// plane starting from the closest to the cathode; `first` is the number of
  /// wireless channels before the regular ones, `second` is the one after them.
  //using WirelessChannelCounts_t
  //= std::array<std::array<std::pair<unsigned int, unsigned int>, 3U>, 2U>;
  

  
  // --- BEGIN -- Readout element information ----------------------------------
  /**
   * @name Readout element information
   * 
   * The geometry and readout data containers have currently no support for
   * resizing and their size is assigned on construction.
   * 
   * Access should happen via the corresponding member functions.
   * 
   */
  /// @{
  
  struct ChannelRange_t
    : std::pair<raw::ChannelID_t, raw::ChannelID_t>
    {
      
      using std::pair<raw::ChannelID_t, raw::ChannelID_t>::pair;

      /// Returns the ID of the first channel in the range.
      constexpr raw::ChannelID_t begin() const { return first; }

      /// Returns the ID of the channel after the last one in the range.
      constexpr raw::ChannelID_t end() const { return second; }
      
      /// Returns whether this range includes the specified `channel`.
      constexpr bool contains(raw::ChannelID_t channel) const
      { return (channel >= begin()) && (channel < end()); }
    };

  
  /// Collected information about TPC sets and readout planes in the geometry.
  struct ReadoutMappingInfo_t {
    /// Number of TPC sets in each cryostat.
    std::vector<unsigned int> fTPCsetCount;
    
    /// All `geo::TPCGeo` objects in each TPC set.
    readout::TPCsetDataContainer<TPCColl_t> fTPCsetTPCs;
    
    /// Number of readout planes in each TPC set.
    readout::TPCsetDataContainer<unsigned int> fROPcount;
    
    /// All `geo::PlaneGeo` objects in each readout plane.
    readout::ROPDataContainer<PlaneColl_t> fROPplanes;
    
    /// The TPC set each TPC belongs to.
    geo::TPCDataContainer<readout::TPCsetID> fTPCtoTPCset;
    
    /// The ROP each wire plane belongs to.
    geo::PlaneDataContainer<readout::ROPID> fPlaneToROP;
    
    ReadoutMappingInfo_t() = default;
    
    void set(
      std::vector<unsigned int>&&                  TPCsetCount,
      readout::TPCsetDataContainer<TPCColl_t>&&    TPCsetTPCs,
      readout::TPCsetDataContainer<unsigned int>&& ROPcount,
      readout::ROPDataContainer<PlaneColl_t>&&     ROPplanes,
      geo::TPCDataContainer<readout::TPCsetID>&&   TPCtoTPCset,
      geo::PlaneDataContainer<readout::ROPID>&&    PlaneToROP
      )
      {
        fTPCsetCount = std::move(TPCsetCount);
        fTPCsetTPCs  = std::move(TPCsetTPCs );
        fROPcount    = std::move(ROPcount   );
        fROPplanes   = std::move(ROPplanes  );
        fTPCtoTPCset = std::move(TPCtoTPCset);
        fPlaneToROP  = std::move(PlaneToROP );
        assert(fTPCsetCount.size() == fTPCsetTPCs.dimSize<0U>());
        assert(fTPCsetCount.size() == fROPcount.dimSize<0U>());
        assert(fTPCsetCount.size() == fROPplanes.dimSize<0U>());
        assert(fTPCsetCount.size() == fTPCtoTPCset.dimSize<0U>());
        assert(fTPCsetCount.size() == fPlaneToROP.dimSize<0U>());
        assert(fTPCsetTPCs.dimSize<1U>() == fROPcount.dimSize<1U>());
        assert(fTPCsetTPCs.dimSize<1U>() == fROPplanes.dimSize<1U>());
        assert(fTPCtoTPCset.dimSize<1U>() == fPlaneToROP.dimSize<1U>());
      } // set()
    
    unsigned int NCryostats() const
      { return fROPplanes.dimSize<0U>(); }
    unsigned int MaxTPCsets() const { return fROPplanes.dimSize<1U>(); }
    unsigned int MaxROPs() const { return fROPplanes.dimSize<2U>(); }
    
    /// Frees the memory and leaves the object unusable until next `set()`.
    void clear()
      {
        fTPCsetCount.clear(); fTPCsetTPCs.clear();
        fROPcount.clear();    fROPplanes.clear();
        fTPCtoTPCset.clear(); fPlaneToROP.clear();
      }
    
    /// Returns whether all the data containers are initialized.
    operator bool() const
      {
        return !fTPCsetCount.empty() && !fTPCsetTPCs.empty()
          && !fROPcount.empty() && !fROPplanes.empty()
          && !fTPCtoTPCset.empty() && !fPlaneToROP.empty();
      }
    
  }; // ReadoutMappingInfo_t
  
  /// Collection of information on one plane.
  struct PlaneInfo_t {
    
    ChannelRange_t fChannelRange; ///< Range of channels covered by the plane.
    readout::ROPID fROPID;        ///< Which readout plane this wire plane belongs to.
    
    /// Returns the range of channels covered by the wire plane.
    constexpr ChannelRange_t const& channelRange() const
      { return fChannelRange; }
    
    /// Returns the ID of the last channel in the range.
    constexpr raw::ChannelID_t firstChannel() const
      { return fChannelRange.begin(); }
    
    /// Returns the ID of the last channel in the range.
    constexpr raw::ChannelID_t lastChannel() const
      { return fChannelRange.end() - 1; }
    
    /// Returns the ID of the channel after the last in the range.
    constexpr raw::ChannelID_t endChannel() const
      { return fChannelRange.end(); }
    
    /// Returns the ID of the readout plane this wire plane belongs to.
    constexpr readout::ROPID ROP() const { return fROPID; }
    
  }; // struct PlaneInfo_t
  
  
  /// Information about TPC sets and readout planes in the geometry.
  ReadoutMappingInfo_t fReadoutMapInfo;
  
  /// Mapping of channels and ROP's.
  geo::dune::vd::ChannelToWireMap fChannelToWireMap;
  
  /// Range of channels covered by each of the wire planes.
  geo::PlaneDataContainer<PlaneInfo_t> fPlaneInfo;
  
   geo::GeoObjectSorterCRU  fSorter;              ///< class to sort geo objects
  
  
  using PlaneType_t = std::size_t; ///< Type for plane type identifier.
  
  /// Identifier for first induction plane type.
  static constexpr PlaneType_t kFirstInductionType = 0U;
  /// Identifier for second induction plane type.
  static constexpr PlaneType_t kSecondInductionType = 1U;
  /// Identifier for collection plane type.
  static constexpr PlaneType_t kCollectionType = 2U;
  /// Identifier for unknown plane type.
  static constexpr PlaneType_t kUnknownType
    = std::numeric_limits<PlaneType_t>::max();
  
  
  // --- BEGIN -- Readout element information access ---------------------------
  /// @name Readout element information access
  /// @{
  
  /// Returns the number of TPC sets in each cryostat.
  std::vector<unsigned int> const& TPCsetCount() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fTPCsetCount; }
  
  /// Returns the number of TPC sets in the specified cryostat `cid`.
  unsigned int TPCsetCount(readout::CryostatID const& cid) const
    { return TPCsetCount()[cid.Cryostat]; }
  
  /// All `geo::TPCGeo` objects in each TPC set
  readout::TPCsetDataContainer<TPCColl_t> const& TPCsetTPCs() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fTPCsetTPCs; }
  
  /// All `geo::TPCID` objects in the specified TPC set `sid`.
  TPCColl_t const& TPCsetTPCs(readout::TPCsetID const& sid) const
    { return TPCsetTPCs()[sid]; }
  
  /// Number of readout planes in each TPC set.
  readout::TPCsetDataContainer<unsigned int> const& ROPcount() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fROPcount; }
  
  /// Number of readout planes in the specified TPC set `sid`.
  unsigned int ROPcount(readout::TPCsetID const& sid) const
    { return ROPcount()[sid]; }
  
  /// All `geo::PlaneGeo` objects in each readout plane
  readout::ROPDataContainer<PlaneColl_t> const& ROPplanes() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fROPplanes; }
  
  /// All `geo::PlaneId` objects in the specified readout plane `rid`.
  PlaneColl_t const& ROPplanes(readout::ROPID const& rid) const
    { return ROPplanes()[rid]; }
  
  /// The TPC set including each TPC.
  geo::TPCDataContainer<readout::TPCsetID> const& TPCtoTPCset() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fTPCtoTPCset; }
  
  /// The readout plane including each wire plane.
  geo::PlaneDataContainer<readout::ROPID> const& PlaneToROP() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fPlaneToROP; }
  
  /// The readout plane the specified wire plane `pid` belongs to.
  readout::ROPID const& PlaneToROP(geo::PlaneID const& pid) const
    { return PlaneToROP()[pid]; }
    
  /// @}
  // --- END -- Readout element information access -----------------------------
  
  
  /// Returns whether the specified cryostat is known to the mapping.
  bool HasCryostat(readout::CryostatID const& cryoid) const;
  
  
  /**
   * @brief Fills the information about readout channel mapping.
   * @param Cryostats the sorted list of cryostats in the detector
   * 
   * 
   * 
   * The readout information must have been already filled
   * (`buildReadoutPlanes()`).
   * 
   */
  void fillChannelToWireMap
    (geo::GeometryData_t::CryostatList_t const& Cryostats);
  
  
  /**
   * @brief Fills information about the TPC set and readout plane structure.
   * @param Cryostats the sorted list of cryostats in the detector
   * 
   */
  void buildReadoutPlanes(geo::GeometryData_t::CryostatList_t const& Cryostats);
  
  
  /**
   * @brief Returns the "type" of readout plane.
   * @param ropid ID of the readout plane to query
   * 
   * Here the "type" refers to the position of the plane in the TPC and as a
   * consequence to its signal type.
   * The type of the readout plane is deduced from the type of the planes it
   * contains.
   * 
   * The returned value is `kFirstInductionType`, `kSecondInductionType`,
   * `kCollectionType` or `kUnknownType`.
   */
  PlaneType_t findPlaneType(readout::ROPID const& ropid) const;


  /// Returns the type of signal on the specified `channel`.
  virtual geo::SigType_t SignalTypeForChannelImpl
    (raw::ChannelID_t const channel) const override;
  
  /// Returns the name of the specified plane type.
  static std::string PlaneTypeName(PlaneType_t planeType);
  
  //
  std::string fLogCategory = "ColdBoxChannelMapAlg";
  
}; // class geo::ColdBoxChannelMapAlg


#endif //

