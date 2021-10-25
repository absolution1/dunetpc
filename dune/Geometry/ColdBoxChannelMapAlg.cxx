/*!
 * \file    ColdBoxChannelMapAlg.h
 * \brief   Channel mapping algorithms for VD ColdBox CRP.
 * \details Adapted from ICARUSChannelMapAlg by G. Petrillo
 * \date    October 25, 2021
 * \author  Vyacheslav Galymov (vgalymov@ipnl.in2p3.fr)
 * \warning Only one CRP is currently supported with 48 deg ind1
 */

#include "ColdBoxChannelMapAlg.h"


// LArSoft libraries
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/ReadoutDataContainers.h"
#include "larcorealg/Geometry/details/extractMaxGeometryElements.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::printBacktrace()
#include "larcorealg/CoreUtils/StdUtils.h" // util::size()
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h" // cet::exception

// C/C++ libraries
#include <vector>
#include <array>
#include <set>
#include <algorithm>  // std::transform(), std::find()
#include <utility>    // std::move()
#include <iterator>   // std::back_inserter()
#include <functional> // std::mem_fn()
#include <tuple>


// -----------------------------------------------------------------------------
dune::ColdBoxChannelMapAlg::ColdBoxChannelMapAlg(Config const& config)
  : fWirelessChannelCounts
    (extractWirelessChannelParams(config.WirelessChannels()))
  , fSorter(getOptionalParameterSet(config.Sorter))
  {}


// -----------------------------------------------------------------------------
void dune::ColdBoxChannelMapAlg::Initialize(geo::GeometryData_t const& geodata)
{
  // This is the only INFO level message we want this object to produce;
  // given the dynamic nature of the channel mapping choice,
  // it's better for the log to have some indication of chosen channel mapping.
  mf::LogInfo("ColdBoxChannelMapAlg")
    << "Initializing ColdBoxChannelMapAlg channel mapping algorithm.";
  
  buildReadoutPlanes(geodata.cryostats);
  fillChannelToWireMap(geodata.cryostats);
  
  MF_LOG_TRACE("ColdBoxChannelMapAlg")
    << "ColdBoxChannelMapAlg::Initialize() completed.";
  
} // dune::ColdBoxChannelMapAlg::Initialize()


// -----------------------------------------------------------------------------
void dune::ColdBoxChannelMapAlg::Uninitialize() {
  
  fReadoutMapInfo.clear();
  fChannelToWireMap.clear();
  fPlaneInfo.clear();
  
} // dune::ColdBoxChannelMapAlg::Uninitialize()


//------------------------------------------------------------------------------
std::vector<geo::WireID> dune::ColdBoxChannelMapAlg::ChannelToWire
  (raw::ChannelID_t channel) const
{
  //
  // input check
  //
  assert(!fPlaneInfo.empty());
  
  //
  // output
  //
  std::vector<geo::WireID> AllSegments;
  
  //
  // find the ROP with that channel
  //
  dune::details::ChannelToWireMap::ChannelsInROPStruct const* channelInfo
    = fChannelToWireMap.find(channel);
  if (!channelInfo) {
    throw cet::exception("Geometry")
      << "dune::ColdBoxChannelMapAlg::ChannelToWire(" << channel
      << "): invalid channel requested (must be lower than "
      << Nchannels() << ")\n";
  }
  
  //
  // find the wire planes in that ROP
  //
  PlaneColl_t const& planes = ROPplanes(channelInfo->ropid);
  
  //
  // associate one wire for each of those wire planes to the channel
  //
  AllSegments.reserve(planes.size()); // this is sometimes (often?) too much
  for (geo::PlaneGeo const* plane: planes) {
    
    geo::PlaneID const& pid = plane->ID();
    ChannelRange_t const& channelRange = fPlaneInfo[pid].channelRange();
    
    if (!channelRange.contains(channel)) continue;
    AllSegments.emplace_back
      (pid, static_cast<geo::WireID::WireID_t>(channel - channelRange.begin()));
    
  } // for planes in ROP
  
  return AllSegments;
  
} // dune::ColdBoxChannelMapAlg::ChannelToWire()


//------------------------------------------------------------------------------
unsigned int dune::ColdBoxChannelMapAlg::Nchannels() const {
  
  return fChannelToWireMap.nChannels();
  
} // dune::ColdBoxChannelMapAlg::Nchannels()


//------------------------------------------------------------------------------
unsigned int dune::ColdBoxChannelMapAlg::Nchannels
  (readout::ROPID const& ropid) const 
{
  dune::details::ChannelToWireMap::ChannelsInROPStruct const* ROPinfo
    = fChannelToWireMap.find(ropid);
  return ROPinfo? ROPinfo->nChannels: 0U;
} // dune::ColdBoxChannelMapAlg::Nchannels(ROPID)


//------------------------------------------------------------------------------
double dune::ColdBoxChannelMapAlg::WireCoordinate
  (double YPos, double ZPos, geo::PlaneID const& planeID) const
{
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ColdBoxChannelMapAlg");
  e << "ColdBoxChannelMapAlg does not support `WireCoordinate()` call."
    "\nPlease update calling software to use geo::PlaneGeo::WireCoordinate()`:"
    "\n";
  
  lar::debug::printBacktrace(e, 4U);
  
  throw e;
} // dune::ColdBoxChannelMapAlg::WireCoordinate()


//------------------------------------------------------------------------------
geo::WireID dune::ColdBoxChannelMapAlg::NearestWireID
  (const TVector3& worldPos, geo::PlaneID const& planeID) const
{
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ColdBoxChannelMapAlg");
  e << "ColdBoxChannelMapAlg does not support `NearestWireID()` call."
    "\nPlease update calling software to use geo::PlaneGeo::NearestWireID()`:"
    "\n";
  
  lar::debug::printBacktrace(e, 3U);
  
  throw e;
} // dune::ColdBoxChannelMapAlg::NearestWireID()


//------------------------------------------------------------------------------
raw::ChannelID_t dune::ColdBoxChannelMapAlg::PlaneWireToChannel
  (geo::WireID const& wireID) const
{
  return fPlaneInfo[wireID].firstChannel() + wireID.Wire;
} // dune::ColdBoxChannelMapAlg::PlaneWireToChannel()


//------------------------------------------------------------------------------
std::set<geo::PlaneID> const& dune::ColdBoxChannelMapAlg::PlaneIDs() const {
  
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ColdBoxChannelMapAlg");
  e << "ColdBoxChannelMapAlg does not support `PlaneIDs()` call."
    "\nPlease update calling software to use geo::GeometryCore::IteratePlanes()`"
    "\n";
  
  lar::debug::printBacktrace(e, 3U);
  
  throw e;
  
} // dune::ColdBoxChannelMapAlg::PlaneIDs()


//------------------------------------------------------------------------------
unsigned int dune::ColdBoxChannelMapAlg::NTPCsets
  (readout::CryostatID const& cryoid) const
{
  return HasCryostat(cryoid)? TPCsetCount(cryoid): 0U;
} // dune::ColdBoxChannelMapAlg::NTPCsets()


//------------------------------------------------------------------------------
/// Returns the largest number of TPC sets any cryostat in the detector has
unsigned int dune::ColdBoxChannelMapAlg::MaxTPCsets() const {
  assert(fReadoutMapInfo);
  return fReadoutMapInfo.MaxTPCsets();
} // dune::ColdBoxChannelMapAlg::MaxTPCsets()


//------------------------------------------------------------------------------
/// Returns whether we have the specified TPC set
/// @return whether the TPC set is valid and exists
bool dune::ColdBoxChannelMapAlg::HasTPCset
  (readout::TPCsetID const& tpcsetid) const
{
  return
    HasCryostat(tpcsetid)? (tpcsetid.TPCset < TPCsetCount(tpcsetid)): false;
} // dune::ColdBoxChannelMapAlg::HasTPCset()


//------------------------------------------------------------------------------
readout::TPCsetID dune::ColdBoxChannelMapAlg::TPCtoTPCset
  (geo::TPCID const& tpcid) const
{
  return tpcid? TPCtoTPCset()[tpcid]: readout::TPCsetID{};
} // dune::ColdBoxChannelMapAlg::TPCtoTPCset()


//------------------------------------------------------------------------------
std::vector<geo::TPCID> dune::ColdBoxChannelMapAlg::TPCsetToTPCs
  (readout::TPCsetID const& tpcsetid) const
{
  std::vector<geo::TPCID> TPCs;
  if (!tpcsetid) return TPCs;
  
  auto const& TPClist = TPCsetTPCs(tpcsetid);
  TPCs.reserve(TPClist.size());
  std::transform(TPClist.begin(), TPClist.end(), std::back_inserter(TPCs),
    std::mem_fn(&geo::TPCGeo::ID)
    );
  return TPCs;
} // dune::ColdBoxChannelMapAlg::TPCsetToTPCs()


//------------------------------------------------------------------------------
geo::TPCID dune::ColdBoxChannelMapAlg::FirstTPCinTPCset
  (readout::TPCsetID const& tpcsetid) const
{
  if (!tpcsetid) return {};
  
  auto const& TPClist = TPCsetTPCs(tpcsetid);
  return TPClist.empty()? geo::TPCID{}: TPClist.front()->ID();
  
} // dune::ColdBoxChannelMapAlg::FirstTPCinTPCset()


//------------------------------------------------------------------------------
unsigned int dune::ColdBoxChannelMapAlg::NROPs
  (readout::TPCsetID const& tpcsetid) const
{
  return HasTPCset(tpcsetid)? ROPcount(tpcsetid): 0U;
} // dune::ColdBoxChannelMapAlg::NROPs()


//------------------------------------------------------------------------------
unsigned int dune::ColdBoxChannelMapAlg::MaxROPs() const {
  assert(fReadoutMapInfo);
  return fReadoutMapInfo.MaxROPs();
} // dune::ColdBoxChannelMapAlg::MaxROPs()

//------------------------------------------------------------------------------
bool dune::ColdBoxChannelMapAlg::HasROP(readout::ROPID const& ropid) const {
  return HasTPCset(ropid)? (ropid.ROP < ROPcount(ropid)): false;
} // dune::ColdBoxChannelMapAlg::HasROP()


//------------------------------------------------------------------------------
  /**
   * @brief Returns the ID of the ROP planeid belongs to, or invalid if none
   * @param planeid ID of the plane
   * @return the ID of the corresponding ROP, or invalid ID when planeid is
   *
   * In this mapping, readout planes and wire planes are mapped one-to-one.
   * The returned value mirrors the plane ID in the readout space.
   * If the plane ID is not valid, an invalid readout plane ID is returned.
   * Note that this check is performed on the validity of the plane ID, that
   * does not necessarily imply that the plane specified by the ID actually
   * exists.
   */
readout::ROPID dune::ColdBoxChannelMapAlg::WirePlaneToROP
  (geo::PlaneID const& planeid) const
{
  return planeid? PlaneToROP(planeid): readout::ROPID{};
} // dune::ColdBoxChannelMapAlg::WirePlaneToROP()


//------------------------------------------------------------------------------
std::vector<geo::PlaneID> dune::ColdBoxChannelMapAlg::ROPtoWirePlanes
  (readout::ROPID const& ropid) const
{
  std::vector<geo::PlaneID> Planes;
  if (!ropid) return Planes;
  
  auto const& PlaneList = ROPplanes(ropid);
  Planes.reserve(PlaneList.size());
  std::transform(PlaneList.begin(), PlaneList.end(), std::back_inserter(Planes),
    std::mem_fn(&geo::PlaneGeo::ID)
    );
  return Planes;
} // dune::ColdBoxChannelMapAlg::ROPtoWirePlanes()


//------------------------------------------------------------------------------
std::vector<geo::TPCID> dune::ColdBoxChannelMapAlg::ROPtoTPCs
  (readout::ROPID const& ropid) const
{
  std::vector<geo::TPCID> TPCs;
  if (!ropid) return TPCs;
  
  auto const& PlaneList = ROPplanes(ropid);
  TPCs.reserve(PlaneList.size());
  std::transform(PlaneList.begin(), PlaneList.end(), std::back_inserter(TPCs),
    std::mem_fn(&geo::PlaneGeo::ID)
    );
  return TPCs;
} // dune::ColdBoxChannelMapAlg::ROPtoTPCs()


//------------------------------------------------------------------------------
readout::ROPID dune::ColdBoxChannelMapAlg::ChannelToROP
  (raw::ChannelID_t channel) const
{
  if (!raw::isValidChannelID(channel)) return {};
  dune::details::ChannelToWireMap::ChannelsInROPStruct const* info
    = fChannelToWireMap.find(channel);
  return info? info->ropid: readout::ROPID{};
} // dune::ColdBoxChannelMapAlg::ChannelToROP()


//------------------------------------------------------------------------------
raw::ChannelID_t dune::ColdBoxChannelMapAlg::FirstChannelInROP
  (readout::ROPID const& ropid) const
{
  if (!ropid) return raw::InvalidChannelID;
  dune::details::ChannelToWireMap::ChannelsInROPStruct const* info
    = fChannelToWireMap.find(ropid);
  return info? info->firstChannel: raw::InvalidChannelID;
} // dune::ColdBoxChannelMapAlg::FirstChannelInROP()


//------------------------------------------------------------------------------
geo::PlaneID dune::ColdBoxChannelMapAlg::FirstWirePlaneInROP
  (readout::ROPID const& ropid) const
{
  if (!ropid) return {};
  PlaneColl_t const& planes = ROPplanes(ropid);
  return planes.empty()? geo::PlaneID{}: planes.front()->ID();
} // dune::ColdBoxChannelMapAlg::FirstWirePlaneInROP()


//------------------------------------------------------------------------------
bool dune::ColdBoxChannelMapAlg::HasCryostat
  (readout::CryostatID const& cryoid) const
{
  assert(fReadoutMapInfo);
  return cryoid.Cryostat < fReadoutMapInfo.NCryostats();
} // dune::ColdBoxChannelMapAlg::HasCryostat()


//------------------------------------------------------------------------------
void dune::ColdBoxChannelMapAlg::fillChannelToWireMap
  (geo::GeometryData_t::CryostatList_t const& Cryostats)
{
  
  //
  // input check
  //
  assert(fReadoutMapInfo);
  assert(!Cryostats.empty());
  
  //
  // output setup
  //
  assert(fPlaneInfo.empty());
  std::array<unsigned int, 3U> maxSizes
    = geo::details::extractMaxGeometryElements<3U>(Cryostats);

  fPlaneInfo.resize(maxSizes[0U], maxSizes[1U], maxSizes[2U]);
  
  
  raw::ChannelID_t nextChannel = 0; // next available channel
  
  // once again we do not have iteration facilities from `geo::GeometryCore`
  // available yet, so we go the nested loop way and bite it
  for (geo::CryostatGeo const& cryo: Cryostats) {
    
    readout::CryostatID const cid { cryo.ID() };
    
    auto const nTPCsets 
      = static_cast<readout::TPCsetID::TPCsetID_t>(TPCsetCount(cid));
    
    for (readout::TPCsetID::TPCsetID_t s: util::counter(nTPCsets)) {
      
      readout::TPCsetID const sid { cid, s };
      
      // select the channel count according to whether the TPC set is even or
      // odd; the selected structure is an array with one element for wire
      // plane signal type (first induction, second induction and collection):
      auto const& TPCsetChannelCounts
        = fWirelessChannelCounts.at(sid.TPCset & 1);
      
      auto const nROPs = static_cast<readout::ROPID::ROPID_t>(ROPcount(sid));
      
      for (readout::ROPID::ROPID_t r: util::counter(nROPs)) {
        
        mf::LogTrace log("ColdBoxChannelMapAlg");
        
        readout::ROPID const rid { sid, r };
        auto const planeType = findPlaneType(rid);
        log << "ROP: " << rid
          << " (plane type: " << PlaneTypeName(planeType) << ")";
        
        auto const& WirelessChannelCounts
          = TPCsetChannelCounts.at(planeType);
        
        PlaneColl_t const& planes = ROPplanes(rid);
        log << " (" << planes.size() << " planes):";
        assert(!planes.empty());
        
        raw::ChannelID_t const firstROPchannel = nextChannel;
        
        auto iPlane = util::begin(planes);
        auto const pend = util::end(planes);
        
        // assign available channels to all wires of the first plane
        nextChannel += WirelessChannelCounts.first + (*iPlane)->Nwires();
        fPlaneInfo[(*iPlane)->ID()] = {
          ChannelRange_t
            { firstROPchannel + WirelessChannelCounts.first, nextChannel },
          rid
          };
        log << " [" << (*iPlane)->ID() << "] "
          << fPlaneInfo[(*iPlane)->ID()].firstChannel()
          << " -- " << fPlaneInfo[(*iPlane)->ID()].lastChannel() << ";";
        
        geo::Point_t lastWirePos = (*iPlane)->LastWire().GetCenter<geo::Point_t>();
        
        while (++iPlane != pend) {
          
          geo::PlaneGeo const& plane = **iPlane;
          
          // find out which wire matches the last wire from the previous plane;
          // if there is no such wire, an exception will be thrown,
          // which is ok to us since we believe it should not happen.
          geo::WireID const lastMatchedWireID
            = plane.NearestWireID(lastWirePos);
          
          /*
          mf::LogTrace("ColdBoxChannelMapAlg")
            << (*std::prev(iPlane))->ID() << " W:" << ((*std::prev(iPlane))->Nwires() - 1)
            << " ending at " << (*std::prev(iPlane))->LastWire().GetEnd<geo::Point_t>()
            << " matched " << lastMatchedWireID
            << " which starts at " << plane.Wire(lastMatchedWireID).GetStart<geo::Point_t>()
            ;
          */
          
          //
          // the first channel in this plane (`iPlane`) is the one associated
          // to the first wire in the plane, which has local wire number `0`;
          // the last channel from the previous plane (`nextChannel - 1`)
          // is associated to the matched wire (`lastMatchedWireID`),
          // which has some wire number (`lastMatchedWireID.Wire`).
          //
          auto const nWires = plane.Nwires();
          raw::ChannelID_t const firstChannel
            = (nextChannel - 1) - lastMatchedWireID.Wire;
          nextChannel = firstChannel + nWires;
          
          fPlaneInfo[plane.ID()] = { { firstChannel, nextChannel }, rid };
          log << " [" << plane.ID() << "] "
            << fPlaneInfo[plane.ID()].firstChannel() << " -- "
            << fPlaneInfo[plane.ID()].lastChannel() << ";";
          
          // update for the next iteration
          lastWirePos = plane.LastWire().GetCenter<geo::Point_t>();
          
        } // while
        
        nextChannel += WirelessChannelCounts.second;
        unsigned int const nChannels = nextChannel - firstROPchannel;
        fChannelToWireMap.addROP(rid, firstROPchannel, nChannels);
        log
          << " => " << nChannels << " channels starting at " << firstROPchannel;
        
      } // for readout plane
      
    } // for TPC set
    
  } // for cryostat
  
  fChannelToWireMap.setEndChannel(nextChannel);
  mf::LogTrace("ColdBoxChannelMapAlg")
    << "Counted " << fChannelToWireMap.nChannels() << " channels.";
  
} // dune::ColdBoxChannelMapAlg::fillChannelToWireMap()


// -----------------------------------------------------------------------------
void dune::ColdBoxChannelMapAlg::buildReadoutPlanes
  (geo::GeometryData_t::CryostatList_t const& Cryostats)
{

  if( Cryostats.size() > 1 ){
    throw cet::exception("Geometry")
      << "dune::ColdBoxChannelMapAlg::buildReadoutPlanes " << Cryostats.size()
      << ": more than one cryostat is currently not supported\n";
  }
  
  if( Cryostats[0].NTPC() != 4 ){
    throw cet::exception("Geometry")
      << "dune::ColdBoxChannelMapAlg::buildReadoutPlanes " << Cryostats[0].NTPC() != 4
      << ": more than four TPCs is currently not supported\n";
  }
  
  // currently do it by hand for CB --> to be generalized for FD
  auto cryo = Cryostats[0];
  unsigned int MaxTPCs = cryo.NTPC();
  unsigned int MaxROPs = 3;
  unsigned int MaxPlanes = 3;
  std::vector<unsigned int> TPCsetCount = {2};
  std::vector<TPCColl_t> AllTPCsInTPCsets;       // TPCs belonging to a set
  AllTPCsInTPCsets.resize( TPCsetCount.size() ); // set index
  // the standard CRU sorter is assumed (fixing this another "todo")
  AllTPCsInTPCsets[0] = {cryo.TPCPtr(0), cryo.TPCPtr(2)};
  AllTPCsInTPCsets[1] = {cryo.TPCPtr(1), cryo.TPCPtr(3)};
  readout::TPCsetDataContainer<TPCColl_t> TPCsetTPCs(1, TPCsetCount.size());
  
  for( auto&& [s, TPCPtrs ]: util::enumerate(AllTPCsInTPCsets) ){
    readout::TPCsetID const tpcsetid
    { cryo.ID(), static_cast<readout::TPCsetID::TPCsetID_t>(s) };
    TPCsetTPCs[tpcsetid] = TPCPtrs;
  }

  //  number of readout planes in each TPC set
  readout::TPCsetDataContainer<unsigned int> ROPcount( TPCsetTPCs.dimSize<0U>(), 
						       TPCsetTPCs.dimSize<1U>(), 0 );
  
  // geo::PlaneGeo objects in each readout plane
  readout::ROPDataContainer<PlaneColl_t> ROPplanes(TPCsetTPCs.dimSize<0U>(),
						   TPCsetTPCs.dimSize<1U>(), MaxROPs);
  
  for (auto [ c, nTPCsets ]: util::enumerate(TPCsetCount)) {
    readout::CryostatID const cryoid
    { static_cast<readout::CryostatID::CryostatID_t>(c) };
    for (readout::TPCsetID::TPCsetID_t s = 0; s < nTPCsets; ++s) {
      readout::TPCsetID const tpcsetid { cryoid, s };
      unsigned int nROP  = 3;
      ROPcount[tpcsetid] = 3;

      auto TPCs = TPCsetTPCs[tpcsetid];
      for( unsigned r = 0; r<nROP;++r ){ // assuming the number of rops == number of planes per TPC
	// get PlaneColl_t for each set for each plane
	PlaneColl_t planes;
	for (geo::TPCGeo const* TPC: TPCs){ 
	  planes.emplace_back( TPC->PlanePtr(r) );
	}
	readout::ROPID const ropid { tpcsetid, r };
	ROPplanes[ ropid ] = std::move( planes );
      } // ROPs/ planes
    } // tpcs in tpcset
  }//cryo, TPCSets
      
  
  //
  // the TPC set each TPC belongs to 
  geo::TPCDataContainer<readout::TPCsetID> TPCtoTPCset( 1, MaxTPCs, {});
  for (auto c
	 : util::counter<geo::CryostatID::CryostatID_t>(TPCsetTPCs.dimSize<0>())){
    
    geo::CryostatID const cid { c };
    for (auto s: util::counter<readout::TPCsetID::TPCsetID_t>(TPCsetTPCs.dimSize<1>()))
      {
	readout::TPCsetID const sid { cid, s };
      	for (geo::TPCGeo const* TPC: TPCsetTPCs[sid]) {
	  assert(TPC && TPC->ID());
	  TPCtoTPCset[TPC->ID()] = sid;
	  MF_LOG_TRACE(fLogCategory) << TPC->ID() << " => " << sid;
	} // for TPCs in TPC set
      } // for TPC sets in cryostat
  } // for cryostats
  
  
  // ROP each wire belongs to
  geo::PlaneDataContainer<readout::ROPID> PlaneToROP;
  for (auto c
    : util::counter<geo::CryostatID::CryostatID_t>(ROPplanes.dimSize<0>())
       )
    {
      geo::CryostatID const cid { c };
      for (auto s: util::counter<readout::TPCsetID::TPCsetID_t>(ROPplanes.dimSize<1>()))
	{
	  readout::TPCsetID const sid { cid, s };
	  for (auto r: util::counter<readout::ROPID::ROPID_t>(ROPplanes.dimSize<2>()))
	    {
	      readout::ROPID const rid { sid, r };
	      for (geo::PlaneGeo const* plane: ROPplanes[rid]) {
		assert(plane && plane->ID());
		fPlaneToROP[plane->ID()] = rid;
		MF_LOG_TRACE(fLogCategory) << plane->ID() << " => " << rid;
	      } // for planes in readout plane
	    } // for readout planes in TPC set
	} // for TPC sets in cryostat
    } // for cryostats
  
  
  // copy to ReadoutMapInfo struct
  fReadoutMapInfo.set(
		      std::move(results).TPCsetCount(), std::move(results).TPCsetTPCs(),
		      std::move(results).ROPcount(), std::move(results).ROPplanes(),
		      std::move(results).TPCtoTPCset(), std::move(results).PlaneToROP()
		      );
  
} // dune::ColdBoxChannelMapAlg::buildReadoutPlanes()


// -----------------------------------------------------------------------------
PlaneType_t dune::ColdBoxChannelMapAlg::findPlaneType(readout::ROPID const& rid) const
{
  /*
   * This implementation is very fragile, relying on the fact that the first
   * induction plane numbers are `kFirstInductionType`, the second induction
   * plane numbers are `kSecondInductionType` and the collection plane numbers
   * are `kCollectionType`. This assumption is not checked anywhere.
   * 
   */
  constexpr std::array PlaneTypes = { // let's C++ figure out type and size
    kFirstInductionType,  // P:0
    kSecondInductionType, // P:1
    kCollectionType       // P:2
  };
  
  PlaneColl_t const& planes = ROPplanes(rid);
  if (planes.empty()) return kUnknownType;
  if (auto const planeNo = planes.front()->ID().Plane; planeNo < PlaneTypes.size())
    return PlaneTypes[planeNo];
  else return kUnknownType;
  
} // dune::ColdBoxChannelMapAlg::findPlaneType()


// ----------------------------------------------------------------------------
geo::SigType_t dune::ColdBoxChannelMapAlg::SignalTypeForChannelImpl
  (raw::ChannelID_t const channel) const
{
  /*
   * We rely on the accuracy of `findPlaneType()` (which is admittedly less than
   * great) to assign signal type accordingly.
   */
  
  dune::details::ChannelToWireMap::ChannelsInROPStruct const* channelInfo
    = fChannelToWireMap.find(channel);
  if (!channelInfo) return geo::kMysteryType;
  
  switch (findPlaneType(channelInfo->ropid)) {
    case kFirstInductionType:
    case kSecondInductionType:
      return geo::kInduction;
    case kCollectionType:
      return geo::kCollection;
    default:
      return geo::kMysteryType;
  } // switch
  
} // dune::ColdBoxChannelMapAlg::SignalTypeForChannelImpl()


/*
// -----------------------------------------------------------------------------
auto dune::ColdBoxChannelMapAlg::extractWirelessChannelParams
  (Config::WirelessChannelStruct const& params) -> WirelessChannelCounts_t
{
  return {
    // even TPC sets (e.g. C:0 S:0)
    std::array{
      std::make_pair(
        params.FirstInductionPreChannels(),
        params.FirstInductionPostChannels()
        ),
      std::make_pair(
        params.SecondInductionEvenPreChannels(),
        params.SecondInductionEvenPostChannels()
        ),
      std::make_pair(
        params.CollectionEvenPreChannels(),
        params.CollectionEvenPostChannels()
        )
    },
    // odd TPC sets (e.g. C:0 S:1)
    std::array{
      std::make_pair(
        params.FirstInductionPreChannels(),
        params.FirstInductionPostChannels()
        ),
      std::make_pair(
        params.SecondInductionOddPreChannels(),
        params.SecondInductionOddPostChannels()
        ),
      std::make_pair(
        params.CollectionOddPreChannels(),
        params.CollectionOddPostChannels()
        )
    }
    };
  
} // dune::ColdBoxChannelMapAlg::extractWirelessChannelParams()
*/

// ----------------------------------------------------------------------------
std::string dune::ColdBoxChannelMapAlg::PlaneTypeName(PlaneType_t planeType) {
  
  using namespace std::string_literals;
  switch (planeType) {
    case kFirstInductionType:  return "first induction"s;
    case kSecondInductionType: return "second induction"s;
    case kCollectionType:      return "collection induction"s;
    case kUnknownType:         return "unknown"s;
    default:
      return "unsupported ("s + std::to_string(planeType) + ")"s;
  } // switch
  
} // dune::ColdBoxChannelMapAlg::PlaneTypeName()


// ----------------------------------------------------------------------------

