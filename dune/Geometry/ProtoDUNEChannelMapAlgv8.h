////////////////////////////////////////////////////////////////////////
/// \file  ProtoDUNEChannelMapAlgv8.h
/// \brief Interface to algorithm class for a specific detector channel mapping in protoDUNE
///
/// \version $Id:  $
/// \author  Laura Paulucci
///
/// 
/// June 2021
///
/// Inherits from the central DuneApaChannelMapAlg, just adding special channel-mapping
/// specific to the photon detectors in protoDUNE.
///
////////////////////////////////////////////////////////////////////////
#ifndef geo_ProtoDUNEChannelMapAlgv8_H
#define geo_ProtoDUNEChannelMapAlgv8_H

#include <map>
#include "cetlib_except/exception.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h" // readout::ROPID, ...
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/Geometry/DuneApaChannelMapAlg.h"


namespace geo{

  class ProtoDUNEChannelMapAlgv8 : public DuneApaChannelMapAlg {

  public:
  
    typedef unsigned int Index;
  
    ProtoDUNEChannelMapAlgv8(const fhicl::ParameterSet& pset);
  
    void Initialize(GeometryData_t const& geodata) override;
    void Uninitialize() override;


    // Photon detector methos being overridden here
    Index NOpChannels(Index NOpDets)                           const override;
    Index MaxOpChannel(Index NOpDets)                          const override;
    Index NOpHardwareChannels(Index opDet)                     const override;
    Index OpChannel(Index detNum, Index channel)               const override;
    Index OpDetFromOpChannel(Index opChannel)                  const override;
    Index HardwareChannelFromOpChannel(Index opChannel)        const override;
    bool  IsValidOpChannel(Index opChannel, Index /*NOpDets*/) const override
      {return fOpDet.count(opChannel); }
  
    // ProtoDUNE-specific
    Index SSPfromOpDet(Index opDet)                const;
    Index FirstSSPChFromOpDet(Index OpDet)         const;
    Index OpChannelFromSSP(Index ssp, Index sspch) const;
    void  PrintChannelMaps()                       const;

  private:
  
    Index fMaxOpChannel;
    Index fNOpChannels;
    std::map<Index, Index> fSSP;
    std::map<Index, Index> fSSPChOne;
    std::map<Index, Index> fOpDet;
    std::map<Index, Index> fHWChannel;

  
  }; // class ProtoDUNEChannelMapAlgv8
  
  //----------------------------------------------------------------------------
  
  
}
#endif

