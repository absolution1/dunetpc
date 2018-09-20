// PdspOnlineChannel.h
//
// David Adams
// June 2018
//
// Tool that converts a protodune SP offline channel number
// to an online index following one of three conventions:
//
// WIB ordered:
//   chanOn = 2560*KAPA + 512*KWIB + 128*KCON + FCHAN
//
// Connector ordered:
//   chanOn = 2560*KAPA + 640*KCON + 128*KWIB + FCHAN
//
// FEMB ordered:
//   chanOn = 2560*KAPA + 128*KFMB + FCHAN
//
// The conventions for KWIB, KCON and KFMB are those shown in the cabling
// figure at https://wiki.dunescience.org/wiki/ProtoDUNE_geometry except
// all are shifted by one to be zero based and the APA digit is dropped
// from the FEMB identifiers. I.e. if the labels in that figure are IWIB,
// ICON and IFMB, then the indices here may be expressed as
//   KWIB = IWIB - 1
//   KCON = ICON -1
//   KFMB = (FMB%100) - 1
//
// PdspChannelMapService is used to obtain these indices.
//
// Here KAPA follows the offline convention, i.e. viewed from above
//   1 3 5
//   0 2 4   beam --->
//
// FCHAN is the channel number in the FEMB (0-127) 16*ASIC = ASIC_channel
// assuming the mapping in DUNE DocDB 4064.
//
// It is expected that FEMB ordering will be used most often so that
//     KAPA = chanOn/2560    identifies APAs (TPC sets) as indicated above
//     FEMB = chanOn/128     is a global FEMB identifier
//     KFMB = chanOn % 128   is the channe # in the FEMB
//
// Configuration parameters:
//   LogLevel - 0=silent, 1=init messages, 2=message every call
//   Ordering - String indicatin the ordering: WIB, connector or FEMB

#ifndef PdspOnlineChannel_H
#define PdspOnlineChannel_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"

class PdspOnlineChannel : public IndexMapTool {

public:

  using Name = std::string;

  PdspOnlineChannel(const fhicl::ParameterSet& ps);

  Index get(Index chanOff) const override;

private:

  // Configuration parameters.
  Index m_LogLevel;
  Name m_Ordering;

  // Derived from configuration.
  bool m_orderByWib;
  bool m_orderByConnector;
  bool m_orderByFemb;

};

DEFINE_ART_CLASS_TOOL(PdspOnlineChannel)

#endif
