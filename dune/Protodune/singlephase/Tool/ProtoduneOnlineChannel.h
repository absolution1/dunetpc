// ProtoduneOnlineChannel.h
//
// David Adams
// August 2018
//
// Tool that converts a protodune SP offline channel number
// to an online number with the convention that
//   chanOn = 128*IFMB + ICHF
//
// IFMB is a global FEMB number 6*IAPA + IFMB_APA
// where IAPA is the offline APA/TPS number and IFMB_APA
// is the geometric FEMB number in the APA running clockwise
// from TPC side to cryostat side.
//
// ICHF is the channel number in the FEMB 8*ASIC + ICHA_ASIC.
//
// The mapping is hardwired here.

#ifndef ProtoduneOnlineChannel_H
#define ProtoduneOnlineChannel_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"

class ProtoduneOnlineChannel : public IndexMapTool {

public:

  ProtoduneOnlineChannel(const fhicl::ParameterSet& ps);

  Index get(Index chanOff) const override;

private:

  Index m_LogLevel;

  Index nwirPlane[4] = {800, 800, 480, 480};
  Index nwirFemb[4] = {40, 40, 48, 48};
  Index uch[40];
  Index vch[40];
  Index zch[48];

};

DEFINE_ART_CLASS_TOOL(ProtoduneOnlineChannel)

#endif
