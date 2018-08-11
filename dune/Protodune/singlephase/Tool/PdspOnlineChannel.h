// PdspOnlineChannel.h
//
// David Adams
// June 2018
//
// Tool that converts a protodune SP offline channel number
// to an online number with the convention that
//   chanOn = 2560*(APA) + 512*WIB_in_APA + 128*FEMB_in_WIB + CHAN_in_FEMB
//
// From this, the APA, global FEMB and FEMB channel are
//           APA = chanOn/2560 + 1
//          FEMB = chanOn/128
//  CHAN_in_FEMB = chanOn % 128
//
// NOte that the FEMB ordering is most likely not geometric and that the
// indices are obtained from PdspChannelMapService which has unusual
// and varying values.

#ifndef PdspOnlineChannel_H
#define PdspOnlineChannel_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"

class PdspOnlineChannel : public IndexMapTool {

public:

  PdspOnlineChannel(const fhicl::ParameterSet& ps);

  Index get(Index chanOff) const override;

private:

  Index m_LogLevel;

};

DEFINE_ART_CLASS_TOOL(PdspOnlineChannel)

#endif
