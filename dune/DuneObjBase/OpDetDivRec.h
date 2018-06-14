////////////////////////////////////////////////////////////////////////
//
// File: EventRecord.h
// Author: Jason Stock (jason.stock@mines.sdsmt.edu)
// This file is a set of classes for writeout of data relative to calibration sources.
// The main purpose of this class is to allow easy fast analysis of calibration
// simulations.
//
////////////////////////////////////////////////////////////////////////

//includes
#ifndef DUNE_DUNEOBJBASE_SIM_H
#define DUNE_DUNEOBJBASE_SIM_H
#include <vector>
#include <TObject.h>

namespace sim {

  class TickDivRec{
    int Photons;
    int crosstalk;
  };

  class ChanDivRec{
    UInt_t OpChanNum;
    std::vector<sim::TickDivRec>;//Should this be index matched? I guess so. After all, this is an expansion record, but this feels dirty to me.
  };

  class OpDetDivRec{
    UInt_t opDetNum;
    std::vector<sim::ChanDivRec> chans;
  };

}
