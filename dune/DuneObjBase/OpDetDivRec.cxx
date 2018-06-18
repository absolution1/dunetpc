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
#include <vector>
#include <TObject.h>
#include "OpDetDivRec.h"

namespace sim {
  void OpDetDivRec::Clear(){
    tick_chans.clear(); 
  }
}
