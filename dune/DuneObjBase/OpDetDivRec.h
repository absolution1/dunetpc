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
#include <map>

namespace sim {

  class ChanDivRec{
    public:
      double tick_photons_frac=0.0;
  };

/*  class TickDivRec{
    public:
      //double sdp_time;//This is a redundancy, but this is an evolving data type, and I am just trying to get it running as quickly as possible for now.
      //No need to know the time. We are index linked to the SDPs times.
      std::map<UInt_t, ChanDivRec> chans;
  };*/

  class OpDetDivRec{
    public:
      //    UInt_t opDetNum; //This is a redundancy
      typedef std::map<UInt_t, ChanDivRec> TickDivChans; //OpChan, ChanDivRec
      std::vector<std::pair<double, TickDivChans>> tick_chans; //Because the btss timePDclockSDPMap gets sorted, we will also need to be able to sort this to make the indicies match. That means we also have to include the time (bad memory redundancy, but in this version unavoidable), so we can also sort this one.
      void Clear();
  };

}
#endif //DUNE_DUNEOBJBASE_SIM_H
