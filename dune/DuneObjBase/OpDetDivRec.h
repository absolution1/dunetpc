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
#include <utility>
#include <algorithm>

namespace sim {
  //these records already exist tied to optical detectors (as they are made tied to the BTRs. No need to be OpDet wise.
  //Remember these are not per sdp smart, just per opchan smart. The normalization is per opchan at each unit time. We will assume the SDPs are split uniformly according to how everything over the time tick was split.

  struct Chan_Phot{
    int opChan;
    double phot;
    Chan_Phot():
      opChan(-1)
    {}
    Chan_Phot(int opChanIn):
      opChan(opChanIn),
      phot(1.0)
    {}
    void AddPhoton(){
      phot+=1.0;
    }
  };
  struct OpDet_Time_Chans{
    //    std::pair<double, std:vector<Chan_Phot>>
    typedef double stored_time_t;
    stored_time_t time;
    OpDet_Time_Chans()
    {}
    OpDet_Time_Chans(stored_time_t& timeIn);
    OpDet_Time_Chans(stored_time_t& timeIn, std::vector<Chan_Phot> photIn);
    std::vector<Chan_Phot> phots;

    std::vector<std::pair<int, double>> GetFracs() const{
      double total=0.0;
      std::vector<std::pair<int, double>> ret;
      for(auto& a : phots){ total += a.phot;}
      for(auto& a : phots){ ret.emplace_back(a.opChan, (a.phot / total));}
      return ret;
    }
    double GetFrac(int opChanIn) const{
      double ret=0.0;
      std::vector<std::pair<int, double>> fracs = GetFracs();
      for(auto pair : fracs){
        if(pair.first==opChanIn){
          ret=pair.second;
          break;
        }
      }
      return ret;
    }
  };

  class OpDetDivRec{
    public:
      typedef std::vector<OpDet_Time_Chans> Tick_Chans_t;
      struct time_slice{
        Tick_Chans_t::const_iterator lower;
        Tick_Chans_t::const_iterator upper;
      };

    private:
      int fOpDetNum; //Move this to be private.
      Tick_Chans_t tick_chans; //Move this to private //This must be filled with emplace to keep it sorted.
    public:
      OpDetDivRec();
      explicit OpDetDivRec(int det);
      Tick_Chans_t const& GetTickChans(){return tick_chans;}
      int OpDetNum() const{ return fOpDetNum; }
      void add_photon(int opchan, OpDet_Time_Chans::stored_time_t pdTime);
      std::vector<std::pair<int, double>> GetFracs(OpDet_Time_Chans::stored_time_t time);
      time_slice GetSlice(OpDet_Time_Chans::stored_time_t low_time, OpDet_Time_Chans::stored_time_t high_time);

      //      double GetFras(OpDet_Time_Chans::stored_time_t, OpChan);
      //      Tick_Chans_t::iterator priv_FindClosestTimeChan(const OpDet_Time_Chans::stored_time_t& pdTime);
      std::pair<OpDetDivRec::Tick_Chans_t::const_iterator, bool> FindClosestTimeChan( OpDet_Time_Chans::stored_time_t pdTime) const;
      void clear();
      template <typename Stream>
        void Dump(Stream&& out, std::string indent, std::string first_indent) const;

      /// Documentation at `Dump(Stream&&, std::string, std::string) const`.
      template <typename Stream>
        void Dump(Stream&& out, std::string indent = "") const
        { Dump(std::forward<Stream>(out), indent, indent); }

    private:
      //  struct CompareByPdTime ;
      struct CompareByPdTime {
        bool operator()
          (OpDet_Time_Chans const& a, OpDet_Time_Chans const& b) const
          {return a.time<b.time;}
        bool operator()
          (OpDet_Time_Chans::stored_time_t const& a, OpDet_Time_Chans const& b)
          {return a<b.time;}
        bool operator()
          (OpDet_Time_Chans const& a, OpDet_Time_Chans::stored_time_t const& b )
          {return a.time<b;}
      };

    private:
      Tick_Chans_t::iterator priv_FindClosestTimeChan(OpDet_Time_Chans::stored_time_t pdTime);
      Tick_Chans_t::const_iterator priv_FindClosestTimeChan( OpDet_Time_Chans::stored_time_t pdTime) const;

  };

}

// -----------------------------------------------------------------------------
// ---  template implementation
// ---
template <class Stream>
void sim::OpDetDivRec::Dump
(Stream&& out, std::string indent, std::string first_indent) const
{
  out << first_indent << "OpDet #" << OpDetNum() << " read " << tick_chans.size()
    << " tick_chans:\n";
  for (const auto& tc: tick_chans) {
    auto time = tc.time;
    out << indent << "  time " << time
      << " with " << tc.phots.size() << " Photon deposits\n";
    for (const auto& photr: tc.phots) {
      out << indent
        << "OpChan: "<<photr.opChan <<" with "<<photr.phot<<" photon detction records\n";
    } // for SDPs
  } // for timePDclocks
} // sim::OpDetBacktrackerRecord::Dump<>()


////////////////////////////////////////////////////////////////////////
#endif //DUNE_DUNEOBJBASE_SIM_H
