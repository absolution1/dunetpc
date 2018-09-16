////////////////////////////////////////////////////////////////////////
//
// File: OpDetDivRec.h
// Author: Jason Stock (jason.stock@mines.sdsmt.edu)
// This file is a set of classes for writeout of data relative to
// OpticalChannels in the DUNE detectors.
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
  //Remember these are not per sdp smart, just per opchan smart. The normalization is per opchan at each unit time. We will assume the SDPs are split uniformly according to how everything over the time time was split.

  struct Chan_Phot{
    int opChan;
    int trackID;
    double phot;
    Chan_Phot():
      opChan(-1)
    {}
    Chan_Phot(int opChanIn, int tid=-999999):
      opChan(opChanIn),
      trackID(tid),
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

    std::vector<std::pair<int, double>> GetFracs(int tid) const{
      double total=0.0;
      std::vector<std::pair<int, double>> ret;
      for(auto& a : phots){ total += (a.trackID==tid)?a.phot:0;}
      for(auto& a : phots){ ret.emplace_back(a.opChan, (a.phot / total));}
      return ret;
    }
    std::vector<std::pair<int, double>> GetFracs() const{
      double total=0.0;
      std::vector<std::pair<int, double>> ret;
      for(auto& a : phots){ total += a.phot;}
      for(auto& a : phots){ ret.emplace_back(a.opChan, (a.phot / total));}
      return ret;
    }
    double GetFrac(int opChanIn, int tid) const{
      double ret=0.0;
      std::vector<std::pair<int, double>> fracs = GetFracs(tid);
      for(auto pair : fracs){
        if(pair.first==opChanIn){
          ret=pair.second;
          break;
        }
      }
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
      typedef std::vector<OpDet_Time_Chans> Time_Chans_t;
      struct time_slice{
        Time_Chans_t::const_iterator lower;
        Time_Chans_t::const_iterator upper;
      };

    private:
      int fOpDetNum; //Move this to be private.
      Time_Chans_t time_chans; //Move this to private //This must be filled with emplace to keep it sorted.
    public:
      OpDetDivRec();
      OpDetDivRec(int det);
      Time_Chans_t const& GetTimeChans(){return time_chans;}
      int OpDetNum() const{ return fOpDetNum; }
      void AddPhoton(int opchan, int tid, OpDet_Time_Chans::stored_time_t pdTime);
      std::vector<std::pair<int, double>> GetFracs(OpDet_Time_Chans::stored_time_t time);
      std::vector<std::pair<int, double>> GetFracs(OpDet_Time_Chans::stored_time_t time, int tid);
      time_slice GetSlice(OpDet_Time_Chans::stored_time_t low_time, OpDet_Time_Chans::stored_time_t high_time);

      //      double GetFras(OpDet_Time_Chans::stored_time_t, OpChan);
      //      Time_Chans_t::iterator priv_FindClosestTimeChan(const OpDet_Time_Chans::stored_time_t& pdTime);
      std::pair<OpDetDivRec::Time_Chans_t::const_iterator, bool> FindClosestTimeChan( OpDet_Time_Chans::stored_time_t pdTime) const;
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
      Time_Chans_t::iterator priv_FindClosestTimeChan(OpDet_Time_Chans::stored_time_t pdTime);
      Time_Chans_t::const_iterator priv_FindClosestTimeChan( OpDet_Time_Chans::stored_time_t pdTime) const;

  };

}

// -----------------------------------------------------------------------------
// ---  template implementation
// ---
template <class Stream>
void sim::OpDetDivRec::Dump
(Stream&& out, std::string indent, std::string first_indent) const
{
  out << first_indent << "OpDet #" << OpDetNum() << " read " << time_chans.size()
    << " time_chans:\n";
  for (const auto& tc: time_chans) {
    auto time = tc.time;
    out << indent << "  time " << time
      << " with " << tc.phots.size() << " Photon deposits\n";
    for (const auto& photr: tc.phots) {
      out << indent
        << "OpChan: "<<photr.opChan <<" with "<<photr.phot<<" photons from particle with TrackID"<<photr.trackID<<"\n";
    } // for SDPs
  } // for timePDclocks
} // sim::OpDetBacktrackerRecord::Dump<>()


////////////////////////////////////////////////////////////////////////
#endif //DUNE_DUNEOBJBASE_SIM_H
