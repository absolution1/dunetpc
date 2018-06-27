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
#include <algorithm>

namespace sim {
  OpDetDivRec::OpDetDivRec():
    fOpDetNum(-1)
  { }
  OpDetDivRec::OpDetDivRec(int det):
    fOpDetNum(det)
  { }

  //int OpDetDivRec::OpDetNum() const { return fOpDetNum; }
  void OpDetDivRec::add_photon(int opchan, OpDet_Time_Chans::stored_time_t time){
    Tick_Chans_t::iterator itr = priv_FindClosestTimeChan(time);
    if(itr == tick_chans.end() || itr->time!=time){
      std::vector<Chan_Phot> cfl;
      cfl.emplace_back(opchan);
      tick_chans.emplace(itr, time, std::move(cfl));
      if(! std::is_sorted(tick_chans.begin(), tick_chans.end(), CompareByPdTime() ) ) //Just to guarantee no funny buisiness in the ordering. This shold generally not be an issue because we should always pass this check. After I convince myself it is always filled correctly, I will remove this check.
        std::sort(tick_chans.begin(), tick_chans.end(), CompareByPdTime());
    }else{
      //for(std::vector<Chan_Phot>::iterator const& cfp : itr->phots){
      for(auto cfp = itr->phots.begin(); cfp!= itr->phots.end(); ++cfp){
        if(cfp!= itr->phots.end() && cfp->opChan != opchan){
          continue;
        }else if(cfp == itr->phots.end()){
          itr->phots.emplace_back(opchan);
        }else{
          cfp->AddPhoton();
        }
      }
    }
    }//End  add_photon

    std::vector<std::pair<int, double>> OpDetDivRec::GetFracs(OpDet_Time_Chans::stored_time_t time){
      std::vector<std::pair<int, double>> ret;
      auto itr = priv_FindClosestTimeChan(time);
      if( itr != tick_chans.end() || itr->time==time ){
        ret = itr->GetFracs();
      }  
      return ret;

    }

    OpDetDivRec::Tick_Chans_t::iterator OpDetDivRec::priv_FindClosestTimeChan(OpDet_Time_Chans::stored_time_t pdTime){
      return std::lower_bound
        (tick_chans.begin(), tick_chans.end(), pdTime, CompareByPdTime() );
    }
    OpDetDivRec::Tick_Chans_t::const_iterator OpDetDivRec::priv_FindClosestTimeChan (OpDet_Time_Chans::stored_time_t pdTime) const{
      return std::lower_bound
        (tick_chans.begin(), tick_chans.end(), pdTime, CompareByPdTime() );
    }
    std::pair<OpDetDivRec::Tick_Chans_t::const_iterator , bool> OpDetDivRec::FindClosestTimeChan(OpDet_Time_Chans::stored_time_t pdTime) const{
      auto ret = priv_FindClosestTimeChan(pdTime);
      bool found=false;
      if(ret!=tick_chans.end()){
        found=true;
      }
      return std::make_pair(ret, found);
    }

    OpDetDivRec::time_slice OpDetDivRec::GetSlice(OpDet_Time_Chans::stored_time_t low_time, OpDet_Time_Chans::stored_time_t high_time){
      OpDetDivRec::time_slice ret;
      ret.lower = priv_FindClosestTimeChan(low_time);
      ret.upper = priv_FindClosestTimeChan(high_time);
      if(ret.lower!=ret.upper && ret.upper < tick_chans.end()){ 
        return ret;
      }else{
        ret.lower=tick_chans.end();
        ret.upper=tick_chans.end();
        return ret;
      }
    }

    OpDet_Time_Chans::OpDet_Time_Chans(stored_time_t& timeIn):
      time(timeIn)
    {}
    OpDet_Time_Chans::OpDet_Time_Chans(stored_time_t& timeIn, std::vector<Chan_Phot> inVec):
      time(timeIn),
      phots(inVec)
    {}
  }
