////////////////////////////////////////////////////////////////////////
//
// File: OpDetDivRec..cxx
// Author: Jason Stock (jason.stock@mines.sdsmt.edu)
// This file is a set of classes for writeout of data relative to 
// OpticalChannels in the DUNE detectors.
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
  void OpDetDivRec::AddPhoton(int opchan, int tid, OpDet_Time_Chans::stored_time_t time){
    Time_Chans_t::iterator itr = priv_FindClosestTimeChan(time);
    if(itr == time_chans.end() || itr->time!=time){
      std::vector<Chan_Phot> cfl;
      cfl.emplace_back(opchan, tid);
      time_chans.emplace(itr, time, std::move(cfl));
      if(! std::is_sorted(time_chans.begin(), time_chans.end(), CompareByPdTime() ) ) //Just to guarantee no funny buisiness in the ordering. This shold generally not be an issue because we should always pass this check. After I convince myself it is always filled correctly, I will remove this check.
        std::sort(time_chans.begin(), time_chans.end(), CompareByPdTime());
    }else{
      for(auto cfp = itr->phots.begin(); cfp!= itr->phots.end(); ++cfp){
        if(cfp!= itr->phots.end() && cfp->opChan != opchan && cfp->trackID!=tid){
          continue;
        }else if(cfp == itr->phots.end()){
          itr->phots.emplace_back(opchan, tid);
        }else{
          cfp->AddPhoton();
          break;
        }
      }
    }
  }//End  AddPhoton

  std::vector<std::pair<int, double>> OpDetDivRec::GetFracs(OpDet_Time_Chans::stored_time_t time){
    std::vector<std::pair<int, double>> ret;
    auto itr = priv_FindClosestTimeChan(time);
    if( itr != time_chans.end() || itr->time==time ){
      ret = itr->GetFracs();
    }  
    return ret;
  }

  std::vector<std::pair<int, double>> OpDetDivRec::GetFracs(OpDet_Time_Chans::stored_time_t time, int tid){
    std::vector<std::pair<int, double>> ret;
    auto itr = priv_FindClosestTimeChan(time);
    if( itr != time_chans.end() || itr->time==time ){
      ret = itr->GetFracs(tid);
    }  
    return ret;
  }

  OpDetDivRec::Time_Chans_t::iterator OpDetDivRec::priv_FindClosestTimeChan(OpDet_Time_Chans::stored_time_t pdTime){
    return std::lower_bound
      (time_chans.begin(), time_chans.end(), pdTime, CompareByPdTime() );
  }
  OpDetDivRec::Time_Chans_t::const_iterator OpDetDivRec::priv_FindClosestTimeChan (OpDet_Time_Chans::stored_time_t pdTime) const{
    return std::lower_bound
      (time_chans.begin(), time_chans.end(), pdTime, CompareByPdTime() );
  }
  std::pair<OpDetDivRec::Time_Chans_t::const_iterator , bool> OpDetDivRec::FindClosestTimeChan(OpDet_Time_Chans::stored_time_t pdTime) const{
    auto ret = priv_FindClosestTimeChan(pdTime);
    bool found=false;
    if(ret!=time_chans.end()){
      found=true;
    }
    return std::make_pair(ret, found);
  }

  OpDetDivRec::time_slice OpDetDivRec::GetSlice(OpDet_Time_Chans::stored_time_t low_time, OpDet_Time_Chans::stored_time_t high_time){
    OpDetDivRec::time_slice ret;
    ret.lower = priv_FindClosestTimeChan(low_time);
    ret.upper = priv_FindClosestTimeChan(high_time);
    if(ret.lower!=ret.upper && ret.upper < time_chans.end()){ 
      return ret;
    }else{
      ret.lower=time_chans.end();
      ret.upper=time_chans.end();
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
