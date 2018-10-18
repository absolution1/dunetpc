#ifndef BEAMDATA_PROTODUNEBEAMEVENT_H
#define BEAMDATA_PROTODUNEBEAMEVENT_H

#include <vector>
#include <bitset>
#include <iostream>
#include <map>
#include "lardataobj/RecoBase/Track.h"
#include "dune/DuneObj/ProtoDUNEBeamSpill.h"

namespace beam
{
  

  class ProtoDUNEBeamEvent{
    public:
      ProtoDUNEBeamEvent(){};
     ~ProtoDUNEBeamEvent(){};

      std::pair< double, double >  GetT0(){ return t0;};

      double            GetFullT0(){ return t0.first + 1.e-9*t0.second; };
      void              SetT0(std::pair< double, double > theT0){ t0 = theT0; };
     
      FBM               GetFBM(std::string);
      void              SetFBMTrigger(std::string, FBM); 

      void               DecodeFibers(std::string);
      double             DecodeFiberTime(std::string, double);
      short              GetFiberStatus(std::string, size_t);
      std::vector<short> GetActiveFibers(std::string);
      double             GetFiberTime(std::string); 
      std::bitset<32>    toBinary(double num){return std::bitset<32>( (uint32_t(num)) );}; 

      void              SetCKov0(CKov theCKov){ CKov0 = theCKov; }; 
      void              SetCKov1(CKov theCKov){ CKov1 = theCKov; };
      short             GetCKov0Status(){ return CKov0.trigger; };
      short             GetCKov1Status(){ return CKov1.trigger; };
      double            GetCKov0Time(){ return CKov0.timeStamp; };
      double            GetCKov1Time(){ return CKov1.timeStamp; };
      double            GetCKov0Pressure(){ return CKov0.pressure; };
      double            GetCKov1Pressure(){ return CKov1.pressure; };


      void              SetTOF0Trigger( std::pair<double,double> theT){ TOF0 = theT; };
      void              SetTOF1Trigger( std::pair<double,double> theT){ TOF1 = theT; }; 
      void              SetTOFChan    ( int theChan )                 { TOFChan = theChan; };

      std::pair< double, double > GetTOF0() { return TOF0; };
      std::pair< double, double > GetTOF1() { return TOF1; };

      double GetFullTOF0() { return TOF0.first + 1.e-9*TOF0.second; };
      double GetFullTOF1() { return TOF1.first + 1.e-9*TOF1.second; };

      double            GetTOF(){ return ( (TOF1.first  - TOF0.first) 
                                 + 1.e-9 * (TOF1.second - TOF0.second) ); };

      int               GetTOFChan(){ return TOFChan; };

      void              AddBeamTrack(recob::Track theTrack){ Tracks.push_back(theTrack);};
      recob::Track      GetBeamTrack(size_t i){ return Tracks[i];};
      size_t            GetNBeamTracks(){return Tracks.size();}
      const std::vector< recob::Track > & GetBeamTracks() const;
      void              ClearBeamTracks(){ Tracks.clear(); };

      void              SetActiveTrigger(size_t theTrigger){ activeTrigger = theTrigger; isMatched = true; };
      bool              CheckIsMatched(){ return isMatched; };
      void              SetUnmatched(){ isMatched = false; };
      size_t            GetActiveTrigger(){ return activeTrigger; };

      void              AddRecoBeamMomentum( double theMomentum ){ RecoBeamMomenta.push_back( theMomentum ); };
      size_t            GetNRecoBeamMomenta()                    { return RecoBeamMomenta.size(); };
      double            GetRecoBeamMomentum( size_t i )          { return RecoBeamMomenta[i]; };  
      void              ClearRecoBeamMomenta()                   { RecoBeamMomenta.clear(); }; 

      void              SetBITrigger(int theTrigger){ BITrigger = theTrigger; };
      int               GetBITrigger(){ return BITrigger; };
      
      void              SetSpillStart(double theSpillStart){ SpillStart = theSpillStart; };
      double            GetSpillStart(){ return SpillStart; };

      void              SetSpillOffset(double theSpillOffset){ SpillOffset = theSpillOffset; };
      double            GetSpillOffset(){ return SpillOffset; };

      void              SetCTBTimestamp(double theCTBTimestamp){ CTBTimestamp = theCTBTimestamp; };
      double            GetCTBTimestamp(){ return CTBTimestamp; };

    private:

      //Time of a coincidence between 2 TOFs
      //Signalling a good particle
      //Known as 'GeneralTrigger'
      //
      std::pair<double,double> t0;

      //Timestamp from the CTB signaling a 
      //Good particle signal was received
      //
      double CTBTimestamp;

      //Set of FBMs
      //Indices: [Monitor in beam]
      std::map<std::string, FBM > fiberMonitors;

      //Set of TOF detectors
      //
      std::pair< double, double > TOF0;
      std::pair< double, double > TOF1;
      int TOFChan;

      //Set of Cerenkov detectors
      //
      CKov CKov0;
      CKov CKov1;

      std::vector<recob::Track> Tracks;
  
      size_t activeTrigger;
      bool   isMatched = false;

      std::vector< double > RecoBeamMomenta;

      int BITrigger;
      double SpillStart;
      double SpillOffset;
  };

  inline const std::vector< recob::Track > & ProtoDUNEBeamEvent::GetBeamTracks() const { return Tracks; }
  

  ////////////Fiber Monitor Access
  inline FBM ProtoDUNEBeamEvent::GetFBM(std::string FBMName){
    FBM dummy; 
    dummy.ID = -1;
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
    
      std::cout << "Error FBM " << FBMName << " not found" << std::endl;
        
      for(auto itF = fiberMonitors.begin(); itF != fiberMonitors.end(); ++itF){
        std::cout << "\t" << itF->first << std::endl; 
      }

      return dummy;
    }

    return fiberMonitors[FBMName]; 
  }

  inline void ProtoDUNEBeamEvent::SetFBMTrigger(std::string FBMName, FBM theFBM){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM not found" << std::endl;
      std::cout << "Setting New FBM " << FBMName << std::endl;
    }

    fiberMonitors[FBMName] = theFBM;
  }



  inline void ProtoDUNEBeamEvent::DecodeFibers(std::string FBMName){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " Not Found"  << std::endl;
      return;
    }
    
    for(int iSet = 0; iSet < 6; ++iSet){
      
      std::bitset<32> theseFibers = toBinary( fiberMonitors[FBMName].fiberData[iSet] );

      for(int  iFiber = 0; iFiber < 32; ++iFiber){      
        fiberMonitors[FBMName].fibers[iSet*32 + iFiber] = theseFibers[iFiber];
      }
    }
  }

  inline double ProtoDUNEBeamEvent::DecodeFiberTime(std::string FBMName, double OffsetTAI){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " not found in list" << std::endl;
      return -1.;
    }
    
    return fiberMonitors[FBMName].timeData[3] - OffsetTAI + fiberMonitors[FBMName].timeData[2]*8.e-9;
  }


  inline short ProtoDUNEBeamEvent::GetFiberStatus(std::string FBMName, size_t iFiber){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " not found in list" << std::endl;
      return -1;          
    }
    if( (iFiber > 191)){
      std::cout << "Please input fiber in range [0,191]" << std::endl;
      return -1;
    }

    return fiberMonitors[FBMName].fibers[iFiber];
  }

  inline std::vector<short> ProtoDUNEBeamEvent::GetActiveFibers(std::string FBMName){
    std::vector<short> active;

    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " not found in list" << std::endl;
      return active;          
    }
    
    for(size_t iF = 0; iF < 192; ++iF){
      if(fiberMonitors[FBMName].fibers[iF]) active.push_back(iF); 
    }

    return active;
  }

  inline double ProtoDUNEBeamEvent::GetFiberTime(std::string FBMName){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM " << FBMName << " not found in list" << std::endl;
      return -1;          
    }

    return fiberMonitors[FBMName].timeStamp;
  }

  /////////////////////////////////


}


#endif
