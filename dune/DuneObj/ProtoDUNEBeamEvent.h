#ifndef BEAMDATA_PROTODUNEBEAM_H
#define BEAMDATA_PROTODUNEBEAM_H

#include <vector>
#include <bitset>
#include <iostream>
#include <map>
#include "lardataobj/RecoBase/Track.h"

namespace beam
{
  
  //Fiber Beam Monitor
  //
  struct FBM{

    //Bitmap for hit fibers in the monitor
//    std::bitset<192> fibers;
    std::array<short,192> fibers;

    //Raw Data from ifbeam
    double fiberData[6];    
    double timeData[4]; 

    double timeStamp;
    
    //ID (position in beamline?) of monitor
    int ID;
  };

  //Cerenkov Threshold Detector
  //
  struct CKov{          
    //Status at time of system trigger (on/off)
    bool trigger;
    
    long long int timeStamp;
  };



  class ProtoDUNEBeamEvent{
    public:
      ProtoDUNEBeamEvent();
      ~ProtoDUNEBeamEvent();

      void              InitFBMs(std::vector<std::string>);
      double            GetT0(size_t);
      size_t            GetNT0(){return t0.size();};
      void              AddT0(double);
     
      FBM               GetFBM(std::string name, size_t theTrigger);
      void              AddFBMTrigger(std::string, FBM); 
      void              ReplaceFBMTrigger(std::string, FBM, size_t); 
      void              DecodeFibers(std::string, size_t);
      double            DecodeFiberTime(std::string, size_t);
      std::array<double,4> ReturnTriggerAndTime(std::string, size_t);
      short             GetFiberStatus(std::string, size_t, size_t);
      std::vector<short> GetActiveFibers(std::string, size_t);
      double            GetFiberTime(std::string, size_t); 
      size_t            GetNFBMTriggers(std::string);
      std::bitset<32>   toBinary(double); 

      void              AddCKov0Trigger(CKov theCKov){ CKov0.push_back(theCKov); }; 
      void              AddCKov1Trigger(CKov theCKov){ CKov1.push_back(theCKov); };
      int               GetNCKov0Triggers(){ return CKov0.size(); };
      int               GetNCKov1Triggers(){ return CKov1.size(); };
      short             GetCKov0Status(size_t);
      short             GetCKov1Status(size_t);
      long long         GetCKov0Time(size_t);
      long long         GetCKov1Time(size_t);

      void              AddTOF0Trigger(double theT){ TOF0.push_back(theT); };
      void              AddTOF1Trigger(double theT){ TOF1.push_back(theT); }; 
      double            GetTOF0(size_t);
      double            GetTOF1(size_t);
      int               GetNTOF0Triggers(){ return TOF0.size(); };
      int               GetNTOF1Triggers(){ return TOF1.size(); };

      void              AddBeamTrack(recob::Track theTrack){ Tracks.push_back(theTrack);};
      recob::Track      GetBeamTrack(size_t i){ return Tracks[i];};
      size_t            GetNBeamTracks(){return Tracks.size();}
      const std::vector< recob::Track > & GetBeamTracks() const;
    private:

      //Time of a coincidence between 2 TOFs
      //Signalling a good particle
      //
      std::vector<double> t0;

      //Set of FBMs
      //Indices: [Monitor in beam]['event' in monitor]
//      std::vector< std::vector < FBM > > fiberMonitors;
     std::map<std::string, std::vector< FBM > > fiberMonitors;
      size_t nFBMs;

      //Set of TOF detectors
      //
      std::vector< double > TOF0;
      std::vector< double > TOF1;

      //Set of Cerenkov detectors
      //
      std::vector< CKov > CKov0;
      std::vector< CKov > CKov1;

      std::vector<recob::Track> Tracks;
  
  };

  inline const std::vector< recob::Track > & ProtoDUNEBeamEvent::GetBeamTracks() const { return Tracks; }
  
  inline double ProtoDUNEBeamEvent::GetT0(size_t trigger){
    if( trigger > t0.size() - 1 ){
      std::cout << "Error. Trigger out of bunds" << std::endl;
      return -1;
    }

    return t0[trigger];
  }
  inline void ProtoDUNEBeamEvent::AddT0(double theT0){
    t0.push_back(theT0);
  }

  ////////////Fiber Monitor Access
  inline FBM ProtoDUNEBeamEvent::GetFBM(std::string FBMName, size_t theTrigger){
    FBM dummy; 
    dummy.ID = -1;
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Error FBM not found" << std::endl;
//      for(size_t i = 0; i < fiberMonitors.find(FBMName); ++i){
        
      for(auto itF = fiberMonitors.begin(); itF != fiberMonitors.end(); ++itF){
        std::cout << "\t" << itF->first << std::endl; 
      }
      return dummy;
    }

    if(theTrigger > t0.size()){
      std::cout << "Error, trigger out of range" << std::endl;
      return dummy;
    }
    return fiberMonitors[FBMName].at(theTrigger); 
  }

  inline void ProtoDUNEBeamEvent::AddFBMTrigger(std::string FBMName, FBM theFBM){
//    if( (FBMName > (nFBMs - 1) ) ){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Error FBM not found" << std::endl;
//      for(size_t i = 0; i < fiberMonitors.find(FBMName); ++i){
        
      for(auto itF = fiberMonitors.begin(); itF != fiberMonitors.end(); ++itF){
        std::cout << "\t" << itF->first << std::endl; 
      }
      return;
    }

    //Check if it's the first time in the monitor. Replace dummy
    if(fiberMonitors[FBMName][0].ID == -1){
      std::cout << "Replacing dummy FBM" << std::endl;
      std::vector<FBM>::iterator theIt = fiberMonitors[FBMName].begin();
      fiberMonitors[FBMName].insert(theIt,theFBM);
      fiberMonitors[FBMName].pop_back();
    }
    else{
      fiberMonitors[FBMName].push_back(theFBM);
    }
  }

  inline void ProtoDUNEBeamEvent::ReplaceFBMTrigger(std::string FBMName, FBM theFBM, size_t theTrigger){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Error FBM not found" << std::endl;
        
      for(auto itF = fiberMonitors.begin(); itF != fiberMonitors.end(); ++itF){
        std::cout << "\t" << itF->first << std::endl; 
      }

      return;
    }

    if( theTrigger > t0.size() - 1 ){
      std::cout << "Trigger out of range of good particles." << std::endl;
      return;
    }
    
//    std::cout << "Replacing FBM Trigger" << std::endl;
    fiberMonitors[FBMName].at(theTrigger) = theFBM;

  }

  inline void ProtoDUNEBeamEvent::DecodeFibers(std::string FBMName, size_t nTrigger){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return;
    }
    if( (nTrigger > fiberMonitors[FBMName].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[FBMName].size() - 1 << "]" << std::endl;
      return;
    }
    
    for(int iSet = 0; iSet < 6; ++iSet){
      
      std::bitset<32> theseFibers = toBinary( fiberMonitors[FBMName][nTrigger].fiberData[iSet] );

      for(int  iFiber = 0; iFiber < 32; ++iFiber){      
        fiberMonitors[FBMName][nTrigger].fibers[iSet*32 + iFiber] = theseFibers[iFiber];
      }
    }
  }

  inline double ProtoDUNEBeamEvent::DecodeFiberTime(std::string FBMName, size_t nTrigger){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM not found in list" << std::endl;
      return -1.;
    }
    if( (nTrigger > fiberMonitors[FBMName].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[FBMName].size() - 1 << "]" << std::endl;
      return -1.;
    }

    
    //timeData[3] -> Event Time in Seconds
    //timeData[2] -> Number of 8ns ticks after that time
//    return fiberMonitors[FBMName][nTrigger].timeData[0];
    return fiberMonitors[FBMName][nTrigger].timeData[3] + fiberMonitors[FBMName][nTrigger].timeData[2]*8.e-9;
  }

  inline std::array<double,4> ProtoDUNEBeamEvent::ReturnTriggerAndTime(std::string FBMName, size_t nTrigger){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "FBM not found in list" << std::endl;
      return {{-1.,-1.,-1.,-1.}};
    }
    if( (nTrigger > fiberMonitors[FBMName].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[FBMName].size() - 1 << "]" << std::endl;
      return {{-1.,-1.,-1.,-1.}};
    }

    return {{fiberMonitors[FBMName][nTrigger].timeData[0], fiberMonitors[FBMName][nTrigger].timeData[1], fiberMonitors[FBMName][nTrigger].timeData[2], fiberMonitors[FBMName][nTrigger].timeData[3]}};
  }

  inline std::bitset<32> ProtoDUNEBeamEvent::toBinary(double num){
    std::bitset<64> mybits( (long(num)) );
    std::bitset<32> upper, lower;
    for(int i = 0; i < 32; ++i){
      lower[i] = mybits[i];
      upper[i] = mybits[i + 32];
    }
    if(upper.any()) std::cout << "WARNING: NONZERO HALF" << std::endl;

    return lower;
  }


  inline short ProtoDUNEBeamEvent::GetFiberStatus(std::string FBMName, size_t nTrigger, size_t iFiber){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1;          
    }
    if( (iFiber > 191)){
      std::cout << "Please input fiber in range [0,191]" << std::endl;
      return -1;
    }
    if( (nTrigger > fiberMonitors[FBMName].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[FBMName].size() - 1 << "]" << std::endl;
      return -1;
    }
    return fiberMonitors[FBMName][nTrigger].fibers[iFiber];
  }

  inline std::vector<short> ProtoDUNEBeamEvent::GetActiveFibers(std::string FBMName, size_t nTrigger){
    std::vector<short> active;

    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return active;          
    }
    if( (nTrigger > fiberMonitors[FBMName].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[FBMName].size() - 1 << "]" << std::endl;
      return active;
    }
    
    for(size_t iF = 0; iF < 192; ++iF){
      if(fiberMonitors[FBMName][nTrigger].fibers[iF]) active.push_back(iF); 
    }

    return active;
  }

  inline double ProtoDUNEBeamEvent::GetFiberTime(std::string FBMName, size_t nTrigger){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1;          
    }
    if( (nTrigger > fiberMonitors[FBMName].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[FBMName].size() - 1 << "]" << std::endl;
      return -1;
    }
    return fiberMonitors[FBMName][nTrigger].timeStamp;
  }

  inline size_t ProtoDUNEBeamEvent::GetNFBMTriggers(std::string FBMName){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1;          
    }
    return fiberMonitors[FBMName].size();
  }
  
  /////////////////////////////////


  ////////////Cerenkov Access
  
  inline short ProtoDUNEBeamEvent::GetCKov0Status(size_t nTrigger){
    if( (nTrigger >= CKov0.size()) ){
      std::cout << "Please input index in range [0," << CKov0.size() - 1 << "]" << std::endl;
      return -1;
    }
    
    return CKov0[nTrigger].trigger;                
  }

  inline short ProtoDUNEBeamEvent::GetCKov1Status(size_t nTrigger){
    if( (nTrigger >= CKov1.size()) ){
      std::cout << "Please input index in range [0," << CKov1.size() - 1 << "]" << std::endl;
      return -1;
    }
    
    return CKov1[nTrigger].trigger;                
  }

  inline long long ProtoDUNEBeamEvent::GetCKov0Time(size_t nTrigger){
    if( (nTrigger >= CKov0.size()) ){
      std::cout << "Please input index in range [0," << CKov0.size() - 1 << "]" << std::endl;
      return -1;
    }
    
    return CKov0[nTrigger].timeStamp;                
  }

  inline long long ProtoDUNEBeamEvent::GetCKov1Time(size_t nTrigger){
    if( (nTrigger >= CKov1.size()) ){
      std::cout << "Please input index in range [0," << CKov1.size() - 1 << "]" << std::endl;
      return -1;
    }
    
    return CKov1[nTrigger].timeStamp;                
  }
  /////////////////////////////////


  ////////////TOF Access
  inline double ProtoDUNEBeamEvent::GetTOF0(size_t nTrigger){
    if( (nTrigger >= TOF0.size()) ){
      std::cout << "Please input index in range [0," << TOF0.size() - 1 << "]" << std::endl;
      return -1;
    }

    return TOF0[nTrigger];
  }

  inline double ProtoDUNEBeamEvent::GetTOF1(size_t nTrigger){
    if( (nTrigger >= TOF1.size()) ){
      std::cout << "Please input index in range [0," << TOF1.size() - 1 << "]" << std::endl;
      return -1;
    }

    return TOF1[nTrigger];
  }
  /////////////////////////////////

  /*
      std::vector< long long int > TOF1;
      std::vector< long long int > TOF2;

      //Set of Cerenkov detectors
      //
      std::vector< CKov > CKov1;
      std::vector< CKov > CKov2;
  */
}


#endif
