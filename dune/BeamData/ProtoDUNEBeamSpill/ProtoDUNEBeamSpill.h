#ifndef BEAMDATA_PROTODUNEBEAMSPILL_H
#define BEAMDATA_PROTODUNEBEAMSPILL_H

#include <vector>
#include <bitset>
#include <iostream>

namespace beamspill
{
  
  //Fiber Beam Monitor
  //
  struct FBM{

    //Bitmap for hit fibers in the monitor
    std::bitset<192> fibers;

    //Raw Data from ifbeam
    double fiberData[6];    
    double timeData[4]; 

    long long int timeStamp;
    
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



  class ProtoDUNEBeamSpill{
    public:
      ProtoDUNEBeamSpill();
     ~ProtoDUNEBeamSpill();

      void              InitFBMs(size_t);
      long long         GetT0(){ return t0; };
      
      void              AddFBMTrigger(size_t, FBM); 
      short             GetFiberStatus(size_t, size_t, size_t);
      long long         GetFiberTime(size_t, size_t); 
      int               GetNFBMTriggers(size_t);
      
      void              AddCKov0Trigger(CKov theCKov){ CKov0.push_back(theCKov); }; 
      void              AddCKov1Trigger(CKov theCKov){ CKov1.push_back(theCKov); };
      int               GetNCKov0Triggers(){ return CKov0.size(); };
      int               GetNCKov1Triggers(){ return CKov1.size(); };
      short             GetCKov0Status(size_t);
      short             GetCKov1Status(size_t);
      long long         GetCKov0Time(size_t);
      long long         GetCKov1Time(size_t);

      void              AddTOF0Trigger(long long theT){ TOF0.push_back(theT); };
      void              AddTOF1Trigger(long long theT){ TOF1.push_back(theT); }; 
      long long         GetTOF0(size_t);
      long long         GetTOF1(size_t);
      int               GetNTOF0Triggers(){ return TOF0.size(); };
      int               GetNTOF1Triggers(){ return TOF1.size(); };

    
    private:

      //First time of anything in the spill
      //
      long long int t0;

      //Set of FBMs
      //Indices: [Monitor in beam]['event' in monitor]
      std::vector< std::vector < FBM > > fiberMonitors;
      size_t nFBMs;

      //Set of TOF detectors
      //
      std::vector< long long int > TOF0;
      std::vector< long long int > TOF1;

      //Set of Cerenkov detectors
      //
      std::vector< CKov > CKov0;
      std::vector< CKov > CKov1;
  
  };



  ////////////Fiber Monitor Access
  inline void ProtoDUNEBeamSpill::AddFBMTrigger(size_t iMonitor, FBM theFBM){
    if( (iMonitor > (nFBMs - 1) ) ){
      std::cout << "Error FBM index out of range" << std::endl;
      return;
    }

    //Check if it's the first time in the monitor. Replace dummy
    if(fiberMonitors[iMonitor][0].ID == -1){
      std::cout << "Replacing dummy FBM" << std::endl;
      std::vector<FBM>::iterator theIt = fiberMonitors[iMonitor].begin();
      fiberMonitors[iMonitor].insert(theIt,theFBM);
      fiberMonitors[iMonitor].pop_back();
    }
    else{
      fiberMonitors[iMonitor].push_back(theFBM);
    }
  }

  inline short ProtoDUNEBeamSpill::GetFiberStatus(size_t iMonitor, size_t nTrigger, size_t iFiber){
    if( (iMonitor > (fiberMonitors.size() - 1)) ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1;          
    }
    if( (iFiber > 191)){
      std::cout << "Please input fiber in range [0,191]" << std::endl;
      return -1;
    }
    if( (nTrigger > fiberMonitors[iMonitor].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[iMonitor].size() - 1 << "]" << std::endl;
      return -1;
    }
    return fiberMonitors[iMonitor][nTrigger].fibers[iFiber];
  }

  inline long long ProtoDUNEBeamSpill::GetFiberTime(size_t iMonitor, size_t nTrigger){
    if(iMonitor > fiberMonitors.size() - 1){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1;          
    }
    if( (nTrigger > fiberMonitors[iMonitor].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[iMonitor].size() - 1 << "]" << std::endl;
      return -1;
    }
    return fiberMonitors[iMonitor][nTrigger].timeStamp;
  }

  inline int ProtoDUNEBeamSpill::GetNFBMTriggers(size_t iMonitor){
    if( (iMonitor > (fiberMonitors.size() - 1)) ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1;          
    }
    return fiberMonitors[iMonitor].size();
  }
  /////////////////////////////////


  ////////////Cerenkov Access
  
  inline short ProtoDUNEBeamSpill::GetCKov0Status(size_t nTrigger){
    if( (nTrigger >= CKov0.size()) ){
      std::cout << "Please input index in range [0," << CKov0.size() - 1 << "]" << std::endl;
      return -1;
    }
    
    return CKov0[nTrigger].trigger;                
  }

  inline short ProtoDUNEBeamSpill::GetCKov1Status(size_t nTrigger){
    if( (nTrigger >= CKov1.size()) ){
      std::cout << "Please input index in range [0," << CKov1.size() - 1 << "]" << std::endl;
      return -1;
    }
    
    return CKov1[nTrigger].trigger;                
  }

  inline long long ProtoDUNEBeamSpill::GetCKov0Time(size_t nTrigger){
    if( (nTrigger >= CKov0.size()) ){
      std::cout << "Please input index in range [0," << CKov0.size() - 1 << "]" << std::endl;
      return -1;
    }
    
    return CKov0[nTrigger].timeStamp;                
  }

  inline long long ProtoDUNEBeamSpill::GetCKov1Time(size_t nTrigger){
    if( (nTrigger >= CKov1.size()) ){
      std::cout << "Please input index in range [0," << CKov1.size() - 1 << "]" << std::endl;
      return -1;
    }
    
    return CKov1[nTrigger].timeStamp;                
  }
  /////////////////////////////////


  ////////////TOF Access
  inline long long ProtoDUNEBeamSpill::GetTOF0(size_t nTrigger){
    if( (nTrigger >= TOF0.size()) ){
      std::cout << "Please input index in range [0," << TOF0.size() - 1 << "]" << std::endl;
      return -1;
    }

    return TOF0[nTrigger];
  }

  inline long long ProtoDUNEBeamSpill::GetTOF1(size_t nTrigger){
    if( (nTrigger >= TOF1.size()) ){
      std::cout << "Please input index in range [0," << TOF1.size() - 1 << "]" << std::endl;
      return -1;
    }

    return TOF1[nTrigger];
  }
  /////////////////////////////////

      std::vector< long long int > TOF1;
      std::vector< long long int > TOF2;

      //Set of Cerenkov detectors
      //
      std::vector< CKov > CKov1;
      std::vector< CKov > CKov2;
  
}


#endif
