#ifndef BEAMDATA_PROTODUNEBEAM_H
#define BEAMDATA_PROTODUNEBEAM_H

#include <vector>
#include <bitset>
#include <iostream>
#include <map>

namespace beam
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



  class ProtoDUNEBeamEvent{
    public:
      ProtoDUNEBeamEvent();
      ~ProtoDUNEBeamEvent();

      void              InitFBMs(size_t);
      long long         GetT0(){ return t0; };
      
      void              AddFBMTrigger(size_t, FBM); 
      void              DecodeFibers(size_t, size_t);
      double            DecodeFiberTime(size_t, size_t);
      short             GetFiberStatus(size_t, size_t, size_t);
      std::vector<short> GetActiveFibers(size_t, size_t);
      long long         GetFiberTime(size_t, size_t); 
      size_t            GetNFBMTriggers(size_t);
      std::bitset<32>   toBinary(double); 

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
//     std::map<std::string, std::vector< FBM > > fiberMonitors;
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
  inline void ProtoDUNEBeamEvent::AddFBMTrigger(size_t iMonitor, FBM theFBM){
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

/* inline void ProtoDUNEBeamEvent::AddFBMTrigger(std::string FBMName, FBM theFBM){
    
    if(!fiberMonitors.count(FBMName)){
      std::cout << "Initializing FBM vector" << std::endl;
      //fiberMonitors[FBMName] = std::vector< FBM >();     
      fiberMonitors.insert(std::map<std::string, std::vector<FBM> >::value_type(FBMName, std::vector<FBM>() )); 
      std::cout << "Done" << std::endl;
    }
//    std::cout << "Size: " << fiberMonitors[FBMName].size() << std::endl;
    fiberMonitors[FBMName].push_back(theFBM);
    std::cout << "Added" << std::endl;
  }*/

  inline void ProtoDUNEBeamEvent::DecodeFibers(size_t iMonitor, size_t nTrigger){
    if( (iMonitor > (fiberMonitors.size() -1) ) ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return;
    }
    if( (nTrigger > fiberMonitors[iMonitor].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[iMonitor].size() - 1 << "]" << std::endl;
      return;
    }
    
    for(int iSet = 0; iSet < 6; ++iSet){
      
      std::bitset<32> theseFibers = toBinary( fiberMonitors[iMonitor][nTrigger].fiberData[iSet] );

      for(int  iFiber = 0; iFiber < 32; ++iFiber){      
        fiberMonitors[iMonitor][nTrigger].fibers[iSet*32 + iFiber] = theseFibers[iFiber];
      }
    }
  }

  inline double ProtoDUNEBeamEvent::DecodeFiberTime(size_t iMonitor, size_t nTrigger){
    if( (iMonitor > (fiberMonitors.size() -1) ) ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1.;
    }
    if( (nTrigger > fiberMonitors[iMonitor].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[iMonitor].size() - 1 << "]" << std::endl;
      return -1.;
    }

    //FOR NOW JUST RETURN THE TRIGGER
    return fiberMonitors[iMonitor][nTrigger].timeData[0];
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


  inline short ProtoDUNEBeamEvent::GetFiberStatus(size_t iMonitor, size_t nTrigger, size_t iFiber){
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

  inline std::vector<short> ProtoDUNEBeamEvent::GetActiveFibers(size_t iMonitor, size_t nTrigger){
    std::vector<short> active;

    if( (iMonitor > (fiberMonitors.size() - 1)) ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return active;          
    }
    if( (nTrigger > fiberMonitors[iMonitor].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[iMonitor].size() - 1 << "]" << std::endl;
      return active;
    }
    
    for(size_t iF = 0; iF < 192; ++iF){
      if(fiberMonitors[iMonitor][nTrigger].fibers[iF]) active.push_back(iF); 
    }

    return active;
  }

  inline long long ProtoDUNEBeamEvent::GetFiberTime(size_t iMonitor, size_t nTrigger){
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

  inline size_t ProtoDUNEBeamEvent::GetNFBMTriggers(size_t iMonitor){
    if( (iMonitor > (fiberMonitors.size() - 1)) ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1;          
    }
    return fiberMonitors[iMonitor].size();
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
  inline long long ProtoDUNEBeamEvent::GetTOF0(size_t nTrigger){
    if( (nTrigger >= TOF0.size()) ){
      std::cout << "Please input index in range [0," << TOF0.size() - 1 << "]" << std::endl;
      return -1;
    }

    return TOF0[nTrigger];
  }

  inline long long ProtoDUNEBeamEvent::GetTOF1(size_t nTrigger){
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
