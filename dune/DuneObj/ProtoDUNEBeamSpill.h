#ifndef BEAMDATA_PROTODUNEBEAMSPILL_H
#define BEAMDATA_PROTODUNEBEAMSPILL_H

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

    std::vector<short> active;
    bool decoded;
  };

  //Cerenkov Threshold Detector
  //
  struct CKov{          
    //Status at time of system trigger (on/off)
    short trigger;
    double pressure;
    double timeStamp;
  };



  class ProtoDUNEBeamSpill{
    public:
      ProtoDUNEBeamSpill();
      ~ProtoDUNEBeamSpill();

      void              InitFBMs(std::vector<std::string>);
      std::pair< double, double >  GetT0(size_t);
      double            GetT0Sec(size_t);
      double            GetT0Nano(size_t);
      size_t            GetNT0(){return t0.size();};
      void              AddT0(std::pair< double, double >);
     
      FBM               GetFBM(std::string name, size_t theTrigger);
      void              AddFBMTrigger(std::string, FBM); 
      void              ReplaceFBMTrigger(std::string, FBM, size_t); 
      void              DecodeFibers(std::string, size_t);
      double            DecodeFiberTime(std::string, size_t, double);
      std::array<double,4> ReturnTriggerAndTime(std::string, size_t);
      short             GetFiberStatus(std::string, size_t, size_t);
      std::vector<short> GetActiveFibers(std::string, size_t);
      double            GetFiberTime(std::string, size_t); 
      size_t            GetNFBMTriggers(std::string);
      std::bitset<32>   toBinary(double); 

//      void              SetCKov0(CKov theCKov){ CKov0 = theCKov; }; 
//      void              SetCKov1(CKov theCKov){ CKov1 = theCKov; };
      void              AddCKov0(CKov theCKov){ CKov0.push_back(theCKov); }; 
      void              AddCKov1(CKov theCKov){ CKov1.push_back(theCKov); };
      CKov              GetCKov0(size_t i){ return CKov0[i]; }; 
      CKov              GetCKov1(size_t i){ return CKov1[i]; };
      size_t            GetNCKov0(){ return CKov0.size(); };
      size_t            GetNCKov1(){ return CKov1.size(); };


      short                         GetCKov0Status(size_t i){ return CKov0[i].trigger; };
      short                         GetCKov1Status(size_t i){ return CKov1[i].trigger; };
      double                        GetCKov0Pressure(size_t i){ return CKov0[i].pressure; };
      double                        GetCKov1Pressure(size_t i){ return CKov1[i].pressure; };
      double                        GetCKov0Time(size_t i){ return CKov0[i].timeStamp; };
      double                        GetCKov1Time(size_t i){ return CKov1[i].timeStamp; };


      void              AddTOF0Trigger( std::pair<double,double> theT){ TOF0.push_back(theT); };
      void              AddTOF1Trigger( std::pair<double,double> theT){ TOF1.push_back(theT); }; 
      void              AddTOFChan(int theChan){ TOFChan.push_back(theChan); };
      std::pair< double, double > GetTOF0(size_t);
      double                      GetTOF0Sec(size_t);
      double                      GetTOF0Nano(size_t);
      std::pair< double, double > GetTOF1(size_t);
      double                      GetTOF1Sec(size_t);
      double                      GetTOF1Nano(size_t);
      double            GetTOF( size_t );
      int               GetTOFChan(size_t);
      int               GetNTOF0Triggers(){ return TOF0.size(); };
      int               GetNTOF1Triggers(){ return TOF1.size(); };

      std::vector< double > GetMultipleTOFs( size_t i ){ return TOFs[i]; };
      std::vector< int    > GetMultipleTOFChans( size_t i ){ return TOFChans[i]; };
      std::vector< size_t > GetUpstreamTriggers( size_t i ){ return UpstreamTriggers[i]; };
      std::vector< size_t > GetDownstreamTriggers( size_t i ){ return DownstreamTriggers[i]; };

      void AddMultipleTOFs( std::vector< double > theTOFs){ TOFs.push_back( theTOFs ); };
      void AddMultipleTOFChans( std::vector< int > theChans){ TOFChans.push_back( theChans ); };
      void AddUpstreamTriggers( std::vector< size_t > theTriggers ){ UpstreamTriggers.push_back( theTriggers ); };
      void AddDownstreamTriggers( std::vector< size_t > theTriggers ){ DownstreamTriggers.push_back( theTriggers ); };

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
      size_t            GetNRecoBeamMomenta() { return RecoBeamMomenta.size(); };
      double            GetRecoBeamMomentum( size_t i ){ return RecoBeamMomenta[i]; };  
      void              ClearRecoBeamMomenta() { RecoBeamMomenta.clear(); }; 

      void              SetBITrigger(int theTrigger){ BITrigger = theTrigger; };
      int               GetBITrigger(){ return BITrigger; };
      
      void              SetSpillStart(double theSpillStart){ SpillStart = theSpillStart; };
      double            GetSpillStart(){ return SpillStart; };

      void              SetSpillOffset(double theSpillOffset){ SpillOffset = theSpillOffset; };
      double            GetSpillOffset(){ return SpillOffset; };

      void              SetCTBTimestamp(double theCTBTimestamp){ CTBTimestamp = theCTBTimestamp; };
      double            GetCTBTimestamp(){ return CTBTimestamp; };

      void              SetMagnetCurrent(double theMagnetCurrent){ MagnetCurrent = theMagnetCurrent; };
      double            GetMagnetCurrent(){ return MagnetCurrent; };


    private:

      //Time of a coincidence between 2 TOFs
      //Signalling a good particle
      //
      std::vector<std::pair<double,double>> t0;

      //Timestamp from the CTB signaling a 
      //Good particle signal was received
      //
      double CTBTimestamp;

      //Set of FBMs
      //Indices: [Monitor in beam]['event' in monitor]
//      std::vector< std::vector < FBM > > fiberMonitors;
      std::map<std::string, std::vector< FBM > > fiberMonitors;
      size_t nFBMs;

      //Set of TOF detectors
      //
      std::vector< std::pair< double, double > > TOF0;
      std::vector< std::pair< double, double > > TOF1;
      std::vector< int > TOFChan;

      //Possible combinations of TOFs
      std::vector< std::vector< double > > TOFs;
      std::vector< std::vector< int    > > TOFChans;
      std::vector< std::vector< size_t > > UpstreamTriggers;
      std::vector< std::vector< size_t > > DownstreamTriggers;

      //Set of Cerenkov detectors
      //
      std::vector< CKov > CKov0;
      std::vector< CKov > CKov1;

      std::vector<recob::Track> Tracks;
  
      size_t activeTrigger;
      bool   isMatched = false;

      std::vector< double > RecoBeamMomenta;

      int BITrigger;
      double SpillStart;
      double SpillOffset;

      double MagnetCurrent;
  };

  inline const std::vector< recob::Track > & ProtoDUNEBeamSpill::GetBeamTracks() const { return Tracks; }
  
  inline std::pair< double, double > ProtoDUNEBeamSpill::GetT0(size_t trigger){
    if( trigger > t0.size() - 1 ){
      std::cout << "Error. Trigger out of bunds" << std::endl;
      return std::make_pair(-1.,-1.);
    }

    return t0[trigger];
  }

  inline double ProtoDUNEBeamSpill::GetT0Sec(size_t trigger){
    if( trigger > t0.size() - 1 ){
      std::cout << "Error. Trigger out of bunds" << std::endl;
      return -1.;
    }

    return t0[trigger].first;
  }

  inline double ProtoDUNEBeamSpill::GetT0Nano(size_t trigger){
    if( trigger > t0.size() - 1 ){
      std::cout << "Error. Trigger out of bunds" << std::endl;
      return -1.;
    }

    return t0[trigger].second;
  }

  inline void ProtoDUNEBeamSpill::AddT0(std::pair< double, double > theT0){
    t0.push_back(theT0);
  }

  ////////////Fiber Monitor Access
  inline FBM ProtoDUNEBeamSpill::GetFBM(std::string FBMName, size_t theTrigger){
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

  inline void ProtoDUNEBeamSpill::AddFBMTrigger(std::string FBMName, FBM theFBM){
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

  inline void ProtoDUNEBeamSpill::ReplaceFBMTrigger(std::string FBMName, FBM theFBM, size_t theTrigger){
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

  inline void ProtoDUNEBeamSpill::DecodeFibers(std::string FBMName, size_t nTrigger){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return;
    }
    if( (nTrigger > fiberMonitors[FBMName].size()) ){
      std::cout << "Please input trigger in range [0," << fiberMonitors[FBMName].size() - 1 << "]" << std::endl;
      return;
    }

    //This always clears the currently active fibers in the FBM.
    fiberMonitors[FBMName][nTrigger].active.clear();
    
    for(int iSet = 0; iSet < 6; ++iSet){
      
      std::bitset<32> theseFibers = toBinary( fiberMonitors[FBMName][nTrigger].fiberData[iSet] );

      for(int  iFiber = 0; iFiber < 32; ++iFiber){      
        fiberMonitors[FBMName][nTrigger].fibers[iSet*32 + iFiber] = theseFibers[iFiber];
        if(theseFibers[iFiber]) fiberMonitors[FBMName][nTrigger].active.push_back(iSet*32 + iFiber);        
      }
    }

    fiberMonitors[FBMName][nTrigger].decoded = true;
  }

  inline double ProtoDUNEBeamSpill::DecodeFiberTime(std::string FBMName, size_t nTrigger, double OffsetTAI){
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
    return fiberMonitors[FBMName][nTrigger].timeData[3] - OffsetTAI + fiberMonitors[FBMName][nTrigger].timeData[2]*8.e-9;
  }

  inline std::array<double,4> ProtoDUNEBeamSpill::ReturnTriggerAndTime(std::string FBMName, size_t nTrigger){
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

  inline std::bitset<32> ProtoDUNEBeamSpill::toBinary(double num){
    std::bitset<64> mybits( (long(num)) );
    std::bitset<32> upper, lower;
    for(int i = 0; i < 32; ++i){
      lower[i] = mybits[i];
      upper[i] = mybits[i + 32];
    }
    if(upper.any()) std::cout << "WARNING: NONZERO HALF" << std::endl;

    return lower;
  }


//  inline std::bitset<32> ProtoDUNEBeamSpill::toBinary(double num){
//    return std::bitset<32>( (uint32_t(num)) );
//  }


  inline short ProtoDUNEBeamSpill::GetFiberStatus(std::string FBMName, size_t nTrigger, size_t iFiber){
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

  inline std::vector<short> ProtoDUNEBeamSpill::GetActiveFibers(std::string FBMName, size_t nTrigger){
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

  inline double ProtoDUNEBeamSpill::GetFiberTime(std::string FBMName, size_t nTrigger){
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

  inline size_t ProtoDUNEBeamSpill::GetNFBMTriggers(std::string FBMName){
    if( fiberMonitors.find(FBMName) == fiberMonitors.end() ){
      std::cout << "Please input monitor in range [0," << fiberMonitors.size() - 1 << "]" << std::endl;
      return -1;          
    }
    return fiberMonitors[FBMName].size();
  }
  /////////////////////////////////

  ////////////TOF Access
  inline std::pair< double, double > ProtoDUNEBeamSpill::GetTOF0(size_t nTrigger){
    if( (nTrigger >= TOF0.size()) ){
      std::cout << "Please input index in range [0," << TOF0.size() - 1 << "]" << std::endl;
      return std::make_pair(-1.,-1.);
    }

    return TOF0[nTrigger];
  }

  inline double ProtoDUNEBeamSpill::GetTOF0Sec(size_t nTrigger){
    if( (nTrigger >= TOF0.size()) ){
      std::cout << "Please input index in range [0," << TOF0.size() - 1 << "]" << std::endl;
      return -1.;
    }

    return TOF0[nTrigger].first;
  }

  inline double ProtoDUNEBeamSpill::GetTOF0Nano(size_t nTrigger){
    if( (nTrigger >= TOF0.size()) ){
      std::cout << "Please input index in range [0," << TOF0.size() - 1 << "]" << std::endl;
      return -1.;
    }

    return TOF0[nTrigger].second;
  }

  inline std::pair< double, double > ProtoDUNEBeamSpill::GetTOF1(size_t nTrigger){
    if( (nTrigger >= TOF1.size()) ){
      std::cout << "Please input index in range [0," << TOF1.size() - 1 << "]" << std::endl;
      return std::make_pair(-1.,-1.);
    }

    return TOF1[nTrigger];
  }

  inline double ProtoDUNEBeamSpill::GetTOF1Sec(size_t nTrigger){
    if( (nTrigger >= TOF1.size()) ){
      std::cout << "Please input index in range [0," << TOF1.size() - 1 << "]" << std::endl;
      return -1.;
    }

    return TOF1[nTrigger].first;
  }

  inline double ProtoDUNEBeamSpill::GetTOF1Nano(size_t nTrigger){
    if( (nTrigger >= TOF1.size()) ){
      std::cout << "Please input index in range [0," << TOF1.size() - 1 << "]" << std::endl;
      return -1.;
    }

    return TOF1[nTrigger].second;
  }

  inline int ProtoDUNEBeamSpill::GetTOFChan(size_t nTrigger){
    if( (nTrigger >= TOF1.size()) ){
      std::cout << "Please input index in range [0," << TOF1.size() - 1 << "]" << std::endl;
      return -1;
    }

    return TOFChan[nTrigger];
  }

  inline double ProtoDUNEBeamSpill::GetTOF(size_t nTrigger){
    if( (nTrigger >= TOF1.size()) ){
      std::cout << "Please input index in range [0," << TOF1.size() - 1 << "]" << std::endl;
      return -1;
    }

    return (TOF1[nTrigger].first - TOF0[nTrigger].first + 1.e-9*(TOF1[nTrigger].second - TOF0[nTrigger].second));
  }
  /////////////////////////////////

}


#endif
