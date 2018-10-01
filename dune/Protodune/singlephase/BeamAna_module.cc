////////////////////////////////////////////////////////////////////////
// Class:       BeamAna
// Plugin Type: producer (art v2_08_03)
// File:        BeamAna_module.cc
//
// Generated at Thu Nov  2 22:57:41 2017 by Jonathan Paley using cetskelgen
// from cetlib version v3_01_01.
//
// Edited by Jake Calcutt
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "IFBeam_service.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
//#include "dune/BeamData/ProtoDUNEBeamSpill/ProtoDUNEBeamSpill.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include <bitset>
#include <iomanip>
#include <utility>

#include "TTree.h"
#include "TH2F.h"
#include "TVectorD.h"
#include "TPolyMarker.h"

namespace proto {
  class BeamAna;
}


class proto::BeamAna : public art::EDProducer {
public:
  explicit BeamAna(fhicl::ParameterSet const & p);
  //  virtual ~BeamAna();

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamAna(BeamAna const &) = delete;
  BeamAna(BeamAna &&) = delete;
  BeamAna & operator = (BeamAna const &) = delete;
  BeamAna & operator = (BeamAna &&) = delete;

  // Required functions.
  void reconfigure(fhicl::ParameterSet const & p);
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;


  std::bitset<sizeof(double)*CHAR_BIT> toBinary(const long num);  

  TVector3 ConvertProfCoordinates(double x, double y, double z, double zOffset);
  void BeamMonitorBasisVectors(); 
  void RotateMonitorVector(TVector3 &vec); 

  void MakeTrack(size_t);
//  void matchCurvedTriggers();
  void GetPairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);
  void GetPairedStraightFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);
  void GetUnpairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);
  double GetPosition(std::string, int);
  TVector3 ProjectToTPC(TVector3, TVector3);
  double GetPairedPosition(std::string, size_t);
  TVector3 TranslateDeviceToDetector(TVector3);
 
  void  InitXBPFInfo(beam::ProtoDUNEBeamEvent *);
  void  parseGeneralXBPF(std::string, uint64_t, size_t);
  void  parseXBPF(uint64_t);
  void  parsePairedXBPF(uint64_t);
  void  parsePairedStraightXBPF(uint64_t);

  void  parseXTOF(uint64_t);
  void  parseXCET(uint64_t);

  void  InitXCETInfo(beam::ProtoDUNEBeamEvent *);

  template<class T> 
  T FetchWithRetries(uint64_t, std::string, int);
   
private:
  
  TTree * fOutTree;
  TTree * fTrackTree;

  std::vector<double> * trackX;
  std::vector<double> * trackY;
  std::vector<double> * trackZ;

  TH2F * fFirstBeamProf2D;
  TH2F * fSecondBeamProf2D;
  std::vector<TH2F *> fBeamProf2D;
  std::map<std::string, TH1F *> fBeamProf1D;
  std::map<std::string, std::vector<short> *> fActiveFibers;
  std::map<std::string, double> fProfTime;
  std::map<std::string, double> fProfTrigger1;
  std::map<std::string, double> fProfTrigger2;
  std::map<std::string, double> fProfTime1;
  std::map<std::string, double> fProfTime2;
  std::map<std::string, TTree*> fProfTree;
  TTree * fGenTrigTree;
  TTree * fXTOF1ATree;
  TTree * fXTOF1BTree;
  TTree * fXTOF2ATree;
  TTree * fXTOF2BTree;
  TTree * fXBTF022702Tree;
  TTree * fMatchedTriggers;
  double matchedGen;
  double matchedTOF1;
  double matchedTOF2;
  std::map<std::string, double> matchedXBPF;
  

  double fGenTrigFrac;
  double fGenTrigCoarse; 
  double fXTOF1AFrac;
  double fXTOF1ACoarse; 
  double fXTOF1BFrac;
  double fXTOF1BCoarse; 
  double fXTOF2AFrac;
  double fXTOF2ACoarse; 
  double fXTOF2BFrac;
  double fXTOF2BCoarse; 
  double fXBTF022702Frac;
  double fXBTF022702Coarse; 
  TH1F * fTOFHist;
  TH1F * fCKovHist;
  recob::Track * theTrack;
  std::vector<recob::Track*> theTracks;
  double eventTime;
  
  TVector3 fBMBasisX;
  TVector3 fBMBasisY;
  TVector3 fBMBasisZ;

  // Declare member data here.
  //bool fLoadFromDB; // unused
  double  fTimeWindow;
  double fTolerance;
  std::string fCSVFileName;
  std::string fBundleName;
  std::string fOutputLabel;
  std::string fInputLabel;
  double fDummyEventTime;
  std::string fURLStr;
  double fValidWindow;
  uint64_t fFixedTime;
  std::vector< uint64_t > fMultipleTimes;
  std::vector< std::string > fDevices;
  std::vector< std::pair<std::string, std::string> > fPairedDevices;
  std::vector< std::pair<std::string, std::string> > fPairedStraightDevices;
  std::map<std::string, std::string > fDeviceTypes;
//  std::vector< std::array<double, 3> > fCoordinates; 
  std::map< std::string, std::array<double,3> > fCoordinates;
  std::map< std::string, std::array<double, 3> > fRotations;
  TVector3 fGlobalDetCoords;
  std::array<double,3> fDetRotation;
  std::map< std::string, double > fFiberDimension;

  int fNRetries;

  // Names of the CERN Beam Devices
  // Names have the form of a prefix + device
  
  // Prefixes for different devices classes:
  // Beam Positions (BPF)
  // Time of flights (TOF)
  // and Cerenkov counters (CET)
  std::string fXBPFPrefix;
  std::string fXTOFPrefix;
  std::string fXCETPrefix;

  // Device Names for each device

  // Time of Flights
  std::string fTOF1;
  std::string fTOF1A, fTOF1B;
  std::string fTOF2;
  std::string fTOF2A, fTOF2B;

  // Cerenkovs
  std::string fCKov1;
  std::string fCKov2;

  double fRotateMonitorXZ;
  double fRotateMonitorYZ;

  double fFirstTrackingProfZ;
  double fSecondTrackingProfZ;
  double fNP04FrontZ;  
  double fBeamX, fBeamY, fBeamZ;

  beam::ProtoDUNEBeamEvent * beamevt;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
};

// Constructor
proto::BeamAna::BeamAna(fhicl::ParameterSet const & p)
{
  // Declare products this module will provide
  produces<std::vector<beam::ProtoDUNEBeamEvent>>();  

  // Configure/Reconfigure
  this->reconfigure(p);
}
// END Constructor
////////////////////////

////////////////////////
// Fetch Method
template <class T> 
T proto::BeamAna::FetchWithRetries(uint64_t time, std::string name, int nRetry){
  T theResult;
  
  uint64_t newTime;
  //Search around time given with nRetries
  for(newTime = time - nRetry; newTime < time + nRetry; ++newTime){
    std::cout << "Trying to grab from folder: " << name << std::endl;
    std::cout << "At Time: " << newTime << std::endl;    
    try{
      theResult = bfp->GetNamedVector(newTime, name);
      std::cout << "Successfully fetched" << std::endl;
      return theResult;
    }
    catch(std::exception e){
      std::cout << "Could not fetch with time " << newTime << std::endl;      
    }
  }
  
  //Try a final time. Let it crash if it doesn't work
  std::cout << "Trying a final time to grab from folder: " << name << std::endl;
  std::cout << "At time: " << newTime << std::endl;
  theResult = bfp->GetNamedVector(newTime, name);
  std::cout << "Successfully fetched" << std::endl;

  return theResult; 
}
// END FeatchWithRetries
////////////////////////

////////////////////////
// Producer Method (reads in the event and derives values)
void proto::BeamAna::produce(art::Event & e)
{
  BeamMonitorBasisVectors();

  //const auto theInfo = e.getValidHandle< std::vector<raw::RDTimeStamp> >("timingrawdecoder");  
//  const auto theInfo = e.getValidHandle< std::vector<raw::ctb::pdspctb> >(fInputLabel);  
//  std::cout << theInfo << std::endl;

  // Open up and read from  the IFBeam Service
  std::cerr << "%%%%%%%%%% Getting ifbeam service handle %%%%%%%%%%" << std::endl;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifb;
  std::cerr << "%%%%%%%%%% Got ifbeam handle %%%%%%%%%%" << std::endl;

  // Read in and cache the beam bundle folder for a specific time
  bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
  std::cerr << "%%%%%%%%%% Got beam folder %%%%%%%%%%" << std::endl;

  // Set the readout window of interest
  bfp->setValidWindow(fValidWindow);
  std::cerr << "%%%%%%%%%% Valid window " << bfp->getValidWindow() << " %%%%%%%%%%" << std::endl;

  // Now for the current event retrieve the time (trigger) of the event
  // This will take the form a 64-bit timestamp with a high and low word
  // Break this up to have access to each word and the long long word
  std::cout <<"Event Time: " << uint64_t(e.time().timeLow()) << std::endl;
  std::cout << "Low: " << e.time().timeLow()  << std::endl;
  std::cout << "High: " << e.time().timeHigh()  << std::endl;

  // Use the bottom word of the time
  eventTime = e.time().timeLow();
  bool usedEventTime = false;
  //Use multiple times provided to fcl
  if( fMultipleTimes.size() ){
    std::cout << "Using multiple times: " << std::endl;
    for(size_t it = 0; it < fMultipleTimes.size(); ++it){
      std::cout << fMultipleTimes[it] << std::endl;
    } // endfor
  }
  //Use singular time provided to fcl
  else if(fFixedTime > 0.){
    std::cout << "Using singular time: " << fFixedTime << std::endl;
    fMultipleTimes.push_back(fFixedTime);
  }
  //Use time of event
  else{
    std::cout <<" Using Event Time: " << uint64_t(e.time().timeLow()) << std::endl;
    fMultipleTimes.push_back(uint64_t(e.time().timeLow()));
    usedEventTime = true;
  } // endif

  //Start getting beam event info

  // Create a new beam event (note the "new" here)  
  beamevt = new beam::ProtoDUNEBeamEvent();

  // Loop over the different particle times 
  for(size_t it = 0; it < fMultipleTimes.size(); ++it){
    std::cout << "Time: " << fMultipleTimes[it] << std::endl;

    // Parse the Time of Flight Counter data for the list
    // of times that we are using
    parseXTOF(fMultipleTimes[it]);
    std::cout << "NGoodParticles: " << beamevt->GetNT0()           << std::endl;
    std::cout << "NTOF0: "          << beamevt->GetNTOF0Triggers() << std::endl;
    std::cout << "NTOF1: "          << beamevt->GetNTOF1Triggers() << std::endl;

    // Parse the Beam postion counter information for the list
    // of time that we are using
    InitXBPFInfo(beamevt);
    parseXBPF(fMultipleTimes[it]);
    parsePairedXBPF(fMultipleTimes[it]);
    parsePairedStraightXBPF(fMultipleTimes[it]);

    std::cout << "NXBPF: " << beamevt->GetNFBMTriggers(fDevices[0]) << std::endl;

    // Loop over the number of TOF counter readouts (i.e. GetNT0())
    for(size_t ip = 0; ip < beamevt->GetNT0(); ++ip){
      std::cout << beamevt->GetT0(ip)   << " "
		<< beamevt->GetTOF0(ip) << " "
		<< beamevt->GetTOF1(ip) << " "
		<< beamevt->GetFiberTime(fDevices[0],ip) << std::endl;

      // Associate the times from TOF-0 and TOF-1 with the master time T0
      matchedGen  = beamevt->GetT0(ip);
      matchedTOF1 = beamevt->GetTOF0(ip);
      matchedTOF2 = beamevt->GetTOF1(ip);

      for(std::map<std::string,double>::iterator itMap = matchedXBPF.begin(); itMap != matchedXBPF.end(); ++itMap){
        itMap->second = beamevt->GetFiberTime(itMap->first, ip);
      }
      
      fMatchedTriggers->Fill();

    }

    for(size_t iTrigger = 0; iTrigger < 10; ++iTrigger){
      MakeTrack(iTrigger);

    }

    for(size_t iTrack = 0; iTrack < theTracks.size(); ++iTrack){    
      beamevt->AddBeamTrack( *(theTracks[iTrack]) ); 
      
      auto thisTrack = theTracks[iTrack]; 
      const recob::TrackTrajectory & theTraj = thisTrack->Trajectory();
      trackX->clear(); 
      trackY->clear(); 
      trackZ->clear(); 
      std::cout << "npositison: " << theTraj.NPoints() << std::endl;
      for(auto const & pos : theTraj.Trajectory().Positions()){
        std::cout << pos.X() << " " << pos.Y() << " " << pos.Z() << std::endl;
        trackX->push_back(pos.X());
        trackY->push_back(pos.Y());
        trackZ->push_back(pos.Z());
      }

      fTrackTree->Fill();
    }

    std::cout << "Added " << beamevt->GetNBeamTracks() << " tracks to the beam event" << std::endl;

//    parseXCET(fMultipleTimes[it]);
  }


 //Setting some dummy Triggers for drawing  


 /*std::map<std::string,double> dummyTriggerTimeLSB = {{"XBPF022707",1.50000e+08},{"XBPF022708",1.50000e+08},{"XBPF022716",1.50002e+08},{"XBPF022717",1.50002e+08}};
 std::map<std::string,double> dummyTriggerTimeMSB = {{"XBPF022707",1.53191e+09},{"XBPF022708",1.53191e+09},{"XBPF022716",1.53191e+09},{"XBPF022717",1.53191e+09}};
 std::map<std::string,double> dummyEventTimeLSB   = {{"XBPF022707",1.50000e+08},{"XBPF022708",1.50000e+08},{"XBPF022716",1.50000e+08},{"XBPF022717",1.50000e+08}};
 std::map<std::string,double> dummyEventTimeMSB   = {{"XBPF022707",1.53191e+09},{"XBPF022708",1.53191e+09},{"XBPF022716",1.53191e+09},{"XBPF022717",1.53191e+09}};


 std::map<std::string,std::array<double,6>> dummyHit     = { {"XBPF022707",{1,0,0,0,0,0}},
                                                            {"XBPF022708",{1,0,0,0,0,0}},
                                                            {"XBPF022716",{0,0,0,0,1,0}},
                                                            {"XBPF022717",{0,0,0,0,1,0}} }; 

 for(size_t ip = 0; ip < fPairedStraightDevices.size(); ++ip){
   beam::FBM fbm;
   fbm.ID = ip;
   std::string name = fPairedStraightDevices[ip].first; 
   fbm.timeData[0] = dummyTriggerTimeLSB[name];
   fbm.timeData[1] = dummyTriggerTimeMSB[name];
   fbm.timeData[2] = dummyEventTimeLSB[name];
   fbm.timeData[3] = dummyEventTimeMSB[name];

   for(int i = 0; i < 6; ++i){
     fbm.fiberData[i] = dummyHit[name][i];
   }

   beamevt->AddFBMTrigger(name,fbm);
   beamevt->DecodeFibers(name,beamevt->GetNFBMTriggers(name) - 1);//Decode the last one

   size_t i = beamevt->GetNFBMTriggers(name) - 1;
   std::cout << name << " has active fibers: ";
   for(size_t iF = 0; iF < beamevt->GetActiveFibers(name,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(name, i)[iF] << " "; 
   std::cout << std::endl;

   fbm = beam::FBM();
   fbm.ID = ip;
   name = fPairedStraightDevices[ip].second; 
   fbm.timeData[0] = dummyTriggerTimeLSB[name];
   fbm.timeData[1] = dummyTriggerTimeMSB[name];
   fbm.timeData[2] = dummyEventTimeLSB[name];
   fbm.timeData[3] = dummyEventTimeMSB[name];

   for(int i = 0; i < 6; ++i){
     fbm.fiberData[i] = dummyHit[name][i];
   }


   beamevt->AddFBMTrigger(name,fbm);
   beamevt->DecodeFibers(name,beamevt->GetNFBMTriggers(name) - 1);//Decode the last one

   i = beamevt->GetNFBMTriggers(name) - 1;
   std::cout << name << " has active fibers: ";
   for(size_t iF = 0; iF < beamevt->GetActiveFibers(name,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(name, i)[iF] << " "; 
   std::cout << std::endl;
 }
 
 std::cout << "Pairing" << std::endl;
 //GetPairedFBMInfo(*beamevt,1.50000e+08);
 GetPairedFBMInfo(*beamevt,fDummyEventTime);
 std::cout << "Paired" << std::endl;

 std::cout << "Unpaired" << std::endl;
 GetUnpairedFBMInfo(*beamevt,fDummyEventTime);
 std::cout << "Unpaired" << std::endl;

 std::cout << "Event stuff" << std::endl;
*/
 std::unique_ptr<std::vector<beam::ProtoDUNEBeamEvent> > beamData(new std::vector<beam::ProtoDUNEBeamEvent>);
 beamData->push_back(beam::ProtoDUNEBeamEvent(*beamevt));
 delete beamevt;
 std::cout << "Putting" << std::endl;
 e.put(std::move(beamData)/*,fOutputLabel*/);

 // Write out the to tree
 fOutTree->Fill();
 std::cout << "Put" << std::endl;
 
 if(usedEventTime) fMultipleTimes.clear();
}
// END BeamAna::produce
////////////////////////

void proto::BeamAna::InitXBPFInfo(beam::ProtoDUNEBeamEvent * beamevt){
  // Places a dummy trigger vector for each device

  // Make a vector the names of each of the devices being readout
  std::vector<std::string> monitors;
  size_t nDev = 0; 

  for(size_t id = 0; id < fDevices.size(); ++id){
    std::string name = fDevices[id];
    std::cout << fXBPFPrefix + fDevices[id] << std::endl;
    std::cout << "At: "      << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;
    std::cout << "Rotated: " << fRotations[name][0]   << " " << fRotations[name][1]   << " " << fRotations[name][2]   << std::endl;

    // Put the current device name on the list
    monitors.push_back(name);
    nDev++;
  }

  for(size_t id = 0; id < fPairedStraightDevices.size(); ++id){

    std::string name = fPairedStraightDevices[id].first;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    std::cout << nDev << std::endl;
    nDev++;
    monitors.push_back(name);

    name = fPairedStraightDevices[id].second;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;
    monitors.push_back(name);
    std::cout << nDev << std::endl;
    nDev++;
  }

  for(size_t id = 0; id < fPairedDevices.size(); ++id){

    std::string name = fPairedDevices[id].first;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    std::cout << nDev << std::endl;
    nDev++;
    monitors.push_back(name);

    name = fPairedDevices[id].second;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    std::cout << nDev << std::endl;
    nDev++;
    monitors.push_back(name);
  }

//  beamevt->InitFBMs(nDev);
  std::cout << "Initializing monitors";
  for(size_t id = 0; id < nDev; ++id){
    std::cout << " " << monitors.at(id);
  }
  std::cout << std::endl;
  beamevt->InitFBMs(monitors);
}
// END BeamAna::InitXBPFInfo
////////////////////////


void proto::BeamAna::parseXTOF(uint64_t time){
  std::cout << "Getting General trigger info " << std::endl;
  std::vector<double> coarseGeneralTrigger = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:coarse[]",fNRetries);
  std::vector<double> fracGeneralTrigger = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:frac[]",fNRetries); 
  std::cout << "Size of coarse,frac: " << coarseGeneralTrigger.size() << " " << fracGeneralTrigger.size() << std::endl; 

  std::cout << "Getting XBTF702 info " << std::endl;
  std::vector<double> coarseXBTF022702 = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/XBTF022702:coarse[]",fNRetries);
  std::vector<double> fracXBTF022702 = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/XBTF022702:frac[]",fNRetries); 
  std::cout << "Size of coarse,frac: " << coarseXBTF022702.size() << " " << fracXBTF022702.size() << std::endl; 
 
  std::cout << "Getting TOF1A info: " << fTOF1 << std::endl;
  std::vector<double> coarseTOF1A = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1A + ":coarse[]",fNRetries);
  std::vector<double> fracTOF1A =   FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1A + ":frac[]",fNRetries);
  std::cout << "Size of coarse,frac: " << coarseTOF1A.size() << " " << fracTOF1A.size() << std::endl; 

  std::cout << "Getting TOF1B info: " << fTOF1 << std::endl;
  std::vector<double> coarseTOF1B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1B + ":coarse[]",fNRetries);
  std::vector<double> fracTOF1B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1B + ":frac[]",fNRetries);
  std::cout << "Size of coarse,frac: " << coarseTOF1B.size() << " " << fracTOF1B.size() << std::endl; 

  std::cout << "Getting TOF2A info: " << fTOF2 << std::endl;
  std::vector<double> coarseTOF2A = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2A + ":coarse[]",fNRetries);
  std::vector<double> fracTOF2A = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2A + ":frac[]",fNRetries);
  std::cout << "Size of coarse,frac: " << coarseTOF2A.size() << " " << fracTOF2A.size() << std::endl; 

  std::cout << "Getting TOF2B info: " << fTOF2 << std::endl;
  std::vector<double> coarseTOF2B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2B + ":coarse[]",fNRetries);
  std::vector<double> fracTOF2B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2B + ":frac[]",fNRetries);
  std::cout << "Size of coarse,frac: " << coarseTOF2B.size() << " " << fracTOF2B.size() << std::endl; 


  std::vector<double> unorderedGenTrigTime;
  std::vector<double> unorderedXBTF022702Time;
  std::vector<double> unorderedTOF1ATime;
  std::vector<double> unorderedTOF1BTime;
  std::vector<double> unorderedTOF2ATime;
  std::vector<double> unorderedTOF2BTime;
    
  for(size_t i = 0; i < 4000; ++i){
    std::cout << i << " "  << 8.*coarseGeneralTrigger[i] + fracGeneralTrigger[i]/512. << std::endl;
    fGenTrigCoarse = coarseGeneralTrigger[i];
    fGenTrigFrac = fracGeneralTrigger[i];
    unorderedGenTrigTime.push_back(fGenTrigCoarse*8. + fGenTrigFrac/512.);

    if (fGenTrigFrac == 0.0) break;
    fGenTrigTree->Fill();

    std::cout << i << " "  << 8.*coarseTOF1A[i] + fracTOF1A[i]/512. << std::endl;
    fXTOF1ACoarse = coarseTOF1A[i];
    fXTOF1AFrac = fracTOF1A[i];
    unorderedTOF1ATime.push_back(fXTOF1ACoarse*8. + fXTOF1AFrac/512.);
    fXTOF1ATree->Fill();
    
    std::cout << i << " "  << 8.*coarseTOF1B[i] + fracTOF1B[i]/512. << std::endl;
    fXTOF1BCoarse = coarseTOF1B[i];
    fXTOF1BFrac = fracTOF1B[i];
    unorderedTOF1BTime.push_back(fXTOF1BCoarse*8. + fXTOF1BFrac/512.);
    fXTOF1BTree->Fill();

    std::cout << i << " "  << 8.*coarseTOF2A[i] + fracTOF2A[i]/512. << std::endl;
    fXTOF2ACoarse = coarseTOF2A[i];
    fXTOF2AFrac = fracTOF2A[i];
    unorderedTOF2ATime.push_back(fXTOF2ACoarse*8. + fXTOF2AFrac/512.);
    fXTOF2ATree->Fill();
    
    std::cout << i << " "  << 8.*coarseTOF2B[i] + fracTOF2B[i]/512. << std::endl;
    fXTOF2BCoarse = coarseTOF2B[i];
    fXTOF2BFrac = fracTOF2B[i];
    unorderedTOF2BTime.push_back(fXTOF2BCoarse*8. + fXTOF2BFrac/512.);
    fXTOF2BTree->Fill();

    std::cout << i << " "  << 8.*coarseTOF2B[i] + fracTOF2B[i]/512. << std::endl;
    fXBTF022702Coarse = coarseXBTF022702[i];
    fXBTF022702Frac =     fracXBTF022702[i];
    unorderedXBTF022702Time.push_back(fXTOF2BCoarse*8. + fXTOF2BFrac/512.);
    fXBTF022702Tree->Fill();
  }

  //Go through the unordered TOF triggers
  //Look for coincidences between TOF1 and TOF2
  //There should only be one match between A and B
  for(size_t iT = 0; iT < unorderedGenTrigTime.size(); ++iT){
    bool found_TOF1 = false;
    bool found_TOF2 = false;

    double the_time = unorderedGenTrigTime[iT];
    double TOF1A_time;
    double TOF1B_time;
    double TOF2A_time;
    double TOF2B_time;

    double the_TOF1 = -1.;
    double the_TOF2 = -1.;

    //Technically out of bounds of the vectors, but it simplifies things
    for(size_t iT2 = 0; iT2 < unorderedGenTrigTime.size(); ++iT2){
      if (iT2 < unorderedTOF1ATime.size()){
        TOF1A_time = unorderedTOF1ATime[iT2];
        if(fabs(the_time - TOF1A_time) < 10000. ){
          found_TOF1 = true; 
          the_TOF1 = TOF1A_time;
        }
      }
      if (iT2 < unorderedTOF1BTime.size()){
        TOF1B_time = unorderedTOF1BTime[iT2];
        if(!found_TOF1 && (fabs(the_time - TOF1B_time) < 10000. )){
          found_TOF1 = true; 
          the_TOF1 = TOF1B_time;
        }
      }

      if (iT2 < unorderedTOF2ATime.size()){
        TOF2A_time = unorderedTOF2ATime[iT2];
        if(fabs(the_time - TOF2A_time) < 10000. ){
          found_TOF2 = true; 
          the_TOF2 = TOF2A_time;
        }
      }
      if (iT2 < unorderedTOF2BTime.size()){
        TOF2B_time = unorderedTOF2BTime[iT2];
        if(!found_TOF2 && (fabs(the_time - TOF2B_time) < 10000. )){
          found_TOF2 = true; 
          the_TOF2 = TOF2B_time;
        }
      }

      if(found_TOF1 && found_TOF2){
        beamevt->AddT0(the_time);
        beamevt->AddTOF0Trigger(the_TOF1);
        beamevt->AddTOF1Trigger(the_TOF2);
        break;
      }

    }
  }

}
// END BeamAna::parseXTOF
////////////////////////


void proto::BeamAna::parseXCET(uint64_t time){
  if(fCKov1 != ""){  
    std::cout << "Getting CKov1 info: " << fCKov1 << std::endl;

    std::vector< double > countsTrigCKov1 = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":countsTrig",fNRetries);
    std::vector< double > countsCKov1     = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":counts",fNRetries);
    std::vector< double > pressureCKov1   = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":pressure",fNRetries);
    std::cout << "countsTrig: " << countsTrigCKov1[0] << std::endl; 
    std::cout << "counts: "     << countsCKov1[0] << std::endl;
    std::cout << "pressure: "   << pressureCKov1[0] << std::endl;
  }
  
  if(fCKov2 != ""){
    std::cout << "Getting CKov2 info: " << fCKov2 << std::endl;

    std::vector< double > countsTrigCKov2 = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":countsTrig",fNRetries);
    std::vector< double > countsCKov2     = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":counts",fNRetries);
    std::vector< double > pressureCKov2   = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":pressure",fNRetries);
    std::cout << "countsTrig: " << countsTrigCKov2[0] << std::endl; 
    std::cout << "counts: "     << countsCKov2[0] << std::endl;
    std::cout << "pressure: "   << pressureCKov2[0] << std::endl;
  }
}
// END BeamAna::parseXCET
////////////////////////



////////////////////////
// 
void proto::BeamAna::parseGeneralXBPF(std::string name, uint64_t time, size_t ID){

  // Retrieve the number of counts in the BPF
  std::vector<double> counts;
  counts = FetchWithRetries< std::vector<double> >(time, fXBPFPrefix + name + ":countsRecords[]",fNRetries);
  std::cout << "Counts: " << counts.size() << std::endl;
  for(size_t i = 0; i < counts.size(); ++i){
    std::cout << counts[i] << std::endl;
  }

  std::vector<double> data;
  data = FetchWithRetries< std::vector<double> >(time, fXBPFPrefix + name + ":eventsData[]",fNRetries);
  std::cout << "Data: " << data.size() << std::endl;

  // If the number of counts is larger than the data we have
  // we bail (without an error???)
  if(counts[1] > data.size()){
    return;
  }
  
  beam::FBM fbm;
  fbm.ID = ID;
    
  //Use this just in case any are out of sync?
  //Shouldn't be, but just to be safe...
  //Helps cut down on time
  std::vector<size_t> leftOvers;      
  for(size_t lo = 0; lo < beamevt->GetNT0(); ++lo){
    leftOvers.push_back(lo);
  }
 
  //Skipping anything > 500, the data seems to be corrupted now
  for(size_t i = 0; (i < counts[1] && i < 500); ++i){      
      std::cout << "Count: " << i << std::endl;
    
    for(int j = 0; j < 10; ++j){
      double theData = data[20*i + (2*j + 1)];
      std::cout << std::setw(15) << theData ;
      if(j < 4){
	fbm.timeData[j] = theData;           
      }else{
	fbm.fiberData[j - 4] = theData;
      } // endif j
    } // endfor j

    // Check the time data for corruption
    if(fbm.timeData[1] < .0000001){
      std::cout << "Skipping bad time" << std::endl;
      continue;
    } // endif timeData
    fbm.timeStamp = fbm.timeData[0]*8.;  // Timestamp is 8x the timeData value
    
    //Go through the valid Good Particles, and emplace the FBM 
    std::cout << "Checking " << beamevt->GetNT0() << " triggers " << leftOvers.size() << std::endl;

    for(std::vector<size_t>::iterator ip = leftOvers.begin(); ip != leftOvers.end(); ++ip){
      // std::cout << "\t" << fbm.timeStamp << " " << beamevt->GetT0(*ip) << std::endl;

      // Compute the time delta between the timeStamp and the T0, see if it's less than 5000
      if( fabs(fbm.timeStamp - beamevt->GetT0(*ip)) < 5000.){
	if(beamevt->GetFBM(name, *ip).ID != -1){
	  std::cout << "Warning: Replacing non-dummy FBM at "
		    << name << " " << *ip << std::endl;
	} // endif GetFBM
	
	 std::cout << "Replacing at timestamp " << fbm.timeStamp << std::endl;
	beamevt->ReplaceFBMTrigger(name, fbm, *ip);
	leftOvers.erase(ip);
	break;
      } // endif time delta calc 
    } // endfor ip

    std::cout << std::endl;
  } // endfor i

  for(size_t i = 0; i < beamevt->GetNFBMTriggers(name); ++i){
    beamevt->DecodeFibers(name,i);
    std::cout << name << " at time: "
	      << beamevt->DecodeFiberTime(name, i) << " has active fibers: ";
    for(size_t iFiber = 0; iFiber < beamevt->GetActiveFibers(name,i).size(); ++iFiber){
      std::cout << beamevt->GetActiveFibers(name, i)[iFiber] << " ";
    }
    std::cout << std::endl;
    
    *fActiveFibers[name] = beamevt->GetActiveFibers(name,i);
    fProfTime[name] = beamevt->DecodeFiberTime(name, i);
    std::cout << beamevt->ReturnTriggerAndTime(name,i)[0] << " "
	      << beamevt->ReturnTriggerAndTime(name,i)[1] << " "
	      << beamevt->ReturnTriggerAndTime(name,i)[2] << " "
	      << beamevt->ReturnTriggerAndTime(name,i)[3] << std::endl;
    fProfTrigger1[name] = beamevt->ReturnTriggerAndTime(name,i)[0];
    fProfTrigger2[name] = beamevt->ReturnTriggerAndTime(name,i)[1];
    fProfTime1[name] = beamevt->ReturnTriggerAndTime(name,i)[2];
    fProfTime2[name] = beamevt->ReturnTriggerAndTime(name,i)[3];
    fProfTree[name]->Fill(); 
    
  } // endfor i

}
// END BeamAna::parseGeneralXBFP
////////////////////////

////////////////////////
// 
void proto::BeamAna::parseXBPF(uint64_t time){
  for(size_t d = 0; d < fDevices.size(); ++d){
    std::string name = fDevices[d];
    std::cout <<"Device: " << name << std::endl;
    parseGeneralXBPF(name, time, d);
  }  
}
// END BeamAna::parseXBFP
////////////////////////

////////////////////////
// 
void proto::BeamAna::parsePairedXBPF(uint64_t time){
  for(size_t d = 0; d < fPairedDevices.size(); ++d){
    std::string name = fPairedDevices[d].first;
    std::cout <<"Device: " << name << std::endl;
    parseGeneralXBPF(name, time, d);

    name = fPairedDevices[d].second;
    std::cout <<"Device: " << name << std::endl;
    parseGeneralXBPF(name, time, d);
  }
}
// END BeamAna::parsePairedXBFP
////////////////////////

////////////////////////
// 
void proto::BeamAna::parsePairedStraightXBPF(uint64_t time){
  for(size_t d = 0; d < fPairedStraightDevices.size(); ++d){
    std::string name = fPairedStraightDevices[d].first;
    std::cout <<"Device: " << name << std::endl;
    parseGeneralXBPF(name, time, d);

    name = fPairedStraightDevices[d].second;
    std::cout <<"Device: " << name << std::endl;
    parseGeneralXBPF(name, time, d);
  }
}
// END BeamAna::parsePairedStrightXBFP
////////////////////////

void proto::BeamAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  

  fTOFHist = tfs->make<TH1F>("TOF","",100,0,100);
  fCKovHist = tfs->make<TH1F>("CKov","",4,0,4);

  fOutTree = tfs->make<TTree>("tree", "lines"); 
  theTrack = new recob::Track();
  eventTime = 0.;
  fOutTree->Branch("Track", &theTrack);
  fOutTree->Branch("Time", &eventTime);

  fTrackTree = tfs->make<TTree>("tracks","");
  trackX = new std::vector<double>;
  trackY = new std::vector<double>;
  trackZ = new std::vector<double>;
  fTrackTree->Branch("X", &trackX);
  fTrackTree->Branch("Y", &trackY);
  fTrackTree->Branch("Z", &trackZ);

  fMatchedTriggers = tfs->make<TTree>("matched","");
  
  fMatchedTriggers->Branch("Gen", &matchedGen);
  fMatchedTriggers->Branch("TOF1", &matchedTOF1);
  fMatchedTriggers->Branch("TOF2", &matchedTOF2);

  for(size_t i = 0; i < fPairedDevices.size(); ++i){
    std::string name = "BeamProf2D_" + fPairedDevices[i].first + "_" + fPairedDevices[i].second;
    std::string title = fPairedDevices[i].first + ", " + fPairedDevices[i].second;
    fBeamProf2D.push_back( tfs->make<TH2F>(name.c_str(),title.c_str(),192,0,192,192,0,192) );

    name = "BeamProf1D_" + fPairedDevices[i].first;
    title = fPairedDevices[i].first;
    fBeamProf1D[fPairedDevices[i].first] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );

    name = "BeamProf1D_" + fPairedDevices[i].second;
    title = fPairedDevices[i].second;
    fBeamProf1D[fPairedDevices[i].second] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );
    
    name = "Fibers_" + fPairedDevices[i].first;
    fActiveFibers[fPairedDevices[i].first] = new std::vector<short>;

    fProfTime[fPairedDevices[i].first] = 0.;
    fProfTrigger1[fPairedDevices[i].first] = 0.;
    fProfTrigger2[fPairedDevices[i].first] = 0.;
    fProfTime1[fPairedDevices[i].first] = 0.;
    fProfTime2[fPairedDevices[i].first] = 0.;

    matchedXBPF[fPairedDevices[i].first] = 0.;
    fMatchedTriggers->Branch((fPairedDevices[i].first).c_str(), &matchedXBPF[fPairedDevices[i].first]);


    fProfTree[fPairedDevices[i].first] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedDevices[i].first]->Branch("time", &fProfTime[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("fibers", &fActiveFibers[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("trigger_1", &fProfTrigger1[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("trigger_2", &fProfTrigger2[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("time_1", &fProfTime1[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("time_2", &fProfTime2[fPairedDevices[i].first]);

    name = "Fibers_" + fPairedDevices[i].second;
    fActiveFibers[fPairedDevices[i].second] = new std::vector<short>;

    fProfTime[fPairedDevices[i].second] = 0.;
    fProfTrigger1[fPairedDevices[i].second] = 0.;
    fProfTrigger2[fPairedDevices[i].second] = 0.;
    fProfTime1[fPairedDevices[i].second] = 0.;
    fProfTime2[fPairedDevices[i].second] = 0.;

    matchedXBPF[fPairedDevices[i].second] = 0.;
    fMatchedTriggers->Branch((fPairedDevices[i].second).c_str(), &matchedXBPF[fPairedDevices[i].second]);

    fProfTree[fPairedDevices[i].second] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedDevices[i].second]->Branch("time", &fProfTime[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("fibers", &fActiveFibers[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("trigger_1", &fProfTrigger1[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("trigger_2", &fProfTrigger2[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("time_1", &fProfTime1[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("time_2", &fProfTime2[fPairedDevices[i].second]);
  }

  for(size_t i = 0; i < fPairedStraightDevices.size(); ++i){
    std::string name = "BeamProf2D_" + fPairedStraightDevices[i].first + "_" + fPairedStraightDevices[i].second;
    std::string title = fPairedStraightDevices[i].first + ", " + fPairedStraightDevices[i].second;
    fBeamProf2D.push_back( tfs->make<TH2F>(name.c_str(),title.c_str(),192,0,192,192,0,192) );

    name = "BeamProf1D_" + fPairedStraightDevices[i].first;
    title = fPairedStraightDevices[i].first;
    fBeamProf1D[fPairedStraightDevices[i].first] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );

    name = "BeamProf1D_" + fPairedStraightDevices[i].second;
    title = fPairedStraightDevices[i].second;
    fBeamProf1D[fPairedStraightDevices[i].second] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );
    
    name = "Fibers_" + fPairedStraightDevices[i].first;
    fActiveFibers[fPairedStraightDevices[i].first] = new std::vector<short>;

    fProfTime[fPairedStraightDevices[i].first] = 0.;
    fProfTrigger1[fPairedStraightDevices[i].first] = 0.;
    fProfTrigger2[fPairedStraightDevices[i].first] = 0.;
    fProfTime1[fPairedStraightDevices[i].first] = 0.;
    fProfTime2[fPairedStraightDevices[i].first] = 0.;

    matchedXBPF[fPairedStraightDevices[i].first] = 0.;
    fMatchedTriggers->Branch((fPairedStraightDevices[i].first).c_str(), &matchedXBPF[fPairedStraightDevices[i].first]);

    fProfTree[fPairedStraightDevices[i].first] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedStraightDevices[i].first]->Branch("time", &fProfTime[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("fibers", &fActiveFibers[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("trigger_1", &fProfTrigger1[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("trigger_2", &fProfTrigger2[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("time_1", &fProfTime1[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("time_2", &fProfTime2[fPairedStraightDevices[i].first]);

    name = "Fibers_" + fPairedStraightDevices[i].second;
    fActiveFibers[fPairedStraightDevices[i].second] = new std::vector<short>;

    fProfTime[fPairedStraightDevices[i].second] = 0.;
    fProfTrigger1[fPairedStraightDevices[i].second] = 0.;
    fProfTrigger2[fPairedStraightDevices[i].second] = 0.;
    fProfTime1[fPairedStraightDevices[i].second] = 0.;
    fProfTime2[fPairedStraightDevices[i].second] = 0.;

    matchedXBPF[fPairedStraightDevices[i].second] = 0.;
    fMatchedTriggers->Branch((fPairedStraightDevices[i].second).c_str(), &matchedXBPF[fPairedStraightDevices[i].second]);

    fProfTree[fPairedStraightDevices[i].second] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedStraightDevices[i].second]->Branch("time", &fProfTime[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("fibers", &fActiveFibers[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("trigger_1", &fProfTrigger1[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("trigger_2", &fProfTrigger2[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("time_1", &fProfTime1[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("time_2", &fProfTime2[fPairedStraightDevices[i].second]);
  }

  for(size_t i = 0; i < fDevices.size(); ++i){
    std::string name = "BeamProf1D_" + fDevices[i];
    std::string title = fDevices[i];
    fBeamProf1D[fDevices[i]] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );

    name = "Fibers_" + fDevices[i];
    fProfTree[fDevices[i]] = ( tfs->make<TTree>(name.c_str(), "XBPF") );

    fProfTrigger1[fDevices[i]] = 0.;
    fProfTrigger2[fDevices[i]] = 0.;
    fProfTime1[fDevices[i]] = 0.;
    fProfTime2[fDevices[i]] = 0.;
    fProfTime[fDevices[i]] = 0.;

    matchedXBPF[fDevices[i]] = 0.;
    fMatchedTriggers->Branch((fDevices[i]).c_str(), &matchedXBPF[fDevices[i]]);

    fActiveFibers[fDevices[i]] = new std::vector<short>;
    fProfTree[fDevices[i]]->Branch("time", &fProfTime[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("fibers", &fActiveFibers[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("trigger_1", &fProfTrigger1[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("trigger_2", &fProfTrigger2[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("time_1", &fProfTime1[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("time_2", &fProfTime2[fDevices[i]]);
  }

  fGenTrigTree = tfs->make<TTree>("GenTrig","");
  fGenTrigTree->Branch("coarse", &fGenTrigCoarse);
  fGenTrigTree->Branch("frac", &fGenTrigFrac);

  fXTOF1ATree = tfs->make<TTree>("TOF1A","");
  fXTOF1ATree->Branch("coarse", &fXTOF1ACoarse);
  fXTOF1ATree->Branch("frac", &fXTOF1AFrac);

  fXTOF1BTree = tfs->make<TTree>("TOF1B","");
  fXTOF1BTree->Branch("coarse", &fXTOF1BCoarse);
  fXTOF1BTree->Branch("frac", &fXTOF1BFrac);

  fXTOF2ATree = tfs->make<TTree>("TOF2A","");
  fXTOF2ATree->Branch("coarse", &fXTOF2ACoarse);
  fXTOF2ATree->Branch("frac", &fXTOF2AFrac);

  fXTOF2BTree = tfs->make<TTree>("TOF2B","");
  fXTOF2BTree->Branch("coarse", &fXTOF2BCoarse);
  fXTOF2BTree->Branch("frac", &fXTOF2BFrac);

  fXBTF022702Tree = tfs->make<TTree>("XBTF022702","");
  fXBTF022702Tree->Branch("coarse", &fXBTF022702Coarse);
  fXBTF022702Tree->Branch("frac", &fXBTF022702Frac);
}

void proto::BeamAna::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endJob()
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fBundleName  = p.get<std::string>("BundleName");
  fOutputLabel = p.get<std::string>("OutputLabel");
  fInputLabel  = p.get<std::string>("InputLabel");
  fURLStr      = p.get<std::string>("URLStr");
  fNRetries    = p.get<int>("NRetries");
  fValidWindow = p.get<double>("ValidWindow");
  fTimeWindow  = p.get<double>("TimeWindow");
  fFixedTime   = p.get<uint64_t>("FixedTime");
  fMultipleTimes = p.get< std::vector<uint64_t> >("MultipleTimes");
  fTolerance   = p.get<double>("Tolerance");

  fDevices     = p.get< std::vector< std::string > >("Devices");    
  fPairedDevices = p.get< std::vector< std::pair<std::string, std::string> > >("PairedDevices");
  fPairedStraightDevices = p.get< std::vector< std::pair<std::string, std::string> > >("PairedStraightDevices");   

  std::vector< std::pair<std::string, std::string> >  tempTypes = p.get<std::vector< std::pair<std::string, std::string> >>("DeviceTypes");
  fDeviceTypes  = std::map<std::string, std::string>(tempTypes.begin(), tempTypes.end() );

  //Location of Device 
  std::vector< std::pair<std::string, std::array<double,3> > > tempCoords = p.get<std::vector< std::pair<std::string, std::array<double,3> > > >("Coordinates");
  fCoordinates = std::map<std::string, std::array<double,3> >(tempCoords.begin(), tempCoords.end());
//  fCoordinates = p.get< std::vector< std::array<double,3> > >("Coordinates");

  //Rotation of Device
  std::vector< std::pair<std::string, std::array<double,3> > > tempRots = p.get< std::vector< std::pair<std::string, std::array<double,3> > > >("Rotations"); 
  fRotations = std::map<std::string, std::array<double,3> >(tempRots.begin(), tempRots.end());
//  fRotations   = p.get< std::vector< std::array<double,3> > >("Rotations"); 

  //Deminsion of Fibers
  std::vector< std::pair<std::string, double> > tempFiberDims = p.get< std::vector<std::pair<std::string,double> > >("Dimension");
  fFiberDimension = std::map<std::string, double>(tempFiberDims.begin(), tempFiberDims.end());
//  fFiberDimension = p.get< std::vector<double> >("Dimension");

  std::array<double,3> detCoords = p.get<std::array<double,3>>("GlobalDetCoords");
  fGlobalDetCoords = TVector3(detCoords[0],detCoords[1],detCoords[2]);

  fDetRotation = p.get<std::array<double,3>>("DetRotation");
  
  //XTOF devices 
  fTOF1 = p.get< std::string >("TOF1");
  fTOF2 = p.get< std::string >("TOF2");

  fTOF1A = fTOF1 + "A";
  fTOF1B = fTOF1 + "B";
  fTOF2A = fTOF2 + "A";
  fTOF2B = fTOF2 + "B";

  fCKov1 = p.get< std::string >("CKov1");
  fCKov2 = p.get< std::string >("CKov2");


  fXBPFPrefix      = p.get<std::string>("XBPFPrefix");
  fXTOFPrefix      = p.get<std::string>("XTOFPrefix");
  fXCETPrefix      = p.get<std::string>("XCETPrefix");
  fDummyEventTime = p.get<double>("DummyEventTime");


  //New parameters to match Leigh's
  fRotateMonitorXZ = p.get<double>("RotateMonitorXZ");
  fRotateMonitorYZ = p.get<double>("RotateMonitorYZ");

  fFirstTrackingProfZ  = p.get<double>("FirstTrackingProfZ");
  fSecondTrackingProfZ = p.get<double>("SecondTrackingProfZ");
  fNP04FrontZ          = p.get<double>("NP04FrontZ");  

  fBeamX               = p.get<double>("BeamX");
  fBeamY               = p.get<double>("BeamY");
  fBeamZ               = p.get<double>("BeamZ");
}

std::bitset<sizeof(long)*CHAR_BIT> proto::BeamAna::toBinary(long num){
   
  std::bitset<sizeof(double)*CHAR_BIT> mybits(num);  
  std::bitset<32> upper, lower;
  for(int i = 0; i < 32; ++i){
    lower[i] = mybits[i];
    upper[i] = mybits[i + 32];   
  }
  if(upper.any()) std::cout << "WARNING: NONZERO HALF" << std::endl;

  return mybits;
}

TVector3 proto::BeamAna::ConvertProfCoordinates(double x, double y, double z, double zOffset){
  double off = fNP04FrontZ - zOffset;

  TVector3 old(x,y,z);

  double newX = x*fBMBasisX.X() + y*fBMBasisY.X() + /*(z-zOffset)*fBMBasisZ.X()*/ + off*fabs(fBMBasisZ.X());
  double newY = x*fBMBasisX.Y() + y*fBMBasisY.Y() + /*(z-zOffset)*fBMBasisZ.Y()*/ + off*fabs(fBMBasisZ.Y());
  double newZ = x*fBMBasisX.Z() + y*fBMBasisY.Z() + /*(z-zOffset)              */ - off*fabs(fBMBasisZ.Z());

  newX += fBeamX*10.;
  newY += fBeamY*10.;
  newZ += fBeamZ*10.;

  TVector3 result(newX/10., newY/10., newZ/10.);
  return result;
}

void proto::BeamAna::BeamMonitorBasisVectors(){
  std::cout << "Rotating" << std::endl;
  fBMBasisX = TVector3(1.,0.,0.);
  fBMBasisY = TVector3(0.,1.,0.);
  fBMBasisZ = TVector3(0.,0.,1.);
  RotateMonitorVector(fBMBasisX);
  std::cout << fBMBasisX.X() << " " << fBMBasisX.Y() << " " << fBMBasisX.Z() << std::endl;
  RotateMonitorVector(fBMBasisY);
  std::cout << fBMBasisY.X() << " " << fBMBasisY.Y() << " " << fBMBasisY.Z() << std::endl;
  RotateMonitorVector(fBMBasisZ);
  std::cout << fBMBasisZ.X() << " " << fBMBasisZ.Y() << " " << fBMBasisZ.Z() << std::endl;
}

void proto::BeamAna::RotateMonitorVector(TVector3 &vec){
  vec.RotateY(fRotateMonitorXZ * TMath::Pi()/180.);
  vec.RotateX(fRotateMonitorYZ * TMath::Pi()/180.);
}

void proto::BeamAna::MakeTrack(size_t theTrigger){
  
  std::cout << "Making Track for time: " << beamevt->GetT0(theTrigger) << std::endl;

  std::string firstUpstreamName = fPairedStraightDevices[0].first;
  std::string secondUpstreamName = fPairedStraightDevices[0].second;
  std::string firstDownstreamName = fPairedStraightDevices[1].first;
  std::string secondDownstreamName = fPairedStraightDevices[1].second;


  //Get the active fibers from the upstream tracking XBPF
  std::vector<short> firstUpstreamFibers = beamevt->GetActiveFibers(firstUpstreamName, theTrigger);
  std::vector<short> secondUpstreamFibers = beamevt->GetActiveFibers(secondUpstreamName, theTrigger);

  std::cout << firstUpstreamName << " has " << firstUpstreamFibers.size() << " active fibers at time " << beamevt->GetFiberTime(firstUpstreamName,theTrigger) << std::endl;
  for(size_t i = 0; i < firstUpstreamFibers.size(); ++i){
    std::cout << firstUpstreamFibers[i] << " ";
  }
  std::cout << std::endl;

  std::cout << secondUpstreamName << " has " << secondUpstreamFibers.size() << " active fibers at time " << beamevt->GetFiberTime(secondUpstreamName,theTrigger) << std::endl;
  for(size_t i = 0; i < secondUpstreamFibers.size(); ++i){
    std::cout << secondUpstreamFibers[i] << " ";
  }
  std::cout << std::endl;
  //////////////////////////////////////////////

  //Get the active fibers from the downstream tracking XBPF
  std::vector<short> firstDownstreamFibers = beamevt->GetActiveFibers(firstDownstreamName, theTrigger);
  std::vector<short> secondDownstreamFibers = beamevt->GetActiveFibers(secondDownstreamName, theTrigger);

  std::cout << firstDownstreamName << " has " << firstDownstreamFibers.size() << " active fibers at time " << beamevt->GetFiberTime(firstDownstreamName,theTrigger) << std::endl;
  for(size_t i = 0; i < firstDownstreamFibers.size(); ++i){
    std::cout << firstDownstreamFibers[i] << " ";
  }
  std::cout << std::endl;

  std::cout << secondDownstreamName << " has " << secondDownstreamFibers.size() << " active fibers at time " << beamevt->GetFiberTime(secondDownstreamName,theTrigger) << std::endl;
  for(size_t i = 0; i < secondDownstreamFibers.size(); ++i){
    std::cout << secondDownstreamFibers[i] << " ";
  }
  std::cout << std::endl;
  //////////////////////////////////////////////

  if( (firstUpstreamFibers.size() < 1) || (secondUpstreamFibers.size() < 1) || (firstDownstreamFibers.size() < 1) || (secondDownstreamFibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not making track" << std::endl;
    return;
  }
  else if( (firstUpstreamFibers.size() > 5) || (secondUpstreamFibers.size() > 5) || (firstDownstreamFibers.size() > 5) || (secondDownstreamFibers.size() > 5) ){
    std::cout << "Warning, too many (>5) active fibers in at least one Beam Profiler. Not making track" << std::endl;
    return;
  }

  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  std::vector< std::pair<size_t, size_t> > upstreamPairedFibers;
  std::vector< std::pair<size_t, size_t> > downstreamPairedFibers;
  std::string firstUpstreamType = fDeviceTypes[firstUpstreamName]; 
  std::string secondUpstreamType = fDeviceTypes[secondUpstreamName]; 
  std::string firstDownstreamType = fDeviceTypes[firstDownstreamName]; 
  std::string secondDownstreamType = fDeviceTypes[secondDownstreamName]; 

  std::vector< TVector3 > upstreamPositions;
  std::vector< TVector3 > downstreamPositions; 
  


  //Pair the upstream fibers together
  std::cout << "Upstream" << std::endl;
  for(size_t iF1 = 0; iF1 < firstUpstreamFibers.size(); ++iF1){
    
    size_t firstFiber = firstUpstreamFibers[iF1];

    for(size_t iF2 = 0; iF2 < secondUpstreamFibers.size(); ++iF2){
      size_t secondFiber = secondUpstreamFibers[iF2];

      std::cout << "Paired: " << firstFiber << " " << secondFiber << std::endl; 
      upstreamPairedFibers.push_back(std::make_pair(firstFiber, secondFiber));

      if (iF2 < secondUpstreamFibers.size() - 1){
        if (secondUpstreamFibers[iF2] == (secondUpstreamFibers[iF2 + 1] - 1)) ++iF2;
      }
    }

    if (iF1 < firstUpstreamFibers.size() - 1){
      if (firstUpstreamFibers[iF1] == (firstUpstreamFibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  for(size_t iF = 0; iF < upstreamPairedFibers.size(); ++iF){
    
    std::pair<size_t,size_t> thePair = upstreamPairedFibers.at(iF);
 
    if(firstUpstreamType == "horiz" && secondUpstreamType == "vert"){
      double xPos = GetPosition(firstUpstreamName, thePair.first);
      double yPos = GetPosition(secondUpstreamName, thePair.second);
      
      std::cout << "normal " << xPos << " " << yPos <<  std::endl;
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
      std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
      upstreamPositions.push_back( posInDet );
    }
    else if(firstUpstreamType == "vert" && secondUpstreamType == "horiz"){
      double yPos = GetPosition(firstUpstreamName, thePair.first);
      double xPos = GetPosition(secondUpstreamName, thePair.second);
      std::cout << "normal " << xPos << " " << yPos <<  std::endl;

      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
      std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
      upstreamPositions.push_back( posInDet );
    }

  }
  
  std::cout << "Downstream" << std::endl;
  //Pair the downstream fibers together
  for(size_t iF1 = 0; iF1 < firstDownstreamFibers.size(); ++iF1){
    
    size_t firstFiber = firstDownstreamFibers[iF1];

    for(size_t iF2 = 0; iF2 < secondDownstreamFibers.size(); ++iF2){
      size_t secondFiber = secondDownstreamFibers[iF2];

      std::cout << "Paired: " << firstFiber << " " << secondFiber << std::endl; 
      downstreamPairedFibers.push_back(std::make_pair(firstFiber, secondFiber));

      if (iF2 < secondDownstreamFibers.size() - 1){
        if (secondDownstreamFibers[iF2] == (secondDownstreamFibers[iF2 + 1] - 1)) ++iF2;
      }
    }

    if (iF1 < firstDownstreamFibers.size() - 1){
      if (firstDownstreamFibers[iF1] == (firstDownstreamFibers[iF1 + 1] - 1)) ++iF1;
    }
  }

  for(size_t iF = 0; iF < downstreamPairedFibers.size(); ++iF){
    
    std::pair<size_t,size_t> thePair = downstreamPairedFibers.at(iF);
 
    if(firstDownstreamType == "horiz" && secondDownstreamType == "vert"){
      double xPos = GetPosition(firstDownstreamName, thePair.first);
      double yPos = GetPosition(secondDownstreamName, thePair.second);

      std::cout << "normal " << xPos << " " << yPos <<  std::endl;
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
      std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
      downstreamPositions.push_back( posInDet );
    }
    else if(firstDownstreamType == "vert" && secondDownstreamType == "horiz"){
      double yPos = GetPosition(firstDownstreamName, thePair.first);
      double xPos = GetPosition(secondDownstreamName, thePair.second);

      std::cout << "normal " << xPos << " " << yPos <<  std::endl;
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
      std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
      downstreamPositions.push_back( posInDet );
    }

  }
 
  //Just for creating tracks
  std::vector< std::vector<double> > dummy;
  std::vector<double> mom(3, util::kBogusD);
  /// 

  for(size_t iU = 0; iU < upstreamPositions.size(); ++iU){
    for(size_t iD = 0; iD < downstreamPositions.size(); ++iD){
      std::vector<TVector3> thePoints;
      thePoints.push_back(upstreamPositions.at(iU));
      thePoints.push_back(downstreamPositions.at(iD));

      //Now project the last point to the TPC face
      thePoints.push_back( ProjectToTPC(thePoints[0],thePoints[1]) );    

      
      std::vector<TVector3> theMomenta;
      //Just push back the unit vector for each point 
      //For now.
      //Eventually, use momentum from curvature?
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );

      recob::Track * tempTrack = new recob::Track(thePoints, theMomenta, dummy, mom, 1);      
      theTracks.push_back(tempTrack);


    }
  }

}

//Gets info from FBMs matching in time
void proto::BeamAna::GetPairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time){

  //This method goes through the paired FBMs and gets 2D hits in their profile that match in time
  //to the input time and fills 2D histograms. 
  //If there are any missing triggers in a pair, then it skips filling the respective hist
  //
  //Need to figure out what to do with multiple active fibers

  std::map<std::string, std::vector<size_t> > goodTriggers;
  for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
    std::string name = fPairedDevices[ip].first;
    std::vector<size_t> triggers;
//    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
//        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;

    name = fPairedDevices[ip].second;
    triggers.clear();
//    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
//        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;
  }
   
  for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
    std::string nameOne = fPairedDevices[ip].first;
    std::string nameTwo = fPairedDevices[ip].second;

    size_t triggerSizeOne = goodTriggers[nameOne].size();
    size_t triggerSizeTwo = goodTriggers[nameTwo].size();

    if(triggerSizeOne < 1){
//      std::cout << "Missing trigger for " << nameOne << std::endl;
    }
    if(triggerSizeTwo < 1){
//      std::cout << "Missing trigger for " << nameTwo << std::endl;
    }

    if(triggerSizeOne > 0 && triggerSizeTwo > 0){
      //std::cout << "Found good triggers for " << nameOne << " " << nameTwo << std::endl;
      size_t triggerOne = goodTriggers[nameOne][0];
      size_t triggerTwo = goodTriggers[nameTwo][0];
      fBeamProf2D[ip]->Fill(beamevt.GetActiveFibers(nameOne,triggerOne)[0], beamevt.GetActiveFibers(nameTwo,triggerTwo)[0]);
    }
  }

}

void proto::BeamAna::GetPairedStraightFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time){

  //This method goes through the paired FBMs and gets 2D hits in their profile that match in time
  //to the input time and fills 2D histograms. 
  //If there are any missing triggers in a pair, then it skips filling the respective hist
  //
  //Need to figure out what to do with multiple active fibers

  std::map<std::string, std::vector<size_t> > goodTriggers;
  for(size_t ip = 0; ip < fPairedStraightDevices.size(); ++ip){
    std::string name = fPairedStraightDevices[ip].first;
    std::vector<size_t> triggers;
    //std::cout << ip << " " << name << " " << fPairedStraightDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;

    name = fPairedStraightDevices[ip].second;
    triggers.clear();
    //std::cout << ip << " " << name << " " << fPairedStraightDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;
  }
   
  for(size_t ip = 0; ip < fPairedStraightDevices.size(); ++ip){
    std::string nameOne = fPairedStraightDevices[ip].first;
    std::string nameTwo = fPairedStraightDevices[ip].second;

    size_t triggerSizeOne = goodTriggers[nameOne].size();
    size_t triggerSizeTwo = goodTriggers[nameTwo].size();

    if(triggerSizeOne < 1){
      //std::cout << "Missing trigger for " << nameOne << std::endl;
    }
    if(triggerSizeTwo < 1){
      //std::cout << "Missing trigger for " << nameTwo << std::endl;
    }

    if(triggerSizeOne > 0 && triggerSizeTwo > 0){
      //std::cout << "Found good triggers for " << nameOne << " " << nameTwo << std::endl;
      size_t triggerOne = goodTriggers[nameOne][0];
      size_t triggerTwo = goodTriggers[nameTwo][0];
      fBeamProf2D[ip]->Fill(beamevt.GetActiveFibers(nameOne,triggerOne)[0], beamevt.GetActiveFibers(nameTwo,triggerTwo)[0]);
    }
  }

}

void proto::BeamAna::GetUnpairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time){
  //This method goes through the unpaired devices, as well as the paired devices individually
  //and gets the info that matches to the inputted time
  //
  //Currently allows for multiple fibers in an event


  //Going through paired devices individually
  for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
    std::string name = fPairedDevices[ip].first;
    //std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iFiber = 0; iFiber < beamevt.GetActiveFibers(name,itN).size(); ++iFiber){
	  fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iFiber]);
	} // endfor 
      }
    }

    name = fPairedDevices[ip].second;
    //std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iFiber = 0; iFiber < beamevt.GetActiveFibers(name,itN).size(); ++iFiber){
	  fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iFiber]);
	} // endfor 
      }
    }
  }

  //Now going through unpaired
  for(size_t id = 0; id < fDevices.size(); ++id){
    std::string name = fDevices[id];
    //std::cout << id << " " << name << " " << fDevices[id] << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iFiber = 0; iFiber < beamevt.GetActiveFibers(name,itN).size(); ++iFiber){
	  fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iFiber]);
	} // endfor
      }
    }
  }

}

TVector3 proto::BeamAna::ProjectToTPC(TVector3 firstPoint, TVector3 secondPoint){
  TVector3 dR = (secondPoint - firstPoint);
  
  double deltaZ = -1.*secondPoint.Z();
  double deltaX = deltaZ * (dR.X() / dR.Z());
  double deltaY = deltaZ * (dR.Y() / dR.Z());

  TVector3 lastPoint = secondPoint + TVector3(deltaX, deltaY, deltaZ);
  return lastPoint;
}

double proto::BeamAna::GetPosition(std::string deviceName, int fiberIdx){
  //NEEDS WORK
  if(fiberIdx > 192){ std::cout << "Please select fiber in range [0,191]" << std::endl; return -1.;}
  double size = fFiberDimension[deviceName];
  //double size = 1.;
  
  //Define 0th fiber as farthest positive. Last fiber is farthest negative. Center is between 96 and 97 
  double pos = size*(96 - fiberIdx) - size/2.;
  return pos;
}

TVector3 proto::BeamAna::TranslateDeviceToDetector(TVector3 globalDeviceCoords){
  //fGlobalDetCoords given by fcl
  //Translate position of device w.r.t. Detector
  TVector3 inDetCoords = globalDeviceCoords - fGlobalDetCoords;
   
  //Rotate into detector coordinates
  inDetCoords.RotateX(fDetRotation[0]);
  inDetCoords.RotateY(fDetRotation[1]);
  inDetCoords.RotateZ(fDetRotation[2]);
  return inDetCoords;
}


DEFINE_ART_MODULE(proto::BeamAna)
