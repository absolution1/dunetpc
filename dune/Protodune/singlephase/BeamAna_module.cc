////////////////////////////////////////////////////////////////////////
// Class:       BeamAna
// Plugin Type: producer (art v2_08_03)
// File:        BeamAna_module.cc
//
// Generated at Thu Nov  2 22:57:41 2017 by Jonathan Paley using cetskelgen
// from cetlib version v3_01_01.
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
  void matchStraightTriggers(beam::ProtoDUNEBeamEvent beamevt,double);
//  void matchCurvedTriggers();
  void GetPairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);
  void GetUnpairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);
  double GetPosition(std::string, size_t);
  TVector3 TranslateDeviceToDetector(TVector3);
 
  void  InitXBPFInfo(beam::ProtoDUNEBeamEvent *);
  void  parseXBPF(uint64_t);
  void  parsePairedXBPF(uint64_t);

  void  parseXTOF(uint64_t);
  void  parseXCET(uint64_t);

  void  InitXCETInfo(beam::ProtoDUNEBeamEvent *);
  
private:
  
  TTree * fOutTree;
  TH2F * fFirstBeamProf2D;
  TH2F * fSecondBeamProf2D;
  std::vector<TH2F *> fBeamProf2D;
  std::map<std::string, TH1F *> fBeamProf1D;
  std::map<std::string, std::vector<short> *> fActiveFibers;
  std::map<std::string, double> fProfTime;
  std::map<std::string, TTree*> fProfTree;
  TH1F * fTOFHist;
  TH1F * fCKovHist;
  recob::TrackTrajectory theTraj;
  recob::Track * theTrack;
  

  // Declare member data here.
  //bool fLoadFromDB; // unused
  double  fTimeWindow;
  double fTolerance;
  std::string fCSVFileName;
  std::string fBundleName;
  double fDummyEventTime;
  std::string fURLStr;
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
  std::string fXBPFPrefix;
  std::string fXTOFPrefix;
  std::string fXCETPrefix;

  std::string fTOF1;
  std::string fTOF2;
  
  std::string fCKov1;
  std::string fCKov2;

  std::map<std::string, size_t> deviceOrder;
  beam::ProtoDUNEBeamEvent * beamevt;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
};


proto::BeamAna::BeamAna(fhicl::ParameterSet const & p)
//  :
//  EDProducer(p)  // ,
 // More initializers here.
{
  produces<std::vector<beam::ProtoDUNEBeamEvent>>();  
  this->reconfigure(p);
}

void proto::BeamAna::produce(art::Event & e)
{

  std::cerr << "%%%%%%%%%% Getting ifbeam service handle %%%%%%%%%%" << std::endl;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifb;
  std::cerr << "%%%%%%%%%% Got ifbeam handle %%%%%%%%%%" << std::endl;

  bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
  std::cerr << "%%%%%%%%%% Got beam folder %%%%%%%%%%" << std::endl;

  //Use multiple times provided to fcl
  if( fMultipleTimes.size() ){
    std::cout << "Using multiple times: " << std::endl;
    for(size_t it = 0; it < fMultipleTimes.size(); ++it){
      std::cout << fMultipleTimes[it] << std::endl;
    }
  }
  //Use singular time provided to fcl
  else if(fFixedTime){
    std::cout << "Using singular time: " << fFixedTime << std::endl;
    fMultipleTimes.push_back(fFixedTime);
  }
  //Use time of event
  else{
    std::cout <<" Using Event Time: " << uint64_t(e.time().timeLow()) << std::endl;
    fMultipleTimes.push_back(uint64_t(e.time().timeLow()));
  }

  //Start getting beam event info
  beamevt = new beam::ProtoDUNEBeamEvent();
  InitXBPFInfo(beamevt);

  for(size_t it = 0; it < fMultipleTimes.size(); ++it){
    std::cout << "Time: " << fMultipleTimes[it] << std::endl;
  //  parseXBPF(fMultipleTimes[it]);
  //  parsePairedXBPF(fMultipleTimes[it]);

    parseXCET(fMultipleTimes[it]);
  }


 //Setting some dummy Triggers for drawing  


 std::map<std::string,double> dummyTriggerTimeLSB = {{"XBPF022707",1.50000e+08},{"XBPF022708",1.50000e+08},{"XBPF022716",1.50002e+08},{"XBPF022717",1.50002e+08}};
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

//   size_t id = deviceOrder[name];

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

//   id = deviceOrder[name];

   beamevt->AddFBMTrigger(name,fbm);
   beamevt->DecodeFibers(name,beamevt->GetNFBMTriggers(name) - 1);//Decode the last one

   i = beamevt->GetNFBMTriggers(name) - 1;
   std::cout << name << " has active fibers: ";
   for(size_t iF = 0; iF < beamevt->GetActiveFibers(name,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(name, i)[iF] << " "; 
   std::cout << std::endl;
 }
 
 std::cout << "Matching" << std::endl;
 matchStraightTriggers(*beamevt,fDummyEventTime);
 std::cout << "Matched" << std::endl;
 
 std::cout << "Pairing" << std::endl;
 //GetPairedFBMInfo(*beamevt,1.50000e+08);
 GetPairedFBMInfo(*beamevt,fDummyEventTime);
 std::cout << "Paired" << std::endl;

 std::cout << "Unpaired" << std::endl;
 GetUnpairedFBMInfo(*beamevt,fDummyEventTime);
 std::cout << "Unpaired" << std::endl;

 std::cout << "Event stuff" << std::endl;
 std::unique_ptr<std::vector<beam::ProtoDUNEBeamEvent> > beamData(new std::vector<beam::ProtoDUNEBeamEvent>);
 beamData->push_back(beam::ProtoDUNEBeamEvent(*beamevt));
 std::cout << "Putting" << std::endl;
 e.put(std::move(beamData));
 std::cout << "Put" << std::endl;
}

void proto::BeamAna::InitXBPFInfo(beam::ProtoDUNEBeamEvent * beamevt){
  
  //Places a dummy trigger vector for each device
  size_t nDev = 0; 
  std::vector<std::string> monitors;


  for(size_t id = 0; id < fDevices.size(); ++id){
    std::cout << fXBPFPrefix + fDevices[id] << std::endl;
    std::string name = fDevices[id];
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;
    std::cout << "Rotated: " << fRotations[name][0] << " " << fRotations[name][1] << " " << fRotations[name][2] << std::endl;
    
    monitors.push_back(name);

    deviceOrder[name] = nDev; 
    nDev++;
  }

  for(size_t id = 0; id < fPairedStraightDevices.size(); ++id){

    std::string name = fPairedStraightDevices[id].first;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    deviceOrder[name] = nDev;
    std::cout << nDev << std::endl;
    nDev++;
    monitors.push_back(name);

    name = fPairedStraightDevices[id].second;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;
    monitors.push_back(name);
    deviceOrder[name] = nDev;
    std::cout << nDev << std::endl;
    nDev++;
  }

  for(size_t id = 0; id < fPairedDevices.size(); ++id){

    std::string name = fPairedDevices[id].first;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    deviceOrder[name] = nDev;
    std::cout << nDev << std::endl;
    nDev++;
    monitors.push_back(name);

    name = fPairedDevices[id].second;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    deviceOrder[name] = nDev;
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

void proto::BeamAna::parseXTOF(uint64_t time){
  std::cout << "Getting TOF1 info: " << fTOF1 << std::endl;
  std::vector<double> dataTOF1 = bfp->GetNamedVector(time, fXTOFPrefix + fTOF1 + ":countsTrig[]");
  std::cout << "Size of countsTrig: " << dataTOF1.size(); 

  std::cout << "Getting TOF2 info: " << fTOF2 << std::endl;
  std::vector<double> dataTOF2 = bfp->GetNamedVector(time, fXTOFPrefix + fTOF2 + ":countsTrig[]");
  std::cout << "Size of countsTrig: " << dataTOF2.size(); 
}

void proto::BeamAna::parseXCET(uint64_t time){
  std::cout << "Getting CKov1 info: " << fCKov1 << std::endl;
  std::vector<double> dataCKov1 = bfp->GetNamedVector(time, fXCETPrefix + fCKov1 + ":countsTrig");
  std::cout << "countsTrig: " << dataCKov1[0] << std::endl; 

  std::cout << "Getting CKov2 info: " << fCKov2 << std::endl;
  std::vector<double> dataCKov2 = bfp->GetNamedVector(time, fXCETPrefix + fCKov2 + ":countsTrig");
  std::cout << "countsTrig: " << dataCKov2[0] << std::endl; 
}

void proto::BeamAna::parseXBPF(uint64_t time){
  for(size_t d = 0; d < fDevices.size(); ++d){
    std::string name = fDevices[d];
//    size_t iDevice = deviceOrder[name];
    std::cout <<"Device: " << name << std::endl;
    std::vector<double> data = bfp->GetNamedVector(time, fXBPFPrefix + name + ":eventsData[]");
    std::vector<double> counts = bfp->GetNamedVector(time, fXBPFPrefix + name + ":countsRecords[]");
    

    std::cout << "Data: " << data.size() << std::endl;
    std::cout << "Counts: " << counts.size() << std::endl;
    for(size_t i = 0; i < counts.size(); ++i){
      std::cout << counts[i] << std::endl;
    }
    if(counts[1] > data.size())continue;

    beam::FBM fbm;
    fbm.ID = d;
         
    for(size_t i = 0; i < counts[1]; ++i){
//      std::cout << "Count: " << i << std::endl;
      for(int j = 0; j < 10; ++j){
        double theData = data[20*i + (2*j + 1)];
//        std::cout << std::setw(15) << theData ;
        if(j < 4){
          fbm.timeData[j] = theData;           
        }
        else{
          fbm.fiberData[j - 4] = theData;
        }
      }
      beamevt->AddFBMTrigger(name ,fbm);
//      std::cout << std::endl;
    }

    for(size_t i = 0; i < beamevt->GetNFBMTriggers(name); ++i){
      beamevt->DecodeFibers(name,i);
      std::cout << name << " at time: " << beamevt->DecodeFiberTime(name, i) << " has active fibers: ";
      for(size_t iF = 0; iF < beamevt->GetActiveFibers(name,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(name, i)[iF] << " "; 
      std::cout << std::endl;
        
      *fActiveFibers[name] = beamevt->GetActiveFibers(name,i);
      fProfTime[name] = beamevt->DecodeFiberTime(name, i);
      fProfTree[name]->Fill(); 

    }
  }
  
}


void proto::BeamAna::parsePairedXBPF(uint64_t time){
  for(size_t d = 0; d < fPairedDevices.size(); ++d){
    std::string name = fPairedDevices[d].first;
//    size_t iDevice = deviceOrder[name];

    std::cout <<"Device: " << name <</* " " << iDevice <<*/ std::endl;
    std::vector<double> data = bfp->GetNamedVector(time, fXBPFPrefix + name + ":eventsData[]");
    std::vector<double> counts = bfp->GetNamedVector(time, fXBPFPrefix + name + ":countsRecords[]");
  

    std::cout << "Data: " << data.size() << std::endl;
    std::cout << "Counts: " << counts.size() << std::endl;
    for(size_t i = 0; i < counts.size(); ++i){
      std::cout << counts[i] << std::endl;
    }
    if(counts[1] > data.size())continue;

    beam::FBM fbm;
    fbm.ID = d;
         
    for(size_t i = 0; i < counts[1]; ++i){
//      std::cout << "Count: " << i << std::endl;
      for(int j = 0; j < 10; ++j){
        double theData = data[20*i + (2*j + 1)];
//        std::cout << std::setw(15) << theData ;
        if(j < 4){
          fbm.timeData[j] = theData;           
        }
        else{
          fbm.fiberData[j - 4] = theData;
        }
      }
      beamevt->AddFBMTrigger(name ,fbm);
//      std::cout << std::endl;
    }

    for(size_t i = 0; i < beamevt->GetNFBMTriggers(name); ++i){
      beamevt->DecodeFibers(name,i);
      std::cout << name << " at time: " << beamevt->DecodeFiberTime(name, i) << " has active fibers: ";
      for(size_t iF = 0; iF < beamevt->GetActiveFibers(name,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(name, i)[iF] << " "; 
      std::cout << std::endl;

      *fActiveFibers[name] = beamevt->GetActiveFibers(name,i);
      fProfTime[name] = beamevt->DecodeFiberTime(name, i);
      fProfTree[name]->Fill(); 
    }

    name = fPairedDevices[d].second;
//    iDevice = deviceOrder[name];
    std::cout <<"Device: " << name << std::endl;
    data = bfp->GetNamedVector(time, fXBPFPrefix + name + ":eventsData[]");
    counts = bfp->GetNamedVector(time, fXBPFPrefix + name + ":countsRecords[]");
 
    std::cout << "Data: " << data.size() << std::endl;
    std::cout << "Counts: " << counts.size() << std::endl;
    for(size_t i = 0; i < counts.size(); ++i){
      std::cout << counts[i] << std::endl;
    }
    if(counts[1] > data.size())continue;


    //Sorry for just repeating this, I'll need more time to refactor it
    //Reset the FBM
    fbm = beam::FBM();
    fbm.ID = d;
         
    for(size_t i = 0; i < counts[1]; ++i){
      //std::cout << "Count: " << i << std::endl;
      for(int j = 0; j < 10; ++j){
        double theData = data[20*i + (2*j + 1)];
        //std::cout << std::setw(15) << theData ;
        if(j < 4){
          fbm.timeData[j] = theData;           
        }
        else{
          fbm.fiberData[j - 4] = theData;
        }
      }
      beamevt->AddFBMTrigger(name ,fbm);
      //std::cout << std::endl;
    }

    for(size_t i = 0; i < beamevt->GetNFBMTriggers(name); ++i){
      beamevt->DecodeFibers(name,i);
      std::cout << name << " at time: " << beamevt->DecodeFiberTime(name, i) << " has active fibers: ";
      for(size_t iF = 0; iF < beamevt->GetActiveFibers(name,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(name, i)[iF] << " "; 
      std::cout << std::endl;

      *fActiveFibers[name] = beamevt->GetActiveFibers(name,i);
      fProfTime[name] = beamevt->DecodeFiberTime(name, i);
      fProfTree[name]->Fill(); 
    }
  }
}

void proto::BeamAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  

  fTOFHist = tfs->make<TH1F>("TOF","",100,0,100);
  fCKovHist = tfs->make<TH1F>("CKov","",4,0,4);

  fOutTree = tfs->make<TTree>("tree", "lines"); 
  theTrack = new recob::Track();
  fOutTree->Branch("Track", &theTrack);

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
    fProfTree[fPairedDevices[i].first] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedDevices[i].first]->Branch("time", &fProfTime[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("fibers", &fActiveFibers[fPairedDevices[i].first]);

    name = "Fibers_" + fPairedDevices[i].second;
    fActiveFibers[fPairedDevices[i].second] = new std::vector<short>;
    fProfTime[fPairedDevices[i].second] = 0.;
    fProfTree[fPairedDevices[i].second] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedDevices[i].second]->Branch("time", &fProfTime[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("fibers", &fActiveFibers[fPairedDevices[i].second]);
  }

  for(size_t i = 0; i < fDevices.size(); ++i){
    std::string name = "BeamProf1D_" + fDevices[i];
    std::string title = fDevices[i];
    fBeamProf1D[fDevices[i]] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );

    name = "Fibers_" + fDevices[i];
    fProfTree[fDevices[i]] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTime[fDevices[i]] = 0.;
    fActiveFibers[fDevices[i]] = new std::vector<short>;
    fProfTree[fDevices[i]]->Branch("time", &fProfTime[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("fibers", &fActiveFibers[fDevices[i]]);
  }
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
  fURLStr      = "";
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

  fCKov1 = p.get< std::string >("CKov1");
  fCKov2 = p.get< std::string >("CKov2");


  fXBPFPrefix      = p.get<std::string>("XBPFPrefix");
  fXTOFPrefix      = p.get<std::string>("XTOFPrefix");
  fXCETPrefix      = p.get<std::string>("XCETPrefix");
  fDummyEventTime = p.get<double>("DummyEventTime");
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

void proto::BeamAna::matchStraightTriggers(beam::ProtoDUNEBeamEvent beamevt, double Time){
  
  std::vector<TVector3> thePoints;
  std::vector<TVector3> theMomenta;
  std::vector< std::vector<double> > dummy;
  std::vector<double> mom(3, util::kBogusD);

  std::map<std::string, size_t> goodTriggers;
  bool foundNext = false;

  //Start with the first in the first pair, this should be the most upstream
  //Next will be the second of the first pair, then the first of the second pair. 
  //Go through each trigger and try to match

  std::string firstName = fPairedStraightDevices[0].first;
  size_t firstIndex = deviceOrder[firstName];
  std::cout << firstName << " " << firstIndex << std::endl;
  for(size_t it = beamevt.GetNFBMTriggers(firstName) - 1; it < beamevt.GetNFBMTriggers(firstName); ++it){
    if ( ( fabs(beamevt.DecodeFiberTime(firstName, it ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(firstName, it ) - Time) >= 0 ) ){
      foundNext = false;
      double firstTime = beamevt.DecodeFiberTime(firstName, it);
      std::cout << "Trigger " << it << " Time " << firstTime << std::endl;

      //Go through the next FBMs and see if they have a good Time
      std::string secondName = fPairedStraightDevices[0].second;
      
      for( size_t it2 = 0; it2 < beamevt.GetNFBMTriggers(secondName); ++it2){
//        std::cout << "\tTime: " << beamevt.DecodeFiberTime(secondName, it2 ) << std::endl;
        if ( ( fabs(beamevt.DecodeFiberTime(secondName, it2 ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(secondName, it2 ) - Time) >= 0 ) ){
          std::cout << "Found it" << std::endl;
          foundNext = true;
          goodTriggers[firstName] = it;
          goodTriggers[secondName] = it2;

          //For now: break after first 
          break;
        }
      }
      if(!foundNext)continue;

//      for(size_t ip = 1; ip < fPairedStraightDevices.size(); ++ip){
        foundNext = false;
//        std::string nextName = fPairedStraightDevices[ip].first;
        std::string nextName = fPairedStraightDevices[1].first;
//        size_t nextIndex = deviceOrder[nextName]; 
        std::cout << "\tNext first " << nextName << /*" " << nextIndex <<*/ std::endl;

        for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(nextName); ++itN){
      //    std::cout << "\tTime: " << beamevt.DecodeFiberTime(nextName, itN ) << std::endl;
          if ( ( fabs(beamevt.DecodeFiberTime(nextName, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(nextName, itN ) - Time) >= 0 ) ){
            std::cout << "Found it" << std::endl;
            foundNext = true;
            goodTriggers[nextName] = itN;
            break;
          }
        }
        if(!foundNext){
          goodTriggers.clear();
          //break;
          continue;
        }
        
        //Move on to second
//        nextName = fPairedStraightDevices[ip].second;
        nextName = fPairedStraightDevices[1].second;
//        nextIndex = deviceOrder[nextName];
        std::cout << "\tNext second " << nextName << /*" " << nextName <<*/ std::endl;
        for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(nextName); ++itN){
       //   std::cout << "\tTime: " << beamevt.DecodeFiberTime(nextName, itN ) << std::endl;
          if ( ( fabs(beamevt.DecodeFiberTime(nextName, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(nextName, itN ) - Time) >= 0 ) ){
            std::cout << "Found it" << std::endl;
            foundNext = true;
            goodTriggers[nextName] = itN;
            break;
          }
        }
        if(!foundNext){
          goodTriggers.clear();
          //break;
          continue;
        }
      //}


      //Found a full track
      if(foundNext){
        std::cout << "Found good track" << std::endl;

        std::vector<double> xPos;      
        std::vector<double> yPos;      

        std::map<std::string, std::vector<double>*> posArray = {{"horiz",&xPos},{"vert",&yPos}};

        for(size_t ip = 0; ip < fPairedStraightDevices.size(); ++ip){
          //Start with the first device
          //Which active Fibers are on?       
          std::string name = fPairedStraightDevices[ip].first;
//          size_t index = deviceOrder[name];
          std::vector<short> active = beamevt.GetActiveFibers(name, goodTriggers[name]);

          //Gets position within the FBM for the first active
          double position = GetPosition(name, active[0]);

          //Checks if vertical or horizontal, pushes to the correct array 
          posArray[fDeviceTypes[name]]->push_back(position); 


          //Do the same for the second device
          name = fPairedStraightDevices[ip].second;
//          index = deviceOrder[name];
          active = beamevt.GetActiveFibers(name, goodTriggers[name]);

          position = GetPosition(name, active[0]);
          posArray[fDeviceTypes[name]]->push_back(position); 

          //Now has an x,y point within detector pair
          std::cout << "Sizes: " << xPos.size() << " " << yPos.size() << std::endl;

          //Translate to device's global coords
          //Rotation? Deal with later
          double devZ = 0.;//fCoordinates[name][2];
          double devY = yPos.at(ip);// + fCoordinates[name][1];
          double devX = xPos.at(ip);// + fCoordinates[name][0];

          TVector3 coordsInDevice(devX,devY,devZ);
          TVector3 rotatedCoords = coordsInDevice;
          rotatedCoords.RotateX(fRotations[name][0]);
          rotatedCoords.RotateY(fRotations[name][1]);
          rotatedCoords.RotateZ(fRotations[name][2]);
          
          TVector3 deviceShift(fCoordinates[name][0],fCoordinates[name][1],fCoordinates[name][2]);

          TVector3 coordsInWorld = rotatedCoords - deviceShift;
               

//          thePoints.push_back( TVector3(xPos.at(ip), yPos.at(ip), fCoordinates[name][2]) );  
//          thePoints.push_back( TranslateDeviceToDetector(TVector3(devX,devY,devZ)) );     
          thePoints.push_back( coordsInWorld );
        //  theMomenta.push_back( TVector3(0.,0.,1.) );
        }

        //Extending track to TPC face
        TVector3 dR( (thePoints[1].X() - thePoints[0].X()), (thePoints[1].Y() - thePoints[0].Y()), (thePoints[1].X() - thePoints[0].Z()));
        
        //2 for now
        theMomenta.push_back(dR.Unit());
        theMomenta.push_back(dR.Unit());

        double dX = dR.X();
        double dY = dR.Y();
        double dZ = dR.Z();

        double tanXZ = dX/dZ;
        double tanYZ = dY/dZ;
        
        double newDZ = 0. - thePoints[1].Z();
        double newX = thePoints[1].X() + newDZ*tanXZ;  
        double newY = thePoints[1].Y() + newDZ*tanYZ;  

        thePoints.push_back( TVector3( newX, newY, 0.)  );          
        theMomenta.push_back( dR.Unit() );

        std::cout << "Making traj" << std::endl;
        *theTrack = recob::Track(thePoints, theMomenta, dummy, mom, 1);
        std::cout << " Done " << std::endl;
        fOutTree->Fill();
      }
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
//    size_t index = deviceOrder[name];
    std::vector<size_t> triggers;
    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;

    name = fPairedDevices[ip].second;
//    index = deviceOrder[name];
    triggers.clear();
    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;
  }
   
  for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
    std::string nameOne = fPairedDevices[ip].first;
    std::string nameTwo = fPairedDevices[ip].second;
//    size_t indexOne = deviceOrder[nameOne];
//    size_t indexTwo = deviceOrder[nameTwo];

    size_t triggerSizeOne = goodTriggers[nameOne].size();
    size_t triggerSizeTwo = goodTriggers[nameTwo].size();

    if(triggerSizeOne < 1){
      std::cout << "Missing trigger for " << nameOne << std::endl;
    }
    if(triggerSizeTwo < 1){
      std::cout << "Missing trigger for " << nameTwo << std::endl;
    }

    if(triggerSizeOne > 0 && triggerSizeTwo > 0){
      std::cout << "Found good triggers for " << nameOne << " " << nameTwo << std::endl;
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
//    size_t name = deviceOrder[name];
    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iF = 0; iF < beamevt.GetActiveFibers(name,itN).size(); ++iF) fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iF]); 
      }
    }

    name = fPairedDevices[ip].second;
//    name = deviceOrder[name];
    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iF = 0; iF < beamevt.GetActiveFibers(name,itN).size(); ++iF) fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iF]); 
      }
    }
  }

  //Now going through unpaired
  for(size_t id = 0; id < fDevices.size(); ++id){
    std::string name = fDevices[id];
//    size_t name = deviceOrder[name];
    std::cout << id << " " << name << " " << fDevices[id] << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN ) - Time)     >= 0 ) ){
        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iF = 0; iF < beamevt.GetActiveFibers(name,itN).size(); ++iF) fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iF]); 
      }
    }
  }

}
double proto::BeamAna::GetPosition(std::string deviceName, size_t iFiber){
  //NEEDS WORK
  if(iFiber > 192){ std::cout << "Please select fiber in range [0,191]" << std::endl; return -1.;}
 // double size = fFiberDimension[deviceName];
  double size = 1.;
  
  //Define 0th fiber as farthest positive. Last fiber is farthest negative. Center is between 96 and 97 
  double pos = size*iFiber + size/2.;
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
