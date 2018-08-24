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
  void matchTriggers(beam::ProtoDUNEBeamEvent beamevt);
  void GetPairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);
  double GetPosition(size_t, size_t);

  void  parseDevices(uint64_t);
  void  parsePairedDevices(uint64_t);

private:
  
  TTree * fOutTree;
  TH2F * fFirstBeamProf2D;
  TH2F * fSecondBeamProf2D;
  TH1F * fTOFHist;
  TH1F * fCKovHist;
  recob::TrackTrajectory theTraj;
  recob::Track theTrack;
  

  // Declare member data here.
  //bool fLoadFromDB; // unused
  double  fTimeWindow;
  double fTolerance;
  std::string fCSVFileName;
  std::string fBundleName;
  std::string fURLStr;
  uint64_t fFixedTime;
  std::vector< uint64_t > fMultipleTimes;
  std::vector< std::string > fDevices;
  std::vector< std::pair<std::string, std::string> > fPairedDevices;
  std::map<std::string, std::string > fDeviceTypes;
//  std::vector< std::array<double, 3> > fCoordinates; 
  std::map< std::string, std::array<double,3> > fCoordinates;
  std::vector< std::array<double, 3> > fRotations;
  std::vector< double > fFiberDimension;
  std::string fPrefix;

  std::map<std::string, size_t> deviceOrder;
  beam::ProtoDUNEBeamEvent * beamevt;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
};


proto::BeamAna::BeamAna(fhicl::ParameterSet const & p)
//  :
//  EDProducer(p)  // ,
 // More initializers here.
{
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

  size_t nDev = 0; 


  for(size_t id = 0; id < fDevices.size(); ++id){
    std::cout << fPrefix + fDevices[id] << std::endl;
    std::string name = fDevices[id];
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;
    std::cout << "Rotated: " << fRotations[id][0] << " " << fRotations[id][1] << " " << fRotations[id][2] << std::endl;
    
    deviceOrder[fDevices[id]] = nDev; 
    nDev++;
  }

  for(size_t id = 0; id < fPairedDevices.size(); ++id){

    std::string name = fPairedDevices[id].first;
    std::cout << fPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    deviceOrder[name] = nDev;
    std::cout << nDev << std::endl;
    nDev++;

    name = fPairedDevices[id].second;
    std::cout << fPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    deviceOrder[name] = nDev;
    std::cout << nDev << std::endl;
    nDev++;
  } 
    
  std::cout << "Using " << nDev << " devices: " << std::endl;
  std::map<std::string, size_t>::iterator itD = deviceOrder.begin();
  for(; itD != deviceOrder.end(); ++itD){
    std::cout << itD->first << " " << itD->second << std::endl;
  }

  //Start getting beam event information
  beamevt = new beam::ProtoDUNEBeamEvent();
  beamevt->InitFBMs(nDev);


  for(size_t it = 0; it < fMultipleTimes.size(); ++it){
    std::cout << "Time: " << fMultipleTimes[it] << std::endl;
    parseDevices(fMultipleTimes[it]);
    parsePairedDevices(fMultipleTimes[it]);
  }


 //Setting some dummy Triggers for drawing  


 std::map<std::string,double> dummyTriggerTimeLSB = {{"XBPF022697_V",1.50000e+08},{"XBPF022697_H",1.50000e+08},{"XBPF022707_V",1.50002e+08},{"XBPF022707_H",1.50002e+08}};
 std::map<std::string,double> dummyTriggerTimeMSB = {{"XBPF022697_V",1.53191e+09},{"XBPF022697_H",1.53191e+09},{"XBPF022707_V",1.53191e+09},{"XBPF022707_H",1.53191e+09}};
 std::map<std::string,double> dummyEventTimeLSB   = {{"XBPF022697_V",1.50000e+08},{"XBPF022697_H",1.50000e+08},{"XBPF022707_V",1.50000e+08},{"XBPF022707_H",1.50000e+08}};
 std::map<std::string,double> dummyEventTimeMSB   = {{"XBPF022697_V",1.53191e+09},{"XBPF022697_H",1.53191e+09},{"XBPF022707_V",1.53191e+09},{"XBPF022707_H",1.53191e+09}};


 std::map<std::string,std::array<double,6>> dummyHit     = { {"XBPF022697_V",{1,0,0,0,0,0}},
                                                            {"XBPF022697_H",{1,0,0,0,0,0}},
                                                            {"XBPF022707_V",{0,0,0,0,1,0}},
                                                            {"XBPF022707_H",{0,0,0,0,1,0}} }; 

 for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
   beam::FBM fbm;
   fbm.ID = ip;
   std::string name = fPairedDevices[ip].first; 
   fbm.timeData[0] = dummyTriggerTimeLSB[name];
   fbm.timeData[1] = dummyTriggerTimeMSB[name];
   fbm.timeData[2] = dummyEventTimeLSB[name];
   fbm.timeData[3] = dummyEventTimeMSB[name];

   for(int i = 0; i < 6; ++i){
     fbm.fiberData[i] = dummyHit[name][i];
   }

   size_t id = deviceOrder[name];

   beamevt->AddFBMTrigger(id,fbm);
   beamevt->DecodeFibers(id,beamevt->GetNFBMTriggers(id) - 1);//Decode the last one

   size_t i = beamevt->GetNFBMTriggers(id) - 1;
   std::cout << name << " has active fibers: ";
   for(size_t iF = 0; iF < beamevt->GetActiveFibers(id,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(id, i)[iF] << " "; 
   std::cout << std::endl;

   fbm = beam::FBM();
   fbm.ID = ip;
   name = fPairedDevices[ip].second; 
   fbm.timeData[0] = dummyTriggerTimeLSB[name];
   fbm.timeData[1] = dummyTriggerTimeMSB[name];
   fbm.timeData[2] = dummyEventTimeLSB[name];
   fbm.timeData[3] = dummyEventTimeMSB[name];

   for(int i = 0; i < 6; ++i){
     fbm.fiberData[i] = dummyHit[name][i];
   }

   id = deviceOrder[name];

   beamevt->AddFBMTrigger(id,fbm);
   beamevt->DecodeFibers(id,beamevt->GetNFBMTriggers(id) - 1);//Decode the last one

   i = beamevt->GetNFBMTriggers(id) - 1;
   std::cout << name << " has active fibers: ";
   for(size_t iF = 0; iF < beamevt->GetActiveFibers(id,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(id, i)[iF] << " "; 
   std::cout << std::endl;
 }

 matchTriggers(*beamevt);
 
 GetPairedFBMInfo(*beamevt,1.50000e+08);

 std::unique_ptr<std::vector<beam::ProtoDUNEBeamEvent> > beamData(new std::vector<beam::ProtoDUNEBeamEvent>);
 beamData->push_back(beam::ProtoDUNEBeamEvent(*beamevt));

 e.put(std::move(beamData));
}

void proto::BeamAna::parseDevices(uint64_t time){
  for(size_t d = 0; d < fDevices.size(); ++d){
    std::string name = fDevices[d];
    size_t iDevice = deviceOrder[name];
    std::cout <<"Device: " << name << std::endl;
    std::vector<double> data = bfp->GetNamedVector(time, fPrefix + name + ":eventsData[]");
    std::vector<double> counts = bfp->GetNamedVector(time, fPrefix + name + ":countsRecords[]");
    

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
      beamevt->AddFBMTrigger(iDevice ,fbm);
//      std::cout << std::endl;
    }

    for(size_t i = 0; i < beamevt->GetNFBMTriggers(iDevice); ++i){
      beamevt->DecodeFibers(iDevice,i);
//      std::cout << name << " has active fibers: ";
 //     for(size_t iF = 0; iF < beamevt->GetActiveFibers(iDevice,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(iDevice, i)[iF] << " "; 
//      std::cout << std::endl;
    }
  }
  
}

void proto::BeamAna::parsePairedDevices(uint64_t time){
  for(size_t d = 0; d < fPairedDevices.size(); ++d){
    std::string name = fPairedDevices[d].first;
    size_t iDevice = deviceOrder[name];

    std::cout <<"Device: " << name << " " << iDevice << std::endl;
    std::vector<double> data = bfp->GetNamedVector(time, fPrefix + name + ":eventsData[]");
    std::vector<double> counts = bfp->GetNamedVector(time, fPrefix + name + ":countsRecords[]");
  

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
      beamevt->AddFBMTrigger(iDevice ,fbm);
//      std::cout << std::endl;
    }

    for(size_t i = 0; i < beamevt->GetNFBMTriggers(iDevice); ++i){
      beamevt->DecodeFibers(iDevice,i);
//      std::cout << name << " has active fibers: ";
//      for(size_t iF = 0; iF < beamevt->GetActiveFibers(iDevice,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(iDevice, i)[iF] << " "; 
//      std::cout << std::endl;
    }

    name = fPairedDevices[d].second;
    iDevice = deviceOrder[name];
    std::cout <<"Device: " << name << std::endl;
    data = bfp->GetNamedVector(time, fPrefix + name + ":eventsData[]");
    counts = bfp->GetNamedVector(time, fPrefix + name + ":countsRecords[]");
 
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
      std::cout << "Count: " << i << std::endl;
      for(int j = 0; j < 10; ++j){
        double theData = data[20*i + (2*j + 1)];
        std::cout << std::setw(15) << theData ;
        if(j < 4){
          fbm.timeData[j] = theData;           
        }
        else{
          fbm.fiberData[j - 4] = theData;
        }
      }
      beamevt->AddFBMTrigger(iDevice ,fbm);
      std::cout << std::endl;
    }

    for(size_t i = 0; i < beamevt->GetNFBMTriggers(iDevice); ++i){
      beamevt->DecodeFibers(iDevice,i);
      std::cout << name << " has active fibers: ";
      for(size_t iF = 0; iF < beamevt->GetActiveFibers(iDevice,i).size(); ++iF)std::cout << beamevt->GetActiveFibers(iDevice, i)[iF] << " "; 
      std::cout << std::endl;
    }
  }
}

void proto::BeamAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  fFirstBeamProf2D = tfs->make<TH2F>("FirstBeamProf2D","",192,0,192,192,0,192);
  fSecondBeamProf2D = tfs->make<TH2F>("SecondBeamProf2D","",192,0,192,192,0,192);

  fTOFHist = tfs->make<TH1F>("TOF","",100,0,100);
  fCKovHist = tfs->make<TH1F>("Cer","",4,0,4);

  fOutTree = tfs->make<TTree>("tree", "lines"); 
  //Need to make this configurable later
  // Implementation of optional member function here.
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

  fDevices     = p.get< std::vector< std::string > >("Devices");    
  fPairedDevices = p.get< std::vector< std::pair<std::string, std::string> > >("PairedDevices");

  std::vector< std::pair<std::string, std::string> >  tempTypes = p.get<std::vector< std::pair<std::string, std::string> >>("DeviceTypes");
  fDeviceTypes  = std::map<std::string, std::string>(tempTypes.begin(), tempTypes.end() );

  std::vector< std::pair<std::string, std::array<double,3> > > tempCoords = p.get<std::vector< std::pair<std::string, std::array<double,3> > > >("Coordinates");

  fTolerance   = p.get<double>("Tolerance");

  //Location of Device 
//  fCoordinates = p.get< std::vector< std::array<double,3> > >("Coordinates");
  fCoordinates = std::map<std::string, std::array<double,3> >(tempCoords.begin(), tempCoords.end());
 

  //Rotation of Device
  fRotations   = p.get< std::vector< std::array<double,3> > >("Rotations"); 

  //Deminsion of Fibers
  fFiberDimension = p.get< std::vector<double> >("Dimension");

  fPrefix      = p.get<std::string>("Prefix");
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

void proto::BeamAna::matchTriggers(beam::ProtoDUNEBeamEvent beamevt){
  
/*  recob::TrackTrajectory::Positions_t thePoints;
  recob::TrackTrajectory::Momenta_t theMomenta;
  recob::TrackTrajectory::Flags_t theFlags; */

  std::vector<TVector3> thePoints;
  std::vector<TVector3> theMomenta;
  std::vector< std::vector<double> > dummy;
  std::vector<double> mom(2, util::kBogusD);

  std::map<std::string, size_t> goodTriggers;
  bool foundNext = false;

  //Start with the first in the first pair, this should be the most upstream
  //Next will be the second of the first pair, then the first of the second pair, and so on. 
  //Go through each trigger and try to match

  std::string firstName = fPairedDevices[0].first;
  size_t firstIndex = deviceOrder[firstName];
  for(size_t it = 0; it < beamevt.GetNFBMTriggers(firstIndex); ++it){
 
    foundNext = false;
    double time = beamevt.DecodeFiberTime(firstIndex, it);

    //Go through the next FBMs and see if they have a good time
    
    std::string secondName = fPairedDevices[0].second;
    size_t secondIndex = deviceOrder[secondName];
    
    for( size_t it2 = 0; it2 < beamevt.GetNFBMTriggers(secondIndex); ++it2){
      if ( ( (beamevt.DecodeFiberTime(secondIndex, it2 ) - time) < /*0.00010e+08*/fTolerance ) && ( (beamevt.DecodeFiberTime(secondIndex, it2 ) - time) >= 0 ) ){
        foundNext = true;
        goodTriggers[firstName] = it;
        goodTriggers[secondName] = it2;

        //For now: break after first 
        break;
      }
    }
    if(!foundNext)continue;

    for(size_t ip = 1; ip < fPairedDevices.size(); ++ip){
      foundNext = false;
      std::string nextName = fPairedDevices[ip].first;
      size_t nextIndex = deviceOrder[nextName]; 

      for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(nextIndex); ++itN){
        if ( ( (beamevt.DecodeFiberTime(nextIndex, itN ) - time) < /*0.00010e+08*/fTolerance ) && ( (beamevt.DecodeFiberTime(nextIndex, itN ) - time) >= 0 ) ){
          foundNext = true;
          goodTriggers[nextName] = itN;
          break;
        }
      }
      if(!foundNext){
        break;
        goodTriggers.clear();
      }
    }
    if(foundNext){
      std::cout << "Found good track" << std::endl;
        

      std::vector<double> xPos;      
      std::vector<double> yPos;      

      std::map<std::string, std::vector<double>*> posArray = {{"vert",&xPos},{"horiz",&yPos}};

      for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
        std::string name = fPairedDevices[ip].first;
        size_t index = deviceOrder[name];
        std::vector<short> active = beamevt.GetActiveFibers(index, goodTriggers[name]);

        double position = GetPosition(index, active[0]);
        posArray[fDeviceTypes[name]]->push_back(position); 

        name = fPairedDevices[ip].second;
        index = deviceOrder[name];
        active = beamevt.GetActiveFibers(index, goodTriggers[name]);

        position = GetPosition(index, active[0]);
        posArray[fDeviceTypes[name]]->push_back(position); 
         
        std::cout << "Sizes: " << xPos.size() << " " << yPos.size() << std::endl;
        //Check if diff size?
        
        thePoints.push_back( TVector3(xPos.at(ip), yPos.at(ip), fCoordinates[name][2]) );  
        theMomenta.push_back( TVector3(0.,0.,1.) );
        //theFlags.push_back( recob::TrackTrajectory::PointFlags_t() );
      }

      std::cout << "Making traj" << std::endl;
      

      //theTraj(thePoints, theMomenta, theFlags, true);      
     // std::cout << "Length: " << theTraj.Length() << std::endl;
     // std::cout << "Sizes: " << xPos.size() << " " << yPos.size() << std::endl;
      theTrack = recob::Track(thePoints, theMomenta, dummy, mom, 1);
      std::cout << " Done " << std::endl;
      fOutTree->Fill();
    }

  }

}

//Gets info from FBMs matching in time
void proto::BeamAna::GetPairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time){

  std::map<std::string, std::vector<size_t> > goodTriggers;


  for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
    std::string name = fPairedDevices[ip].first;
    size_t index = deviceOrder[name];
    std::vector<size_t> triggers;
    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(index); ++itN){
      if ( ( (beamevt.DecodeFiberTime(index, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(index, itN ) - Time)     >= 0 ) ){
        std::cout << "Found Good Time" << name << " " << beamevt.DecodeFiberTime(index, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;

    name = fPairedDevices[ip].second;
    index = deviceOrder[name];
    triggers.clear();
    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(index); ++itN){
      if ( ( (beamevt.DecodeFiberTime(index, itN ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(index, itN ) - Time)     >= 0 ) ){
        std::cout << "Found Good Time" << name << " " << beamevt.DecodeFiberTime(index, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;
  }
  

  std::string nameOne = fPairedDevices[0].first;
  std::string nameTwo = fPairedDevices[0].second;
  size_t indexOne = deviceOrder[nameOne];
  size_t indexTwo = deviceOrder[nameTwo];

  size_t triggerOne = goodTriggers[nameOne][0];
  size_t triggerTwo = goodTriggers[nameTwo][0];
  //for(size_t iF = 0; iF < beamevt.GetActiveFibers(indexOne,triggerOne).size(); ++iF){
  //}
  fFirstBeamProf2D->Fill(beamevt.GetActiveFibers(indexOne,triggerOne)[0], beamevt.GetActiveFibers(indexTwo,triggerTwo)[0]);

  nameOne = fPairedDevices[1].first;
  nameTwo = fPairedDevices[1].second;
  indexOne = deviceOrder[nameOne];
  indexTwo = deviceOrder[nameTwo];
  triggerOne = goodTriggers[nameOne][0];
  triggerTwo = goodTriggers[nameTwo][0];
  //for(size_t iF = 0; iF < beamevt.GetActiveFibers(indexOne,triggerOne).size(); ++iF){
  //}
  fSecondBeamProf2D->Fill(beamevt.GetActiveFibers(indexOne,triggerOne)[0], beamevt.GetActiveFibers(indexTwo,triggerTwo)[0]);
}

double proto::BeamAna::GetPosition(size_t iDevice, size_t iFiber){
  if(iFiber > 192){ std::cout << "Please select fiber in range [0,191]" << std::endl; return -1.;}
//  if(iDevice > fDevices.size() -1 ){ std::cout << "Device out of range " << std::endl; return -1.;}
 // double size = fFiberDimension[iDevice];
  double size = 1.;
  
  //Define 0th fiber as origin. Go to middle of the fiber
  double pos = size*iFiber + size/2.;
  return pos;
}

DEFINE_ART_MODULE(proto::BeamAna)
