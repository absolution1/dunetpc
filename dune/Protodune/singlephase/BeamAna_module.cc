////////////////////////////////////////////////////////////////////////
// Class:       BeamAna
// Plugin Type: analyzer (art v2_08_03)
// File:        BeamAna_module.cc
//
// Generated at Thu Nov  2 22:57:41 2017 by Jonathan Paley using cetskelgen
// from cetlib version v3_01_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
//#include "dune/BeamData/ProtoDUNEBeamSpill/ProtoDUNEBeamSpill.h"
#include "dune/DuneObj/ProtoDUNEBeamSpill.h"
#include <bitset>
#include <iomanip>

#include "TTree.h"
#include "TPolyLine3D.h"

namespace proto {
  class BeamAna;
}


class proto::BeamAna : public art::EDAnalyzer {
public:
  explicit BeamAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamAna(BeamAna const &) = delete;
  BeamAna(BeamAna &&) = delete;
  BeamAna & operator = (BeamAna const &) = delete;
  BeamAna & operator = (BeamAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run const & r) override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p);
  std::bitset<sizeof(double)*CHAR_BIT> toBinary(const long num);  
  void matchTriggers(beamspill::ProtoDUNEBeamSpill spill);
  double GetPosition(size_t, size_t);

private:
  
  TTree * fOutTree;
  TPolyLine3D * theLine;

  // Declare member data here.
  //bool fLoadFromDB; // unused
  double  fTimeWindow;
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
};


proto::BeamAna::BeamAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  this->reconfigure(p);
}

void proto::BeamAna::analyze(art::Event const & e)
{

  std::cerr << "%%%%%%%%%% Getting ifbeam service handle %%%%%%%%%%" << std::endl;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifb;
  std::cerr << "%%%%%%%%%% Got ifbeam handle %%%%%%%%%%" << std::endl;

  std::unique_ptr<ifbeam_ns::BeamFolder> bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
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

  size_t nDev = fDevices.size();
  
  std::cout << "Using " << nDev << " devices: " << std::endl;
  for(size_t id = 0; id < fDevices.size(); ++id){
    std::cout << fPrefix + fDevices[id] << std::endl;
//    std::cout << "At: " << fCoordinates[id][0] << " " << fCoordinates[id][1] << " " << fCoordinates[id][2] << std::endl;
    std::cout << "Rotated: " << fRotations[id][0] << " " << fRotations[id][1] << " " << fRotations[id][2] << std::endl;
  }

  
 
  

  //Start getting spill information
  beamspill::ProtoDUNEBeamSpill * spill = new beamspill::ProtoDUNEBeamSpill();
  spill->InitFBMs(nDev);


  for(size_t it = 0; it < fMultipleTimes.size(); ++it){
    std::cout << "Time: " << fMultipleTimes[it] << std::endl;
    for(size_t d = 0; d < fDevices.size(); ++d){
      std::cout <<"Device: " << fDevices[d] << std::endl;
      std::vector<double> data = bfp->GetNamedVector(fMultipleTimes[it], fPrefix + fDevices[d] + ":eventsData[]");
      std::vector<double> counts = bfp->GetNamedVector(fMultipleTimes[it], fPrefix + fDevices[d] + ":countsRecords[]");
    
  
      std::cout << "Data: " << data.size() << std::endl;
      std::cout << "Counts: " << counts.size() << std::endl;
      for(size_t i = 0; i < counts.size(); ++i){
        std::cout << counts[i] << std::endl;
      }
      if(counts[1] > data.size())continue;

      beamspill::FBM fbm;
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
        spill->AddFBMTrigger(d,fbm);
        std::cout << std::endl;
      }

      for(size_t i = 0; i < spill->GetNFBMTriggers(d); ++i){
        spill->DecodeFibers(d,i);
        std::cout << fDevices[d] << " has active fibers: ";
        for(size_t iF = 0; iF < spill->GetActiveFibers(d,i).size(); ++iF)std::cout << spill->GetActiveFibers(d, i)[iF] << " "; 
        std::cout << std::endl;
      }

    }
  }

/*  for(size_t id = 0; id < 6; ++id){
    std::cout << "FBM: " << fDevices[id] << std::endl;
    std::cout << "N Triggers: " << spill->GetNFBMTriggers(id) << std::endl;
    
    for(size_t it = 0; it < spill->GetNFBMTriggers(id); ++it){
      std::cout << "Trigger: " << it <<std::endl;
      std::cout << "Hit Fibers: ";
      for(size_t iF = 0; iF < 192; ++iF){
        if( spill->GetFiberStatus(id,it,iF) ){
          std::cout << iF << " ";
        }
      }
      std::cout << std::endl;
    }
  }
*/

 //Setting some dummy Triggers for drawing  
/*
 double dummyTriggerTimeLSB[6] = {1.50000e+08,1.50000e+08,1.50002e+08,1.50002e+08,1.50004e+08,1.50004e+08};
 double dummyTriggerTimeMSB[6] = {1.53191e+09,1.53191e+09,1.53191e+09,1.53191e+09,1.53191e+09,1.53191e+09};
 double dummyEventTimeLSB[6]   = {1.50000e+08,1.50000e+08,1.50000e+08,1.50000e+08,1.50000e+08,1.50000e+08};
 double dummyEventTimeMSB[6]   = {1.53191e+09,1.53191e+09,1.53191e+09,1.53191e+09,1.53191e+09,1.53191e+09};

 double dummyHit[6][6]     = { {1,0,0,0,0,0},
                               {1,0,0,0,0,0},
                               {0,0,1,0,0,0},
                               {0,0,1,0,0,0},
                               {0,0,0,0,1,0},
                               {0,0,0,0,1,0}}; //first index: device, second LSB->MSB

 for(size_t id = 0; id < fDevices.size(); ++id){
   beamspill::FBM fbm;
   fbm.ID = id;

   fbm.timeData[0] = dummyTriggerTimeLSB[id];
   fbm.timeData[1] = dummyTriggerTimeMSB[id];
   fbm.timeData[2] = dummyEventTimeLSB[id];
   fbm.timeData[3] = dummyEventTimeMSB[id];

   for(int i = 0; i < fDevices.size(); ++i){
     fbm.fiberData[i] = dummyHit[id][i];
   }

   spill->AddFBMTrigger(id,fbm);
   spill->DecodeFibers(id,spill->GetNFBMTriggers(id) - 1);//Decode the last one
   for(size_t iF = 0; iF < 192; ++iF){
     if( spill->GetFiberStatus(id,spill->GetNFBMTriggers(id) - 1,iF) ){
       std::cout << iF << " ";
     }
   }
   std::cout << std::endl;
 }

 matchTriggers(*spill);
*/
}

void proto::BeamAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  fOutTree = tfs->make<TTree>("tree", "lines"); 
  //Need to make this configurable later
  theLine = new TPolyLine3D(3);
  fOutTree->Branch("Line", &theLine);
  // Implementation of optional member function here.
}

void proto::BeamAna::beginRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endJob()
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endSubRun(art::SubRun const & sr)
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

void proto::BeamAna::matchTriggers(beamspill::ProtoDUNEBeamSpill spill){
  
  std::vector<size_t> goodTriggers;
  bool foundNext = false;

  //Start with the earliest device. Should be most upstream
  //Go through each trigger and try to match
  for(size_t it = 0; it < spill.GetNFBMTriggers(0); ++it){

    double time = spill.DecodeFiberTime(0, it);
//    std::cout << time << std::endl;

    //Go through the next FBMs and see if they have a good time
    for(size_t id = 1; id < fDevices.size(); ++id){
      foundNext = false;

      for (size_t it2 = 0; it2 < spill.GetNFBMTriggers(id); ++it2){
//        std::cout << (spill.DecodeFiberTime(id, it2 ) - time) << std::endl;

        if ( ( (spill.DecodeFiberTime(id, it2 ) - time) < 0.00010e+08 ) && ( (spill.DecodeFiberTime(id, it2 ) - time) >= 0 ) ){

//          std::cout << "Found matching trigger in " << id << " " << spill.DecodeFiberTime(id, it2) << std::endl; 
          foundNext = true;
          if(goodTriggers.empty()) goodTriggers.push_back(it);
          goodTriggers.push_back(it2);

          //For now: break after first one
          break;
        }       
      }
      //Break after missing a device in the sequence
      if(!foundNext){
        goodTriggers.clear();
        break;
      }      
    }
    //Using found next from last one
    if(foundNext){
      std::cout << "Found good track" << std::endl;
        
      if(theLine->GetN()){
        //Reset it
        *theLine = TPolyLine3D(3);
      }

      std::vector<double> xPos;      
      std::vector<double> yPos;      
      for(size_t ig = 0; ig < goodTriggers.size(); ++ig){
        std::vector<short> active = spill.GetActiveFibers(ig, goodTriggers[ig]);
        std::cout << "Trigger " << ig << " has " << active.size() << " active fibers "<< active[0] << std::endl;        
        //just get the first
        double position = GetPosition(ig, active[0]);                
        if(!(ig % 2)){ //Even numbers
          xPos.push_back(position);
        }
        else{ //Odd numbers
          yPos.push_back(position);         
        }
      }
      std::cout << "Sizes: " << xPos.size() << " " << yPos.size() << std::endl;
      //Check if diff size?
      for(size_t ip = 0; ip < xPos.size(); ++ip){
//        theLine->SetPoint(ip, xPos.at(ip), yPos.at(ip), fCoordinates[2*ip][2]);
      }

      fOutTree->Fill();
    }
  }

  
}

double proto::BeamAna::GetPosition(size_t iDevice, size_t iFiber){
  if(iFiber > 192){ std::cout << "Please select fiber in range [0,191]" << std::endl; return -1.;}
  if(iDevice > fDevices.size() -1 ){ std::cout << "Device out of range " << std::endl; return -1.;}
  double size = fFiberDimension[iDevice];
  
  //Define 0th fiber as origin. Go to middle of the fiber
  double pos = size*iFiber + size/2.;
  return pos;
}

DEFINE_ART_MODULE(proto::BeamAna)
