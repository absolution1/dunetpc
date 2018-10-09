////////////////////////////////////////////////////////////////////////
// Class:       BeamDQM
// Plugin Type: analyzer (art v2_08_03)
// File:        BeamDQM_module.cc
//
// Generated at Thu Nov  2 22:57:41 2017 by Jonathan Paley using cetskelgen
// from cetlib version v3_01_01.
//
// Edited by Jake Calcutt
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
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNETimeStamp.h"
#include <bitset>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <limits>

#include "TTree.h"
#include "TH2F.h"
#include "TVectorD.h"
#include "TPolyMarker.h"

namespace proto {
  class BeamDQM;
}



class proto::BeamDQM : public art::EDAnalyzer {
public:
  explicit BeamDQM(fhicl::ParameterSet const & p);
  //  virtual ~BeamDQM();

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamDQM(BeamDQM const &) = delete;
  BeamDQM(BeamDQM &&) = delete;
  BeamDQM & operator = (BeamDQM const &) = delete;
  BeamDQM & operator = (BeamDQM &&) = delete;

  // Required functions.
  void reconfigure(fhicl::ParameterSet const & p);
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  
  uint64_t GetCTBInfo(art::Event const &);
  void     MomentumSpec(size_t);
  double   MomentumCosTheta(double, double, double);
  void     parseXCET(uint64_t);

  template<class T> 
  T FetchWithRetries(uint64_t, std::string, int);
   
private:
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;

  TTree * fOutTree;
  TH1F * fCKovHist;
  int HLTWord;
  long long int HLTTS;
  int BeamOn;
  int BITrigger;
  int C1;
  int C2;

  uint64_t prev_fetch_time;
  long long int prev_event_time;

  int eventNum;
  int runNum;
  int subRunNum;
  double CKov1Pressure;
  double CKov2Pressure;
  double CKov1Efficiency;
  double CKov2Efficiency;
  
  double  fTimeWindow;
  std::string fBundleName;
  std::string fURLStr;
  double fValidWindow;
  uint64_t fFixedTime;
  std::vector< uint64_t > fMultipleTimes;


  int fNRetries;

  // Names of the CERN Beam Devices
  // Names have the form of a prefix + device
  
  // Prefixes for different devices classes:
  // Beam Positions (BPF)
  // Time of flights (TOF)
  // and Cerenkov counters (CET)
  std::string fXCETPrefix;

  // Device Names for each device

  // Cerenkovs
  std::string fCKov1;
  std::string fCKov2;


  art::ServiceHandle<ifbeam_ns::IFBeam> ifb;
  uint64_t validTimeStamp;

  double L1=1.980, L2=1.69472, L3=2.11666;
  double magnetLen, magnetField;
  std::vector<double> current;

  //Hardware Parameters for magnetic field stuff
  double mag_P1 = 5.82044830e-3;
  double mag_P2 = 0.;
  double mag_P3 = -4.68880000e-6;
  double mag_P4 = 324.573967;

};

// Constructor
proto::BeamDQM::BeamDQM(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
  // Configure/Reconfigure
  this->reconfigure(p);
}
// END Constructor
////////////////////////

////////////////////////
// Fetch Method
template <class T> 
T proto::BeamDQM::FetchWithRetries(uint64_t time, std::string name, int nRetry){
  T theResult;
  
  std::cout << std::endl;
  uint64_t newTime;
  //Search at and above time given with nRetries
  //Will later search below, just in case the event time is actually greater than
  for(newTime = time; newTime < time + nRetry; ++newTime){
    std::cout << "Trying to grab from folder: " << name << std::endl;
    std::cout << "At Time: " << newTime << std::endl;    
    try{
      theResult = bfp->GetNamedVector(newTime, name);
      std::cout << "Successfully fetched" << std::endl;
      prev_fetch_time = newTime;
      return theResult;
    }
    catch(std::exception e){
      std::cout << "Could not fetch with time " << newTime << std::endl;      
    }
  }
  //Now search below
  for(newTime = time - 1; newTime > time - nRetry - 1; --newTime){
    std::cout << "Trying to grab from folder: " << name << std::endl;
    std::cout << "At Time: " << newTime << std::endl;    
    try{
      theResult = bfp->GetNamedVector(newTime, name);
      std::cout << "Successfully fetched" << std::endl;
      prev_fetch_time = newTime;
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
  std::cout << std::endl;
  return theResult; 
}
// END FetchWithRetries
////////////////////////


//Gets the Timing and CTB raw-decoded info.
//Finds the triggers, and looks for a valid trigger word
//(i.e. coming from beam)
//
//Returns the timestamp of the high level trigger.
uint64_t proto::BeamDQM::GetCTBInfo(art::Event const & e){

  std::cout << std::endl;
  std::cout << "Getting Raw Decoder Info" << std::endl;

  auto const CTBHandle = e.getValidHandle< std::vector< raw::ctb::pdspctb > >("ctbrawdecoder:daq");
  std::cout << "CTB valid? " << CTBHandle.isValid() << std::endl;
  if(CTBHandle.isValid()){
    auto const & CTB = (*CTBHandle)[0];

    bool noHLT = true;

    std::cout << "NTriggers: " << CTB.GetNTriggers() << std::endl;
    for (size_t nTrig = 0; nTrig < CTB.GetNTriggers(); ++nTrig){

      raw::ctb::Trigger ctbTrig = CTB.GetTrigger(nTrig);      
      uint32_t  theType  = ctbTrig.word_type;
      ULong64_t theWord  = ctbTrig.trigger_word;
      ULong64_t theTS    = ctbTrig.timestamp;

      if (theType == 2){
        std::cout << "Found the High Level Trigger" << std::endl;
        
        HLTWord = theWord;
        HLTTS = theTS;
        
        //The High Level Trigger consists of 8 bits
        //HLT7 -> HLT0
        std::bitset<8> theHLT(theWord);
        std::cout << "High Level Trigger: " << theHLT << std::endl;
        
        //HLT5 corresponds to excluding Low Level Triggers
        //from the Beamline
        //So return 0, we'll skip this event
        if (theHLT[5]) {         
          std::cout << "HLT 5 activated. Excluding beam events. Skipping this event" << std::endl;
          break;
        }
        else if (theHLT[0] && (theHLT.count() == 1)) {
          std::cout << "Only HLT 0 activated. This is just a random trigger. No beamline info was activated." << std::endl
                    << "Skipping this Event." << std::endl;
          break;
        }
        else{
          noHLT = false;
          std::cout << "Found valid beam event." << std::endl;
          break;
        }
      }       
    }

    if(noHLT){
      //This should never happen but...
      std::cout << "Error! No High Level Trigger Found!" << std::endl;
      return 0;
    }
    else{

      //Now check the channel statuses        
      std::cout << "ChStatuses: " << CTB.GetNChStatuses() << std::endl;
      raw::ctb::ChStatus theStatus = CTB.GetChStatuse(0);
      
      uint32_t the_beam_hi    = theStatus.beam_hi; 
      uint32_t the_beam_lo    = theStatus.beam_lo; 
      ULong64_t the_timestamp = theStatus.timestamp; 
      
      std::cout << "Timestamp : " << the_timestamp << std::endl;
      if( the_timestamp > ULong64_t(HLTTS) ){
        std::cout << "Found Channel status > HLT timestamp. Skipping" << std::endl;
//        continue;
      }

      std::cout << "beam_hi   : " << std::bitset<5>(the_beam_hi) << std::endl; 
      std::cout << "beam_lo   : " << std::bitset<4>(the_beam_lo) << std::endl; 

      std::bitset<4> beam_lo(the_beam_lo);
      std::bitset<5> beam_hi(the_beam_hi);

      BeamOn     =  beam_lo[0];
      BITrigger  =  beam_lo[1];
      C1         =  beam_lo[3];
      C2         =  beam_hi[0];


      std::cout << "%%%%Decoding the beam channels%%%" << std::endl;
      std::cout << "Beam On:    " << BeamOn    << std::endl
                << "BI Trigger: " << BITrigger << std::endl
                << "C1:         " << C1        << std::endl
                << "C2:         " << C2        << std::endl
                << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;        

      fCKovHist->Fill(C1 + 2*C2);

      return HLTTS;
    }
    
  }
  std::cout << "Error! Invalid CTB Handle!" << std::endl;
  std::cout << std::endl;
  return 0;

}



////////////////////////
// Analyzer Method (reads in the event and derives values)
void proto::BeamDQM::analyze(art::Event const & e)
{
  //Reset 
  HLTWord = 0;
  HLTTS = 0;
  BeamOn     = -1;
  BITrigger  = -1;
  C1         = -1;
  C2         = -1;

  eventNum = e.event();
  runNum = e.run();
  subRunNum = e.subRun();
   
  validTimeStamp = GetCTBInfo(e);
  std::cout << std::endl;
  if(validTimeStamp){

    



/*
    // Read in and cache the beam bundle folder for a specific time
    bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
    std::cerr << "%%%%%%%%%% Got beam folder %%%%%%%%%%" << std::endl;

    // Set the readout window of interest
    bfp->setValidWindow(fValidWindow);
    std::cerr << "%%%%%%%%%% Valid window " << bfp->getValidWindow() << " %%%%%%%%%%" << std::endl;

    std::cout << std::endl;
    // Now for the current event retrieve the time (trigger) of the event
    // This will take the form a 64-bit timestamp with a high and low word
    // Break this up to have access to each word and the long long word
    std::cout <<"Event Time: " << uint64_t(e.time().timeLow()) << std::endl;
    std::cout << "Low: " << e.time().timeLow()  << std::endl;
    std::cout << "High: " << e.time().timeHigh()  << std::endl;
    std::cout << std::endl;
    //Need to use either for different versions
    if( e.time().timeHigh() == 0 ) eventTime = e.time().timeLow();
    else eventTime = e.time().timeHigh();

    //Use multiple times provided to fcl
    if( fMultipleTimes.size() ){
      std::cout << "Using multiple times: " << std::endl;
      for(size_t it = 0; it < fMultipleTimes.size(); ++it){
        std::cout << fMultipleTimes[it] << std::endl;
      }
    }
    //Use singular time provided to fcl
    else if(fFixedTime > 0.){
      std::cout << "Using singular time: " << fFixedTime << std::endl;
      fMultipleTimes.push_back(fFixedTime);
    }
    //Use time of event
    else{
      std::cout <<" Using Event Time: " << uint64_t(eventTime) << std::endl;
      fMultipleTimes.push_back(uint64_t(eventTime));
      usedEventTime = true;
    } 
*/
    //Start getting beam event info
    // Loop over the different particle times 
/*    for(size_t it = 0; it < fMultipleTimes.size(); ++it){
      std::cout << "Time: " << fMultipleTimes[it] << std::endl;
      uint64_t fetch_time;
      if( abs( (long long)(fMultipleTimes[it]) - (long long)(prev_fetch_time) ) < 18 ){
        fetch_time = prev_fetch_time;
      }
      else{
        fetch_time = fMultipleTimes[it];
      }
   
      parseXCET(fetch_time);
*/    
 }
  

 // Write out the to tree
 fOutTree->Fill();
 current.clear();
// if(usedEventTime) fMultipleTimes.clear();
}
// END BeamDQM::analyze
////////////////////////


void proto::BeamDQM::parseXCET(uint64_t time){
  if(fCKov1 != ""){  
    std::cout << "Getting CKov1 info: " << fCKov1 << std::endl;

    std::vector< double > countsTrigCKov1 = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":countsTrig",fNRetries);
    std::vector< double > countsCKov1     = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":counts",fNRetries);
    std::vector< double > pressureCKov1   = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":pressure",fNRetries);
    std::cout << "countsTrig: " << countsTrigCKov1[0] << std::endl; 
    std::cout << "counts: "     << countsCKov1[0] << std::endl;
    std::cout << "pressure: "   << pressureCKov1[0] << std::endl;

    CKov1Pressure = pressureCKov1[0];
    if(countsTrigCKov1[0] > 0){
      CKov1Efficiency = countsCKov1[0]/countsTrigCKov1[0];
    }
  }
  
  if(fCKov2 != ""){
    std::cout << "Getting CKov2 info: " << fCKov2 << std::endl;

    std::vector< double > countsTrigCKov2 = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":countsTrig",fNRetries);
    std::vector< double > countsCKov2     = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":counts",fNRetries);
    std::vector< double > pressureCKov2   = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":pressure",fNRetries);
    std::cout << "countsTrig: " << countsTrigCKov2[0] << std::endl; 
    std::cout << "counts: "     << countsCKov2[0] << std::endl;
    std::cout << "pressure: "   << pressureCKov2[0] << std::endl;

    CKov2Pressure = pressureCKov2[0];
    if(countsTrigCKov2[0] > 0){
      CKov2Efficiency = countsCKov2[0]/countsTrigCKov2[0];
    }
  }

/*  beam::CKov CKov1Status, CKov2Status;

  double triggerTime = beamevt->GetT0( beamevt->GetActiveTrigger() ).first;
  triggerTime += 1.e-9*beamevt->GetT0( beamevt->GetActiveTrigger() ).second;

  CKov1Status.timeStamp = triggerTime;
  CKov1Status.pressure  = CKov1Pressure;
  CKov1Status.trigger   = C1;
  beamevt->SetCKov0( CKov1Status );

  CKov2Status.timeStamp = triggerTime;
  CKov2Status.pressure  = CKov2Pressure;
  CKov2Status.trigger   = C2;
  beamevt->SetCKov1( CKov2Status );
*/
}
// END BeamDQM::parseXCET
////////////////////////

////////////////////////
void proto::BeamDQM::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fCKovHist = tfs->make<TH1F>("CKov","",4,0,4);

  fOutTree = tfs->make<TTree>("tree", "lines"); 
  //eventTime = 0;
  eventNum = 0;
  runNum = 0;
  subRunNum = 0;
  CKov1Pressure = 0.;
  CKov2Pressure = 0.;
  CKov1Efficiency = 0.;
  CKov2Efficiency = 0.;
  fOutTree->Branch("Event", &eventNum);
  fOutTree->Branch("Run",   &runNum);
  fOutTree->Branch("Subrun", &subRunNum);
  fOutTree->Branch("Pressure1", &CKov1Pressure);
  fOutTree->Branch("Eff1", &CKov1Efficiency);
  fOutTree->Branch("Pressure2", &CKov2Pressure);
  fOutTree->Branch("Eff2", &CKov2Efficiency);
}

void proto::BeamDQM::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fBundleName  = p.get<std::string>("BundleName");
  fURLStr      = p.get<std::string>("URLStr");
  fNRetries    = p.get<int>("NRetries");

  fValidWindow = p.get<double>("ValidWindow");
  fTimeWindow  = p.get<double>("TimeWindow");

  fFixedTime   = p.get<uint64_t>("FixedTime");
  fMultipleTimes = p.get< std::vector<uint64_t> >("MultipleTimes");

  //For Momentum Spectrometry////
//  firstBPROF1  = p.get< std::string >("FirstBPROF1");
//  secondBPROF1 = p.get< std::string >("SecondBPROF1");
//  BPROF2       = p.get< std::string >("BPROF2");
//  BPROF3       = p.get< std::string >("BPROF3");
//  fBeamBend    = p.get< double >("BeamBend");
/*  L1           = 2.004;//(m)
  L2           = 1.718*cos(fBeamBend);//(m)
  L3           = 2.728*cos(fBeamBend);//(m)
*/
  magnetLen    = 1.;//(m)
  magnetField  = 1000.;//()
  ///////////////////////////////


  fCKov1 = p.get< std::string >("CKov1");
  fCKov2 = p.get< std::string >("CKov2");
  fXCETPrefix      = p.get<std::string>("XCETPrefix");
}

void proto::BeamDQM::MomentumSpec(size_t theTrigger){
  
/*  std::cout << "Doing momentum spectrometry for trigger " << beamevt->GetFullT0(theTrigger) << std::endl;

  //Get the active fibers from the upstream tracking XBPF
  std::string firstBPROF1Type    = fDeviceTypes[firstBPROF1]; 
  std::string secondBPROF1Type   = fDeviceTypes[secondBPROF1]; 
  std::vector<short> BPROF1Fibers;
  std::string BPROF1Name;

  if (firstBPROF1Type == "horiz" && secondBPROF1Type == "vert"){
    BPROF1Fibers = beamevt->GetActiveFibers(firstBPROF1, theTrigger);
    BPROF1Name = firstBPROF1;

    std::cout << firstBPROF1 << " has " << BPROF1Fibers.size() << " active fibers at time " << beamevt->GetFiberTime(firstBPROF1,theTrigger) << std::endl;
    for(size_t i = 0; i < BPROF1Fibers.size(); ++i){
      std::cout << BPROF1Fibers[i] << " ";
    }
    std::cout << std::endl;

  }
  else if(secondBPROF1Type == "horiz" && firstBPROF1Type == "vert"){
    BPROF1Fibers = beamevt->GetActiveFibers(secondBPROF1, theTrigger);
    BPROF1Name = secondBPROF1;

    std::cout << secondBPROF1 << " has " << BPROF1Fibers.size() << " active fibers at time " << beamevt->GetFiberTime(secondBPROF1,theTrigger) << std::endl;
    for(size_t i = 0; i < BPROF1Fibers.size(); ++i){
      std::cout << BPROF1Fibers[i] << " ";
    }
    std::cout << std::endl;
  }
  else{
    std::cout << "Error: type is not correct" << std::endl;
    return;
  }


  //////////////////////////////////////////////

  if( (BPROF1Fibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not checking momentum" << std::endl;
    return;
  }
  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  std::vector< short > strippedFibers; 
  std::cout << "BPROF1" << std::endl;
  for(size_t iF1 = 0; iF1 < BPROF1Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF1Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[iF1] == (BPROF1Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
  
  double X1, X2, X3;

  //X-direction convention is reversed for the spectrometer   
  X1 = -1.*GetPosition(BPROF1Name, strippedFibers[0])/1.E3;

  //BPROF2////
  //
  std::vector<short>  BPROF2Fibers = beamevt->GetActiveFibers(BPROF2, theTrigger);
  if( (BPROF2Fibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not checking momentum" << std::endl;
    return;
  }
  std::cout << BPROF2 << " has " << BPROF2Fibers.size() << " active fibers at time " << beamevt->GetFiberTime(BPROF2,theTrigger) << std::endl;
  for(size_t i = 0; i < BPROF2Fibers.size(); ++i){
    std::cout << BPROF2Fibers[i] << " ";
  }
  std::cout << std::endl;

  strippedFibers.clear(); 
  std::cout << "BPROF2" << std::endl;
  for(size_t iF1 = 0; iF1 < BPROF2Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF2Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF2Fibers.size() - 1){
      if (BPROF2Fibers[iF1] == (BPROF2Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  X2 = -1.*GetPosition(BPROF2, strippedFibers[0])/1.E3; 
  //X2 = X2*cos(TMath::Pi() - fBeamBend/1000.);
  ////////////

  //BPROF3////
  //
  std::vector<short>  BPROF3Fibers = beamevt->GetActiveFibers(BPROF3, theTrigger);
  if( (BPROF3Fibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not checking momentum" << std::endl;
    return;
  }
  std::cout << BPROF3 << " has " << BPROF3Fibers.size() << " active fibers at time " << beamevt->GetFiberTime(BPROF3,theTrigger) << std::endl;
  for(size_t i = 0; i < BPROF3Fibers.size(); ++i){
    std::cout << BPROF3Fibers[i] << " ";
  }
  std::cout << std::endl;

  strippedFibers.clear(); 
  std::cout << "BPROF3" << std::endl;
  for(size_t iF1 = 0; iF1 < BPROF3Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF3Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF3Fibers.size() - 1){
      if (BPROF3Fibers[iF1] == (BPROF3Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  X3 = -1.*GetPosition(BPROF3, strippedFibers[0])/1.E3; 
  //X3 = X3*cos(TMath::Pi() - fBeamBend/1000.);
  ////////////
  
  std::cout << "fBeamBend " << fBeamBend  << std::endl;
  double cosTheta = MomentumCosTheta(X1,X2,X3);
  std::cout <<"CosTheta " << cosTheta << std::endl;
  std::cout <<"Theta " << acos(cosTheta) << std::endl;
  double LB = mag_P1*fabs(current[0]);
  double deltaI = fabs(current[0]) - mag_P4;
  if(deltaI>0) LB+= mag_P3*deltaI*deltaI;

  double momentum = 299792458*LB/(1.E9 * acos(cosTheta));
  fMomentumHist->Fill(momentum);
  if( (BPROF1Fibers.size() == 1) && (BPROF2Fibers.size() == 1) && (BPROF3Fibers.size() == 1) ){
    std::cout << "Filling Cut Momentum Spectrum" << std::endl;
    fCutMomentum->Fill(momentum);

    beamevt->AddRecoBeamMomentum(momentum);
  }

  std::cout << "Momentum: " << 299792458*LB/(1.E9 * acos(cosTheta)) << std::endl; 


  std::cout << "Getting all trio-wise hits" << std::endl;
  std::cout << "N1,N2,N3 " << BPROF1Fibers.size()
            << " "         << BPROF2Fibers.size() 
            << " "         << BPROF3Fibers.size() << std::endl;
  for(size_t i1 = 0; i1 < BPROF1Fibers.size(); ++i1){
    double x1,x2,x3;

    x1 = -1.*GetPosition(BPROF1Name, BPROF1Fibers[i1])/1.E3;
    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Add .5 mm
        x1 += .0005;
      }
    }

    for(size_t i2 = 0; i2 < BPROF2Fibers.size(); ++i2){
      x2 = -1.*GetPosition(BPROF2, BPROF2Fibers[i2])/1.E3;
      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Add .5 mm
          x2 += .0005;
        }
      }

      for(size_t i3 = 0; i3 < BPROF3Fibers.size(); ++i3){
        std::cout << "\t" << i1 << " " << i2 << " " << i3 << std::endl;
        x3 = -1.*GetPosition(BPROF3, BPROF3Fibers[i3])/1.E3;
        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Add .5 mm
            x3 += .0005;
          }
        }

        double cosTheta_full = MomentumCosTheta(x1,x2,x3);        
        double momentum_full = 299792458*LB/(1.E9 * acos(cosTheta_full));
        fFullMomentum->Fill(momentum_full);

        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Skip the next
            ++i3;
          }
        }
        
      }

      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Skip the next
          ++i2;
        }
      }
    }

    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Skip the next
        ++i1;
      }
    }
  }
  */
}

double proto::BeamDQM::MomentumCosTheta(double X1, double X2, double X3){
/*
  double a = cos(fBeamBend)*(X3*L2 - X2*L3) / (L3 - L2);

  double cosTheta = (a + X1)*( (L3 - L2)*tan(fBeamBend) + (X2 - X3)*cos(fBeamBend) )+ (L3 - L2)*L1 ;

  double denomTerm1, denomTerm2, denom; 
  denomTerm1 = sqrt( L1*L1 + (a + X1)*(a + X1) );
  denomTerm2 = sqrt( 
               (L3 - L2)*(L3 - L2) 
             + TMath::Power( ( (L3 - L2)*tan(fBeamBend) 
                             + (X2 - X3)*cos(fBeamBend) ), 2 ) );
                            
  denom = denomTerm1*denomTerm2;
  cosTheta = cosTheta / denom;

  //
 
  double a = (X2*L3 - X3*L2) / (L3 - L2);
   
  double cosTheta = (a - X1)*(X3 - X2)*cos(fBeamBend) + (L3 - L2 )*( L1 + (a - X1)*tan(fBeamBend) );
    
  double denomTerm1, denomTerm2, denom; 
  denomTerm1 = sqrt( L1*L1 + (a - X1)*(a - X1) );
  denomTerm2 = sqrt( 
               (L3 - L2)*(L3 - L2) 
             + TMath::Power( ( (L3 - L2)*tan(fBeamBend) 
                             + (X3 - X2)*cos(fBeamBend) ), 2 ) );
  
  denom = denomTerm1*denomTerm2;
  cosTheta = cosTheta / denom;
 */
/*
    //double a = L2*tan(fBeamBend) + X2*cos(fBeamBend) - ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) )*( L2 - X2*sin(fBeamBend) );
    double a =  ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) )*( L2 - X2*sin(fBeamBend) );
    a = a / (L3 - X3*sin(fBeamBend) - L2 + X2*sin(fBeamBend) );
    a = L2*tan(fBeamBend) + X2*cos(fBeamBend) - a;
  
    double numTerm = (a - X1)*(L3 - L2)*tan(fBeamBend) + (a - X1)*(X3 - X2)*cos(fBeamBend) + L1*( (L3 - L2) - (X3 - X2)*sin(fBeamBend) );
  
    double denomTerm1, denomTerm2, denom;
    denomTerm1 = sqrt( L1*L1 + (a - X1)*(a - X1) );
    denomTerm2 = sqrt( TMath::Power( ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) ),2)
                     + TMath::Power( ( (L3 - L2)                - (X3 - X2)*sin(fBeamBend) ),2) );
    denom = denomTerm1 * denomTerm2; 

    double cosTheta = numTerm/denom;*/
   double cosTheta = 0.;
  return cosTheta;
}

DEFINE_ART_MODULE(proto::BeamDQM)
