////////////////////////////////////////////////////////////////////////
// Class:       CheckTiming
// Module Type: Analyzer
// File:        CheckTiming_module.cc
//
// Module written by Karl Warburton to compare the timings for RCE, SSP
//  and PTB data streams. 
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "utilities/UnpackFragment.h"

#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "lbne-raw-data/Overlays/SSPFragment.hh"
#include "lbne-raw-data/Overlays/anlTypes.hh"
#include "lbne-raw-data/Overlays/PennMilliSlice.hh"
#include "lbne-raw-data/Overlays/PennMicroSlice.hh"
#include "lbne-raw-data/Overlays/PennMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragments.hh"

#include "CheckTime.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"

namespace DAQToOffline {
  class CheckTiming;
}

class DAQToOffline::CheckTiming : public art::EDAnalyzer {
public:
  explicit CheckTiming(fhicl::ParameterSet const & p);
  virtual ~CheckTiming();

  void analyze(art::Event const & evt) override;
  void beginJob();
  void endJob();
  void printParameterSet();

private:


  
  std::string fRCEFragType, fRCERawDataLabel, fRCEOutputDataLabel, fRCEChannelMapFile;
  std::string fSSPFragType, fSSPRawDataLabel, fSSPOutputDataLabel;
  std::string fPTBFragType, fPTBRawDataLabel, fPTBOutputDataLabel, fPTBMapFile, fPTBMapDir;
  double fNOvAClockFrequency; //MHz
  bool fDebug;

  std::map<int,int> fRCEChannelMap;
  std::map<int,int> fPTBMap;

  TTree* fTree;
  int EvNum;
  long long RCETime, SSPTime, PTBTime;
  long long RCE_PTB_diff, RCE_SSP_diff, SSP_PTB_diff;
  int NumADCs, ConsistRCE;

  int RunNumber, nevts, nSSP, nRCE, nPTB, nSSPPayloads, nRCEPayloads, nConsistRCEPayloads, nPTBPayloads;
};

DAQToOffline::CheckTiming::CheckTiming(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
    //---------------------------RCE--------------------------------------
  ,fRCEFragType        ( pset.get<std::string>("RCEFragType"))
  ,fRCERawDataLabel    ( pset.get<std::string>("RCERawDataLabel"))
  ,fRCEOutputDataLabel ( pset.get<std::string>("RCEOutputDataLabel"))
  ,fRCEChannelMapFile  ( pset.get<std::string>("RCEChannelMapFile"))
    //---------------------------SSP--------------------------------------
  ,fSSPFragType        ( pset.get<std::string>("SSPFragType"))
  ,fSSPRawDataLabel    ( pset.get<std::string>("SSPRawDataLabel"))
  ,fSSPOutputDataLabel ( pset.get<std::string>("SSPOutputDataLabel"))
    //---------------------------PTB--------------------------------------
  ,fPTBFragType        ( pset.get<std::string>("PTBFragType"))
  ,fPTBRawDataLabel    ( pset.get<std::string>("PTBRawDataLabel"))
  ,fPTBOutputDataLabel ( pset.get<std::string>("PTBOutputDataLabel"))
  ,fPTBMapFile         ( pset.get<std::string>("PTBMapFile"))
  ,fPTBMapDir          ( pset.get<std::string>("PTBMapDir"))
    //--------------------------------------------------------------------
  ,fNOvAClockFrequency ( pset.get<double>("NOvAClockFrequency")) // in MHz
  ,fDebug         ( pset.get<bool>("Debug"))
{
}

DAQToOffline::CheckTiming::~CheckTiming() {
}

void DAQToOffline::CheckTiming::beginJob() {

  RunNumber = nevts = nSSP = nRCE = nPTB = nSSPPayloads = nRCEPayloads = nConsistRCEPayloads = nPTBPayloads = 0;

  if(fDebug) printParameterSet();
  
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("TimingCheck","TimingCheck");
  fTree->Branch("EvNum"       , &EvNum       , "EvNum/I");
  fTree->Branch("RCETime"     , &RCETime     , "RCETime/I");
  fTree->Branch("SSPTime"     , &SSPTime     , "SSPTime/I");
  fTree->Branch("PTBTime"     , &PTBTime     , "PTBTime/I");
  fTree->Branch("RCE_SSP_diff", &RCE_SSP_diff,"RCE_SSP_diff/I");
  fTree->Branch("RCE_PTB_diff", &RCE_PTB_diff,"RCE_PTB_diff/I");
  fTree->Branch("SSP_PTB_diff", &SSP_PTB_diff,"SSP_PTB_diff/I");
  fTree->Branch("NumADCs"     , &NumADCs     , "NumADCs/I");
  fTree->Branch("ConsistRCE"  , &ConsistRCE  , "ConsistRCE/I");
}

void DAQToOffline::CheckTiming::printParameterSet(){

  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
  std::cout << "RCE's" << std::endl;
  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
  std::cout << "fFragType: " << fRCEFragType << std::endl;
  std::cout << "fRawDataLabel: " << fRCERawDataLabel << std::endl;
  std::cout << "fOutputDataLabel: " << fRCEOutputDataLabel << std::endl;
  std::cout << "fDebug: ";
  if(fDebug) std::cout << "true" << std::endl;
  else std::cout << "false" << std::endl;

  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
  std::cout << "SSP's" << std::endl;
  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
  std::cout << "fFragType: " << fSSPFragType << std::endl;
  std::cout << "fRawDataLabel: " << fSSPRawDataLabel << std::endl;
  std::cout << "fOutputDataLabel: " << fSSPOutputDataLabel << std::endl;

  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
  std::cout << "PTB's" << std::endl;
  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
  std::cout << "fFragType: " << fPTBFragType << std::endl;
  std::cout << "fRawDataLabel: " << fPTBRawDataLabel << std::endl;
  std::cout << "fOutputDataLabel: " << fPTBOutputDataLabel << std::endl;

  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
}

void DAQToOffline::CheckTiming::analyze(art::Event const & evt)
{
  if (nevts == 0 ) {
    RunNumber = evt.run();
    std::cout << "Got runNumber it is " << RunNumber << std::endl;
  }

  ++nevts;
  RCETime = 0, SSPTime = 0, PTBTime = 0;
  ConsistRCE = -1;
  NumADCs = 0;
  EvNum = evt.event();
  // ------------- RCE Section ------------------------
  bool RCEPresent = true;
  art::Handle<artdaq::Fragments> RCErawFragments;
  evt.getByLabel(fRCERawDataLabel, fRCEFragType, RCErawFragments);

  try { RCErawFragments->size(); }
  catch(std::exception e) {
    std::cout << "WARNING: Raw RCE data not found in event " << evt.event() << std::endl;
    RCEPresent = false;
  }

  if (RCEPresent) {
    ++nRCE;
    if(!RCErawFragments.isValid()){
      std::cerr << "Run: " << evt.run() << ", SubRun: " << evt.subRun()	<< ", Event: " << evt.event() << " is NOT VALID" << std::endl;
      throw cet::exception("RCErawFragments NOT VALID");
      return;
    }
    artdaq::Fragments const& rawFragmentsRCE = *RCErawFragments;
    DAQToOffline::GetRCEFirstTimestamp ( rawFragmentsRCE, ConsistRCE, NumADCs, RCETime );
  } //RCEPresent
  std::cout << "Got RCE start time, it is " << RCETime << std::endl;
  if (ConsistRCE > -1 ) ++nRCEPayloads;
  if (ConsistRCE == 1 ) ++nConsistRCEPayloads;

  // ------------- SSP Section ------------------------
  bool SSPPresent = true;
  art::Handle<artdaq::Fragments> SSPrawFragments;
  evt.getByLabel(fSSPRawDataLabel, fSSPFragType, SSPrawFragments);

  try { SSPrawFragments->size(); }
  catch(std::exception e) {
    mf::LogWarning("SSPToOffline") << "WARNING: Raw SSP data not found in event " << evt.event();
    SSPPresent = false;
  }

  if (SSPPresent) {
    ++nSSP;
    if(!SSPrawFragments.isValid()){
      mf::LogError("SSPToOffline") << "Run: " << evt.run() << ", SubRun: " << evt.subRun() << ", Event: " << evt.event() << " is NOT VALID";
      throw cet::exception("raw NOT VALID");
      return;
    }
    
    artdaq::Fragments const& rawFragmentsSSP = *SSPrawFragments;
    DAQToOffline::GetSSPFirstTimestamp ( rawFragmentsSSP, nSSPPayloads, SSPTime );
    std::cout << "Got SSP start time, it is " << SSPTime << std::endl;
  } // SSP Present
  
  // ------------- PTB Section ------------------------
  bool PTBPresent = true;
  art::Handle<artdaq::Fragments> PTBrawFragments;
  evt.getByLabel(fPTBRawDataLabel, fPTBFragType, PTBrawFragments);
  
  try { PTBrawFragments->size(); }
  catch(std::exception e) {
    mf::LogWarning("PTBToOffline") << "WARNING: Raw PTB data not found in event " << evt.event();
    PTBPresent = false;
  }

  if (PTBPresent) {
    ++nPTB;
    if(!PTBrawFragments.isValid()){
      mf::LogError("PTBToOffline") << "Run: " << evt.run() << ", SubRun: " << evt.subRun() << ", Event: " << evt.event() << " is NOT VALID";
      throw cet::exception("raw NOT VALID");
      return;
    }

    artdaq::Fragments const& rawFragmentsPTB = *PTBrawFragments;
    DAQToOffline::GetPTBFirstTimestamp( rawFragmentsPTB, nPTBPayloads, PTBTime );
    std::cout << "Got PTB start time, it is " << PTBTime << std::endl;
  } // PTB Present
  // ------------------------ NOW TO COMPARE ALL THESE NUMBERS -------------------------------------------
  
  if (nevts%1000 == 0) 
    std::cout << "Looking at event " << evt.event() << " it had " << RCETime << " " << SSPTime << " " << PTBTime << " " << NumADCs << std::endl;

  RCE_PTB_diff = RCETime - PTBTime;
  RCE_SSP_diff = RCETime - SSPTime;
  SSP_PTB_diff = SSPTime - PTBTime;
  
  fTree->Fill();
}

void DAQToOffline::CheckTiming::endJob()
{
  std::cout << "IDENTIFIER: Run: " << RunNumber << " has " << nevts << " events in total, " 
	    << nSSP << " have SSP data and " << nSSPPayloads << " payloads, " 
	    << nRCE << " have RCE data and " << nRCEPayloads << " payloads " << nConsistRCEPayloads << " of which had consistent numbers of RCE payloads, "
	    << nPTB << " have PTB data and " << nPTBPayloads << " payloads." << std::endl;
}

DEFINE_ART_MODULE(DAQToOffline::CheckTiming)
