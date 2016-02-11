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

#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"

//lbne-artdaq includes
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "lbne-raw-data/Overlays/SSPFragment.hh"
#include "lbne-raw-data/Overlays/anlTypes.hh"
#include "lbne-raw-data/Overlays/PennMilliSlice.hh"
#include "lbne-raw-data/Overlays/PennMicroSlice.hh"
#include "artdaq-core/Data/Fragments.hh"

//larsoft includes
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "Geometry/Geometry.h"

#include "tpcFragmentToRawDigits.h"
#include "PennToOffline.h"
#include "SSPReformatterAlgs.h"
#include "utilities/UnpackFragment.h"

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
  SSPReformatterAlgs sspReform;
  std::string fPTBFragType, fPTBRawDataLabel, fPTBOutputDataLabel, fPTBMapFile, fPTBMapDir;
  double fNOvAClockFrequency; //MHz
  bool fDebug;

  int Event = 0;
  std::map<int,int> fRCEChannelMap;
  std::map<int,int> fPTBMap;

  TTree* fTree;
  int EvNum;
  long long RCETime, SSPTime, PTBTime;
  long long RCE_PTB_diff, RCE_SSP_diff, SSP_PTB_diff;
  int NumADCs, ConsistRCE;
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
  ,sspReform           ( pset.get<fhicl::ParameterSet>("SSPReformatter"))
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

  if(fDebug) printParameterSet();

  BuildTPCChannelMap  (fRCEChannelMapFile, fRCEChannelMap);
  BuildPTBChannelMap(fPTBMapDir, fPTBMapFile, fPTBMap);
  
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
  ++Event;
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
    if(!RCErawFragments.isValid()){
      std::cerr << "Run: " << evt.run() << ", SubRun: " << evt.subRun()	<< ", Event: " << evt.event() << " is NOT VALID" << std::endl;
      throw cet::exception("RCErawFragments NOT VALID");
      return;
    }
    artdaq::Fragments const& rawFragmentsRCE = *RCErawFragments;
    for ( size_t fragIndex=0; fragIndex < rawFragmentsRCE.size(); ++fragIndex ) {
      int ThisADCcount = 0;
      const artdaq::Fragment &singleFragment = rawFragmentsRCE[fragIndex];
      lbne::TpcMilliSliceFragment millisliceFragment(singleFragment);
      auto numMicroSlices = millisliceFragment.microSliceCount();
      for (unsigned int i_micro=0;i_micro<numMicroSlices;i_micro++) { // Loop through all MicroSlices
	std::unique_ptr <const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(i_micro);
	auto numNanoSlices = microSlice->nanoSliceCount();
	if (numNanoSlices) {
	  ConsistRCE = 1;
	  ThisADCcount += numNanoSlices;
	} // If numNanoSlices.
	if ( fragIndex==0 && i_micro==0 ) {
	  //std::cout << "Getting RCE time for event " << evt.event();
	  if ( microSlice->nanoSliceCount() ) { // If this MicroSlice has any NanoSlices
	    RCETime = microSlice->nanoSlice(0)->nova_timestamp();
	    //std::cout << ", taking RCETime from first nanoslice " << RCETime << std::endl;
	  } // NanoSlice
	  else {
	    RCETime = microSlice->software_message();
	    //std::cout << ", taking RCETime from header " << RCETime << std::endl;
	  }
	} // Looking at first MicroSlice of first Fragment
      } // MicroSlice
      if ( fragIndex == 0 ) {
	NumADCs = ThisADCcount;
      } else {
	if ( ThisADCcount != NumADCs ) ConsistRCE = 0;
      }
      //std::cout << "Looking at event " << evt.event() << ", fragment " << fragIndex << ", it has " << ThisADCcount << ", InconsistRCEs? " << InconsistRCE << std::endl;
    } // rawFragments.size()
    //std::cout << "Got RCE start time, it is " << RCETime << std::endl;
  } //RCEPresent
  
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
    if(!SSPrawFragments.isValid()){
      mf::LogError("SSPToOffline") << "Run: " << evt.run() << ", SubRun: " << evt.subRun() << ", Event: " << evt.event() << " is NOT VALID";
      throw cet::exception("raw NOT VALID");
      return;
    }
    
    artdaq::Fragments const& rawFragmentsSSP = *SSPrawFragments;
    if ( rawFragmentsSSP.size() ) {
      const auto& frag(rawFragmentsSSP[0]);
      const SSPDAQ::MillisliceHeader* meta=0;
      if(frag.hasMetadata())
	{
	  meta = &(frag.metadata<lbne::SSPFragment::Metadata>()->sliceHeader);
	  
	  //std::cout << "=== SSP Metadata, Start time " << meta->startTime << ", End time " << meta->endTime << " Packet length " << meta->length << " Number of triggers " << meta->nTriggers << "===" << std::endl;
	  SSPTime = meta->startTime;
	}
    }
    //std::cout << "Got SSP start time, it is " << SSPTime << std::endl;
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
    if(!PTBrawFragments.isValid()){
      mf::LogError("PTBToOffline") << "Run: " << evt.run() << ", SubRun: " << evt.subRun() << ", Event: " << evt.event() << " is NOT VALID";
      throw cet::exception("raw NOT VALID");
      return;
    }

    lbne::PennMicroSlice::Payload_Header    *word_header    = nullptr;
    lbne::PennMicroSlice::Payload_Timestamp *FirstTimestamp = nullptr;
    uint32_t payload_index = 0;
    uint16_t counter, trigger, timestamp, payloadCount;
    uint8_t* payload_data = nullptr;
    
    if (PTBrawFragments->size()) {
      const auto& frag((*PTBrawFragments)[0]);
      lbne::PennMilliSliceFragment msf(frag);
      
      payloadCount = msf.payloadCount(counter, trigger, timestamp);
      
      while (payload_index < uint32_t(payloadCount-1) && FirstTimestamp == nullptr) {
	payload_data = msf.get_next_payload(payload_index,word_header);
	switch(word_header->data_packet_type) {
	case lbne::PennMicroSlice::DataTypeTimestamp:
	  FirstTimestamp = reinterpret_cast<lbne::PennMicroSlice::Payload_Timestamp*>(payload_data);
	  break;
	default:
	  break;
	}
      }
    }
    PTBTime = FirstTimestamp->nova_timestamp;
    //std::cout << "Got PTB start time, it is " << PTBTime << std::endl;
  } // PTB Present
  // ------------------------ NOW TO COMPARE ALL THESE NUMBERS -------------------------------------------
  
  if (Event%1000 == 0) std::cout << "Looking at event " << Event << " it had " << RCETime << " " << SSPTime << " " << PTBTime << " " << NumADCs << std::endl;

  RCE_PTB_diff = RCETime - PTBTime;
  RCE_SSP_diff = RCETime - SSPTime;
  SSP_PTB_diff = SSPTime - PTBTime;
  
  fTree->Fill();
}

void DAQToOffline::CheckTiming::endJob()
{
}

DEFINE_ART_MODULE(DAQToOffline::CheckTiming)
