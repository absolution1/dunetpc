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
  SSPReformatterAlgs     sspReform;
  std::string fPTBFragType, fPTBRawDataLabel, fPTBOutputDataLabel, fPTBChannelMapFile;
  double fNOvAClockFrequency; //MHz
  bool fUseChannelMap;
  bool fDebug;

  raw::Compress_t fCompression;   // compression type to use
  unsigned int    fZeroThreshold; // Zero suppression threshold

  std::map<int,int> fRCEChannelMap, fSSPChannelMap, fPTBChannelMap;
  std::pair <std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::CounterWordSize> >,
	     std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::TriggerWordSize> > > PrevTimeStampWords;

  TTree* fTree;
  int EvNum;
  int64_t RCETime, SSPTime, PTBTime;
  int NumADCs, AvNumADCs;
  bool InconsistRCE;
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
  ,fPTBChannelMapFile  ( pset.get<std::string>("PTBChannelMapFile"))
    //--------------------------------------------------------------------
  ,fNOvAClockFrequency ( pset.get<double>("NOvAClockFrequency")) // in MHz
  ,fUseChannelMap ( pset.get<bool>("UseChannelMap"))
  ,fDebug         ( pset.get<bool>("Debug"))
{
}

DAQToOffline::CheckTiming::~CheckTiming() {
}

void DAQToOffline::CheckTiming::beginJob() {

  fZeroThreshold=0;
  fCompression=raw::kNone;
  if(fDebug) printParameterSet();

  if (fUseChannelMap) {
    BuildTPCChannelMap  (fRCEChannelMapFile, fRCEChannelMap);
    BuildPTBChannelMap  (fPTBChannelMapFile, fPTBChannelMap);
  }

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("TimingCheck","TimingCheck");
  fTree->Branch("EvNum", &EvNum, "EvNum/I");
  fTree->Branch("RCETime", &RCETime, "RCETime/I");
  fTree->Branch("SSPTime", &SSPTime, "SSPTime/I");
  fTree->Branch("PTBTime", &PTBTime, "PTBTime/I");
  fTree->Branch("NumADCs", &NumADCs, "NumADCs/I");
  fTree->Branch("AvNumADCs", &AvNumADCs, "AvNumADCs/I");
  fTree->Branch("InconsistRCE", &InconsistRCE, "InconsistRCE/B");
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

  RCETime = 0, SSPTime = 0, PTBTime = 0;
  InconsistRCE = false;
  NumADCs = AvNumADCs = 0;
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
	std::unique_ptr <const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(0);
	auto numNanoSlices = microSlice->nanoSliceCount();
	ThisADCcount += numNanoSlices;
	if ( fragIndex==0 && i_micro==0 ) {
	  if ( microSlice->nanoSliceCount() ) { // If this MicroSlice has any NanoSlices
	    RCETime = microSlice->nanoSlice(0)->nova_timestamp();
	    //std::cout << "Taking RCETime from first nanoslice." << std::endl;
	  } // NanoSlice
	  else {
	    RCETime = microSlice->software_message();
	    //std::cout << "Taking RCETime from header." << std::endl;
	  }
	} // Looking at first MicroSlice of first Fragment
      } // MicroSlice
      if ( fragIndex == 0 ) {
	NumADCs = ThisADCcount;
      } else {
	if ( ThisADCcount != NumADCs ) InconsistRCE = true;
      }
      AvNumADCs += ThisADCcount;
    } // rawFragments.size()
    AvNumADCs = AvNumADCs / rawFragmentsRCE.size();
    std::cout << "Got RCE start time, it is " << RCETime << std::endl;
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
    std::cout << "Got SSP start time, it is " << SSPTime << std::endl;
    std::vector<raw::OpDetWaveform> waveforms = sspReform.SSPFragmentToOpDetWaveform(*SSPrawFragments);
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
    
    auto triggers = PennFragmentToExternalTrigger(*PTBrawFragments, fPTBChannelMap, PrevTimeStampWords);
    
  } // PTB Present
  // ------------------------ NOW TO COMPARE ALL THESE NUMBERS -------------------------------------------
  
  std::cout << RCETime << " " << SSPTime << " " << PTBTime << " " << NumADCs << std::endl;
  fTree->Fill();
}

void DAQToOffline::CheckTiming::endJob()
{
}

DEFINE_ART_MODULE(DAQToOffline::CheckTiming)
