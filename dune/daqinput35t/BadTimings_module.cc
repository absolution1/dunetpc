//////////////////////////////////////////////////////////////////////////
// Class:       BadTimings
// Module type: analyser
// File:        BadTimings_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), August 2016
//
// Short analyser to find any timing disparities between detector
// components. Written because of the observed offset between RCEs and
// PTB in 35ton.
//////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// lbne-raw-data includes
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "lbne-raw-data/Overlays/PennMilliSlice.hh"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "tpcFragmentToRawDigits.h"
#include "utilities/UnpackFragment.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "PennToOffline.h"

#include "TH1I.h"
#include "TTree.h"

namespace DAQToOffline {
  class BadTimings;
}

class DAQToOffline::BadTimings : public art::EDAnalyzer {

public:

  explicit BadTimings(const fhicl::ParameterSet& pset);

  void analyze(const art::Event& evt);

private:

  std::string fTPCLabel, fTPCInstance;
  std::string fPTBLabel, fPTBInstance;

  lbne::TpcNanoSlice::Header::nova_timestamp_t fRCETriggerTimestamp;
  lbne::PennMilliSlice::Header::timestamp_t fPTBTriggerTimestamp;
  // unsigned long long fRCETriggerTimestamp;
  // unsigned long long fPTBTriggerTimestamp;

  art::ServiceHandle<geo::Geometry> fGeometry;
  art::ServiceHandle<art::TFileService> fTFS;

  std::string fPTBMapFile;
  std::map<int,int> fPTBMap;

  std::vector<bool> fWithinTrigger;

  TTree* fTree;
  int fRun, fEvent, fRCE;
  int fDiffTriggerTimestamps;
  int fMicrosliceWithTrigger;

  TH1I* hDiffTriggerTimestamps;
  TH1I* hTriggerStart;

};

DAQToOffline::BadTimings::BadTimings(const fhicl::ParameterSet& pset) : EDAnalyzer(pset) {
  fTPCLabel = pset.get<std::string>("TPCLabel");
  fTPCInstance = pset.get<std::string>("TPCInstance");
  fPTBLabel = pset.get<std::string>("PTBLabel");
  fPTBInstance = pset.get<std::string>("PTBInstance");
  fPTBMapFile = pset.get<std::string>("PTBMapFile");

  fWithinTrigger = std::vector<bool>(16,false);

  hDiffTriggerTimestamps = fTFS->make<TH1I>("PTBRCEDiffTimestamps",";RCE Trigger Timestamp - PTB Trigger Timestamp (NOvA ticks);",100,0,2000);
  hTriggerStart = fTFS->make<TH1I>("TriggerStart",";Microslice number containing trigger;",15,0,15);
  hTriggerStart->GetXaxis()->SetNdivisions(15);
  hTriggerStart->GetXaxis()->CenterLabels();

  fTree = fTFS->make<TTree>("BadTimings","BadTimings");
  fTree->Branch("Run",                  &fRun);
  fTree->Branch("Event",                &fEvent);
  fTree->Branch("RCE",                  &fRCE);
  // fTree->Branch("RCETriggerTimestamp",  &fRCETriggerTimestamp);
  // fTree->Branch("PTBTriggerTimestamp",  &fPTBTriggerTimestamp);
  fTree->Branch("DiffTriggerTimestamps",&fDiffTriggerTimestamps);
  fTree->Branch("MicrosliceWithTrigger",&fMicrosliceWithTrigger);

  DAQToOffline::BuildPTBChannelMap("", fPTBMapFile, fPTBMap);
}

void DAQToOffline::BadTimings::analyze(const art::Event& evt) {

  fRun = evt.run();
  fEvent = evt.event();

  // Get the TPC data out of the event
  art::Handle<artdaq::Fragments> RCERawFragments;
  evt.getByLabel(fTPCLabel, fTPCInstance, RCERawFragments);

  // Get the counter data out of the event
  art::Handle<artdaq::Fragments> PTBRawFragments;
  evt.getByLabel(fPTBLabel, fPTBInstance, PTBRawFragments);

  // Look at PTB
  for (size_t fragIndex = 0; fragIndex < PTBRawFragments->size(); ++fragIndex) {
    const artdaq::Fragment &singleFragment = (*PTBRawFragments)[fragIndex];
    lbne::PennMilliSliceFragment msf(singleFragment);

    lbne::PennMicroSlice::Payload_Header *word_header = nullptr;
    //lbne::PennMicroSlice::Payload_Counter *word_p_counter = nullptr;
    lbne::PennMicroSlice::Payload_Trigger *word_p_trigger = nullptr;
    lbne::PennMicroSlice::Payload_Timestamp *previous_timestamp = nullptr;
    lbne::PennMicroSlice::Payload_Header *future_timestamp_header = nullptr;
    lbne::PennMicroSlice::Payload_Timestamp *future_timestamp = nullptr;
    uint8_t* payload_data = nullptr;
    uint32_t payload_index = 0;

    uint16_t counter, trigger, timestamp, payloadCount;
    payloadCount = msf.payloadCount(counter, trigger, timestamp);

    if (trigger != 1) {
      //std::cout << "Event has no trigger in PTB data" << std::endl;
      return;
    }

    while (payload_index < uint32_t(payloadCount-1)) {

      payload_data = msf.get_next_payload(payload_index,word_header);
      if (payload_data == nullptr)
        continue;

      if (word_header->data_packet_type == lbne::PennMicroSlice::DataTypeTimestamp)
	previous_timestamp = reinterpret_cast<lbne::PennMicroSlice::Payload_Timestamp*>(payload_data);
      if (word_header->data_packet_type != lbne::PennMicroSlice::DataTypeTrigger)
	continue;

      word_p_trigger = reinterpret_cast<lbne::PennMicroSlice::Payload_Trigger*>(payload_data);
      if (!word_p_trigger->has_muon_trigger()) {
	std::cout << "Trigger is not from muon" << std::endl;
	return;
      }

      if (previous_timestamp != nullptr)
	fPTBTriggerTimestamp = word_header->get_full_timestamp_post(previous_timestamp->nova_timestamp);
      else {
	if (future_timestamp == nullptr) {
	  future_timestamp = reinterpret_cast<lbne::PennMicroSlice::Payload_Timestamp*>(msf.get_next_timestamp(future_timestamp_header));
	  if (future_timestamp == nullptr) {
	    std::cerr << "CAN'T FIND PTB TIMESTAMP WORDS IN MILLISLICE FRAGMENT!!! Logic will fail." << std::endl;
	    return;
	  }
	}
	fPTBTriggerTimestamp = word_header->get_full_timestamp_pre(future_timestamp->nova_timestamp);
      }

    } // loop over payloads

  } // fragments

  // Look at RCE
  // Create a map containing (fragmentID, fragIndex) for the event, will be used to check if each channel is present
  std::map<unsigned int,unsigned int> mapFragID;
  for (size_t fragIndex = 0; fragIndex < RCERawFragments->size(); ++fragIndex) {
    const artdaq::Fragment &singleFragment = (*RCERawFragments)[fragIndex];
    unsigned int fragmentID = singleFragment.fragmentID();
    mapFragID.insert(std::pair<unsigned int, unsigned int>(fragmentID,fragIndex));
  }

  int maxRCEs = 16;

  // Loop over the RCEs
  for (std::map<unsigned int,unsigned int>::iterator rceIt = mapFragID.begin(); rceIt != mapFragID.end() and std::distance(mapFragID.begin(),rceIt) < maxRCEs; ++rceIt) {

    fRCE = rceIt->first - 100;

    // Get millislice
    const artdaq::Fragment &singleFragment = (*RCERawFragments)[rceIt->second];
    lbne::TpcMilliSliceFragment millisliceFragment(singleFragment);
    //uint64_t millisliceTimestamp = millisliceFragment.microSlice(0)->software_message() & (((long int)0x1 << 56) - 1);

    // Get microslices
    auto numMicroSlices = millisliceFragment.microSliceCount();
    for (unsigned int i_micro = 0; i_micro < numMicroSlices; i_micro++) {

      std::unique_ptr<const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(i_micro);
      uint64_t rceTriggerTimestamp = microSlice->software_message() & (((long int)0x1 << 56) - 1);

      // Get nanoslices
      auto numNanoSlices = microSlice->nanoSliceCount();

      if (!fWithinTrigger[fRCE] and numNanoSlices) {

	// New trigger
	fWithinTrigger[fRCE] = true;
	std::cout << std::endl << "New trigger (RCE " << fRCE << ") on microslice " << i_micro << " (event " << evt.event() << ")" << std::endl;
	//uint64_t nanoSliceTimestamp = microSlice->nanoSlice(0)->nova_timestamp();

	fRCETriggerTimestamp = rceTriggerTimestamp;
	fMicrosliceWithTrigger = (rceTriggerTimestamp - microSlice->nanoSlice(0)->nova_timestamp())/*time*/ * (15.625/500)/*nova->tpc tick*/ * 1/1000/* nano->micro*/;
	fDiffTriggerTimestamps = fRCETriggerTimestamp - fPTBTriggerTimestamp;

	hDiffTriggerTimestamps->Fill(fDiffTriggerTimestamps);
	hTriggerStart->Fill(fMicrosliceWithTrigger);

	fTree->Fill();

      }
      else if (fWithinTrigger[fRCE] and !numNanoSlices) {

	// End of trigger
	fWithinTrigger[fRCE] = false;
	std::cout << "End of trigger (RCE " << fRCE << ") on microslice " << i_micro << " (event " << evt.event() << ")" << std::endl;

      }

    } // microslices

  } // rces

  return;

}

DEFINE_ART_MODULE(DAQToOffline::BadTimings)




  // // Create a map containing (fragmentID, fragIndex) for the event, will be used to check if each channel is present
  // std::map<unsigned int,unsigned int> mapFragID;
  // for (size_t fragIndex = 0; fragIndex < RCERawFragments->size(); ++fragIndex) {
  //   const artdaq::Fragment &singleFragment = (*RCERawFragments)[fragIndex];
  //   unsigned int fragmentID = singleFragment.fragmentID();
  //   mapFragID.insert(std::pair<unsigned int, unsigned int>(fragmentID,fragIndex));
  // }

  // fTimestampSet = false;

  // for (size_t chan = 0; chan < fGeometry->Nchannels(); ++chan) {

  //   // Each channel is uniquely identified by (fragmentID, group, sample) in an online event
  //   unsigned int fragmentID = UnpackFragment::getFragIDForChan(chan);
  //   //unsigned int sample = UnpackFragment::getNanoSliceSampleForChan(chan);

  //   if (mapFragID.find(fragmentID) == mapFragID.end())
  //     continue;

  //   unsigned int fragIndex = mapFragID[fragmentID];

  //   // Get millislice
  //   const artdaq::Fragment &singleFragment = (*RCERawFragments)[fragIndex];
  //   lbne::TpcMilliSliceFragment millisliceFragment(singleFragment);

  //   // Get microslices
  //   auto numMicroSlices = millisliceFragment.microSliceCount();
  //   for (unsigned int i_micro = 0; i_micro < numMicroSlices; i_micro++) {

  //     std::unique_ptr<const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(i_micro);
  //     uint64_t novaTimestamp = microSlice->software_message() & (((long int)0x1 << 56) - 1);

  //     // Get nanoslices
  //     auto numNanoSlices = microSlice->nanoSliceCount();

  //     // Only interested here in the first microslice in a trigger
  //     if (fragIndex == 0 and i_micro == 0 and numNanoSlices) {
  //     }

  //     // Get timestamps for each RCE
  //     if (chan % 128 == 0) {
  // 	if (!fWithinTrigger && numNanoSlices) {
  // 	  fWithinTrigger = true;
  // 	  uint64_t nanoSliceTimestamp = microSlice->nanoSlice(0)->nova_timestamp();
  // 	  //std::cout << "Timestamps: nova " << novaTimestamp << " and first nanoslice " << nanoSliceTimestamp << std::endl;
  // 	  fRCENanosliceTimestamp = nanoSliceTimestamp;
  // 	  fRCETriggerTimestamp = novaTimestamp;
  // 	  if (!fTimestampSet || nanoSliceTimestamp < fRCENanosliceTimestamp) {
  // 	    fRCENanosliceTimestamp = nanoSliceTimestamp;
  // 	    fTimestampSet = true;
  // 	  }
  // 	  std::cout << "Fragment " << fragmentID << ", microslice " << i_micro << ": " << std::endl << "First nanoslice: " << fRCENanosliceTimestamp << ", RCE trigger " << fRCETriggerTimestamp << " and PTB trigger " << fPTBTriggerTimestamp << " (final timestamp is " << fRCENanosliceTimestamp << ")" << std::endl;
  // 	}
  // 	else if (fWithinTrigger && !numNanoSlices) {
  // 	  fWithinTrigger = false;
  // 	}
  //     }

  //   } // microslices

  // } // channels
