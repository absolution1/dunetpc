////////////////////////////////////////////////////////////////////////
// File:   CheckTime.cc
// Author: Karl Warburton (Feb 2016)
//
// Utility to provide methods for checking the timestamps of the 
// DAQ components.
////////////////////////////////////////////////////////////////////////

#include "CheckTime.h"
#include <map>
#include <iostream>
#include <iomanip>
//=======================================================================================
void DAQToOffline::GetRCEFirstTimestamp( artdaq::Fragments const& Fragments, int &ConsistRCE, int &NumADCs, long long &RCETime ) {
  
  for ( size_t fragIndex=0; fragIndex < Fragments.size(); ++fragIndex ) {
    int ThisADCcount = 0;
    const artdaq::Fragment &singleFragment = Fragments[fragIndex];
    lbne::TpcMilliSliceFragment millisliceFragment(singleFragment);
    auto numMicroSlices = millisliceFragment.microSliceCount();
    for (unsigned int i_micro=0;i_micro<numMicroSlices;i_micro++) { // Loop through all MicroSlices
      std::unique_ptr <const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(i_micro);
      auto numNanoSlices = microSlice->nanoSliceCount();
      if (numNanoSlices) {
	ConsistRCE = 1;
	ThisADCcount += numNanoSlices;
	//if (ThisADCcount == 1000 )
	//std::cout << "Some microslice things..."
	//	    << " Soft Message " << microSlice->software_message()
	//	    << ", Firm Message " << microSlice->firmware_message()
	//	    << ", Type ID " << microSlice->type_id()
	//	    << ", sequence ID " << microSlice->sequence_id()
	//	    << ", SoftTrig? " << microSlice->softTrig()
	//	    << ", Ext Trig? " << microSlice->extTrig()  << std::endl;
      } // If numNanoSlices.
      if ( fragIndex==0 && i_micro==0 ) {
	// MWallbank 8/30/16: fixed this bug. The 8 msb need to be masked as they represent the external trigger counter
	uint64_t timestamp = microSlice->software_message() & (((long int)0x1 << 56) - 1);
	if ( microSlice->nanoSliceCount() ) { // If this MicroSlice has any NanoSlices
	  RCETime = microSlice->nanoSlice(0)->nova_timestamp();
	} // NanoSlice
	else {
	  RCETime = timestamp;
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
  return;
}
//=======================================================================================
void DAQToOffline::GetSSPFirstTimestamp( artdaq::Fragments const& Fragments, int &nSSPPayloads, long long &SSPTime ) {
  if ( Fragments.size() ) {
    const auto& frag(Fragments[0]);
    const SSPDAQ::MillisliceHeader* meta=0;
    if(frag.hasMetadata()) {
      ++nSSPPayloads;
      
      meta = &(frag.metadata<lbne::SSPFragment::Metadata>()->sliceHeader);
      //std::cout << "=== SSP Metadata, Start time " << meta->startTime << ", End time " << meta->endTime << " Packet length " << meta->length << " Number of triggers " << meta->nTriggers << "===" << std::endl;
      SSPTime = meta->startTime;
    }
  }
  //std::cout << "Got SSP start time, it is " << SSPTime << std::endl;
  return;
}
//=======================================================================================
void DAQToOffline::GetPTBFirstTimestamp( artdaq::Fragments const& PTBrawFragments, int &nPTBPayloads, long long &PTBTime ) {

  if (PTBrawFragments.size()) {
    ++nPTBPayloads;
    lbne::PennMicroSlice::Payload_Header    *word_header    = nullptr;
    lbne::PennMicroSlice::Payload_Timestamp *FirstTimestamp = nullptr;
    uint32_t payload_index = 0;
    uint16_t counter, trigger, timestamp, payloadCount;
    uint8_t* payload_data = nullptr;
    
    //const auto& frag((PTBrawFragments)[0]);
    const artdaq::Fragment frag = PTBrawFragments[0];
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
    PTBTime = FirstTimestamp->nova_timestamp;  
  }
  return;
}
//=======================================================================================
  
