////////////////////////////////////////////////////////////////////////
// File:   PennToOffline.cc
// Author: Karl Warburton (Oct 2015)
//
// Utility to provide methods for reformatting the raw online PENN data
// into a structure which can be used in the splitter. 
// Heavily uses lbne-artdaq/OnlineMonitoring/DataReformatter.cxx as a 
// template (written by D. Brailsford)
////////////////////////////////////////////////////////////////////////

#include "PennToOffline.h"
#include <map>
#include <iostream>
//=======================================================================================
std::vector<raw::ExternalTrigger> 
DAQToOffline::PennFragmentToExternalTrigger( artdaq::Fragments const& Fragments, int const& PTBIgnoreBit, std::map<int,int> const& PTBChannelMap,
					     std::pair<std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::CounterWordSize> >,
					     std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::TriggerWordSize> > > &PrevTimeStampWords) {
  
  std::vector<raw::ExternalTrigger> ExternTrigs;
  unsigned int numFragments = Fragments.size();
  
  std::vector<std::bitset<TypeSizes::CounterWordSize> > fCounterBits;
  std::vector<int> fCounterTimes;
  std::vector<std::bitset<TypeSizes::TriggerWordSize> > fMuonTriggerBits;
  std::vector<int> fMuonTriggerTimes;

  lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t LastCounterTimestamp = PrevTimeStampWords.first.first;
  std::bitset<TypeSizes::CounterWordSize>                      LastCounterWord      = PrevTimeStampWords.first.second;
  lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t LastTriggerTimestamp = PrevTimeStampWords.second.first;
  std::bitset<TypeSizes::TriggerWordSize>                      LastTriggerWord      = PrevTimeStampWords.second.second;;
    
  for ( size_t idx = 0; idx < numFragments; ++idx ) {
    const auto& frag(Fragments[idx]);
    
    //Get the PennMilliSliceFragment from the artdaq fragment
    lbne::PennMilliSliceFragment msf(frag);
    
    //Get the number of each payload type in the millislice
    lbne::PennMilliSlice::Header::payload_count_t n_frames, n_frames_counter, n_frames_trigger, n_frames_timestamp;
    n_frames = msf.payloadCount(n_frames_counter, n_frames_trigger, n_frames_timestamp);
    //std::cout << n_frames << " " << n_frames_counter << " " << n_frames_trigger << " " << n_frames_timestamp << std::endl;
    
    //Now we need to grab the payload information in the millislice
    lbne::PennMicroSlice::Payload_Header::data_packet_type_t type;
    lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t timestamp;
    uint8_t* payload_data;
    size_t payload_size;
    
    //Loop over the payloads
    for (uint32_t ip = 0; ip < n_frames; ip++) {
      //Get the actual payload
      payload_data = msf.payload(ip, type, timestamp, payload_size);
      unsigned int payload_type = (unsigned int) type;
      //std::cout << "Payload " << ip << ", payload size " << payload_size << ", timestamp " << (int)timestamp << ", payload_type " << payload_type << std::endl;

      if ( ip == 0 ) {
	if ( (int)timestamp - (int)LastTriggerTimestamp < 0) {
	  std::cout << "Got a smaller timestamp " << timestamp << " than last trigger so resetting." << std::endl;
	  //std::cout << PrevTimeStampWords.first.first << " " << PrevTimeStampWords.first.second << " " << PrevTimeStampWords.second.first << " " << PrevTimeStampWords.second.second << std::endl;
	  LastTriggerTimestamp = 0; PrevTimeStampWords.second.first = 0; PrevTimeStampWords.second.second = 0;
	}
	if ( (int)timestamp - (int)LastCounterTimestamp < 0) {
	  std::cout << "Got a smaller timestamp " << timestamp << " than last counter so resetting." << std::endl;
	  //std::cout << PrevTimeStampWords.first.first << " " << PrevTimeStampWords.first.second << " " << PrevTimeStampWords.second.first << " " << PrevTimeStampWords.second.second << std::endl;
	  LastCounterTimestamp = 0; PrevTimeStampWords.first.first = 0; PrevTimeStampWords.first.second = 0;
	}
      }
      //Loop over the words in the payload
      if (!((payload_data != nullptr) && payload_size)) continue;
      switch (payload_type){
      case 1:
	//Dealing with counter word
	//std::cout << "Collecting Counter word!" << std::endl;
	CollectCounterBits(payload_data, payload_size, fCounterBits, fCounterTimes, timestamp, PrevTimeStampWords, PTBIgnoreBit);
	break;
      case 2:
	//Dealing with trigger word
	//std::cout << "Collecting Trigger word!" << std::endl;
	CollectMuonTrigger(payload_data, payload_size, fMuonTriggerBits, fMuonTriggerTimes, timestamp, PrevTimeStampWords, PTBIgnoreBit);
      default:
	break;
      }
    }
  } // Loop over fragments...
  
  //std::cout << "Looked at all fragments. fCounterTimes.size() = " <<  fCounterTimes.size() << ", fCounterBits.size() " << fCounterBits.size() 
  //<< ", fMuonTriggerTimes.size() = " << fMuonTriggerTimes.size() << ", fMuonTriggerBits.size() " << fMuonTriggerBits.size() << std::endl;
  
  // *************** Now I want to loop over all coutners / triggers and find when they are 'turned on' ****************
  
  //We need to loop through the requested counter / trigger to check and obey the following laws - Michelle Oct 2015
  //  A) An ignore bit is 400 ns ( 7 ticks ) - this is the default fcl param size passed to this function.
  //  B) Once a counter is turned on, ignore any other times it turns on before ignore bit has passed.
  //       This is something to do with how the signal decays, meaning that you can get ghost hits after initial hit.
  //  C) Only make a trigger when the counter turns on. 
  //       If ignoreBit passes and still on, ignore all bits until it turns off and then on again.
  //  D) If don't know whether a counter was on from the previous event, turn on the ignore bit.
 
  for (int counter_index=0; counter_index<97; ++counter_index ) {
    //Use the Online to Offline map.
    int offlineCounter = -1;
    if (PTBChannelMap.size() == 0) offlineCounter = counter_index;
    else offlineCounter = PTBChannelMap.at(counter_index);
    
    //To satisfy law D we have to assume that counter was previously on, so that we can't trigger on the first word.
    bool counter_previously_on = true;
    int LastCounter = -1e7;
    if ( LastCounterTimestamp != 0 ) {
      counter_previously_on = LastCounterWord[offlineCounter];
      LastCounter = (int)LastCounterTimestamp;
    }

    for (std::vector<std::bitset<TypeSizes::CounterWordSize> >::const_iterator countIt = fCounterBits.begin(); countIt != fCounterBits.end(); countIt++){
      //Get the counter status
      bool counter_currently_on = (*countIt)[offlineCounter];
      //We need the array index to fetch the timestamp of the payload
      int index = std::distance((std::vector<std::bitset<TypeSizes::CounterWordSize> >::const_iterator)fCounterBits.begin(), countIt);
      int this_time = fCounterTimes.at(index);
      
      //std::cout << "Looking at counter " << counter_index << ", corresponding to bit " << offlineCounter << ", bitset " << index << " correpsonding to time "
      //	<< fCounterTimes.at(index) << ". Counter was on? " << counter_previously_on << ", and now? " << counter_currently_on << std::endl;
      
      // If the first counter word was on, we want to ignore the first 400 ticks, as per law D.
      if ( std::distance((std::vector<std::bitset<TypeSizes::CounterWordSize> >::const_iterator)fCounterBits.begin(), countIt) == 0 && counter_currently_on && LastCounter==-1e7) {
	LastCounter = fCounterTimes.at(index);
      }
      //If the counter was previously on and it is still on, just continue
      if (counter_previously_on && counter_currently_on) 
	{ continue; }
      //The counter WAS on but it is now off so record the counter as switching off and move on
      else if (counter_previously_on && !counter_currently_on){
	counter_previously_on = false;
	continue;
      }
      //The counter was off and it is still off so very little to do here.  Continue!
      else if (!counter_previously_on && !counter_currently_on)
	{ continue; }
      //The counter has switched on!!!!!
      else if (!counter_previously_on && counter_currently_on){
	counter_previously_on = true; //Record that it is now on and record some numbers
	//If this bit is outside an ignoreBit then make an ExternalTrigger
	//std::cout << "Counter just turned on at timestamp " << this_time << ", last turned on at " << LastCounter << std::endl;
	if ( this_time-LastCounter > PTBIgnoreBit ){
	  //std::cout << "Making an external trigger!!" << std::endl;
	  raw::ExternalTrigger counter( counter_index, this_time );
	  ExternTrigs.push_back(counter);
	  LastCounter = fCounterTimes.at(index);
	  //std::cout << "Pushing back a counter" << std::endl;
	} 
      }
      else std::cout<<"ERROR IN PTBFORMATTER'S LOGIC IN COUNTER ANALYSIS"<<std::endl; //We should never get here
    } // Loop over counter word size   
  } // Loop over channels.

  //============================ Now loop through triggers ============================
  for ( int trigger_index=0; trigger_index<4; ++trigger_index ) { 

    bool counter_previously_on = true;
    int LastTrigger = -1e7;
    if ( LastTriggerTimestamp != 0 ) {
      counter_previously_on = LastTriggerWord[trigger_index];
      LastTrigger = (int)LastTriggerTimestamp;
    }
        
    int FillIndex = trigger_index;
    if ( trigger_index == 0      ) FillIndex = 111; // East lower, West upper coincidence.....Muon Trigger Pattern 1
    else if ( trigger_index == 1 ) FillIndex = 113; // North lower, South upper coincidence...Muon Trigger Pattern 2
    else if ( trigger_index == 2 ) FillIndex = 112; // North upper, South lower coincidence...Muon Trigger Pattern 4
    else if ( trigger_index == 3 ) FillIndex = 110; // The 'telescope' coincidence............Muon Trigger Pattern 8

    for (std::vector<std::bitset<TypeSizes::TriggerWordSize> >::const_iterator countIt = fMuonTriggerBits.begin(); countIt != fMuonTriggerBits.end(); countIt++){
      bool counter_currently_on = (*countIt)[trigger_index];
      //We need the array index to fetch the timestamp of the payload
      int index = std::distance((std::vector<std::bitset<TypeSizes::TriggerWordSize> >::const_iterator)fMuonTriggerBits.begin(), countIt);
      int this_time = fMuonTriggerTimes.at(index);
      
      //std::cout << "Looking at trigger " << trigger_index << ", bitset " << index << " correpsonding to time " << fMuonTriggerTimes.at(index)
      //	<< ". Trigger was on? " << counter_previously_on << ", and now? " << counter_currently_on << std::endl;

      // If the first counter word was on, we want to ignore the first 400 ticks, as per law D.
      if ( std::distance((std::vector<std::bitset<TypeSizes::TriggerWordSize> >::const_iterator)fMuonTriggerBits.begin(), countIt) == 0 && counter_currently_on && LastTrigger==-1e7) {
	LastTrigger = fMuonTriggerTimes.at(index);
      }

      //If the counter was previously on and it is still on, just continue
      if (counter_previously_on && counter_currently_on) 
	{ continue; }
      //The counter WAS on but it is now off so record the counter as switching off and move on
      else if (counter_previously_on && !counter_currently_on){
	counter_previously_on = false;
	continue;
      }
      //The counter was off and it is still off so very little to do here.  Continue!
      else if (!counter_previously_on && !counter_currently_on)
	{ continue; }
      //The counter has switched on!!!!!
      else if (!counter_previously_on && counter_currently_on){
	counter_previously_on = true; //Record that it is now on and record some numbers
	//If this bit is outside an ignoreBit then make an ExternalTrigger
	//std::cout << "Trigger " << this_time << " - " << LastTrigger << " " << PTBIgnoreBit << std::endl;
	if ( this_time-LastTrigger > PTBIgnoreBit ){
	  //std::cout << "Making a trigger :D " << std::endl;
	  raw::ExternalTrigger counter( FillIndex, this_time );
	  ExternTrigs.push_back(counter);
	  LastTrigger = fMuonTriggerTimes.at(index);
	} 
      }
      else std::cout<<"ERROR IN PTBFORMATTER'S LOGIC IN COUNTER ANALYSIS"<<std::endl; //We should never get here
    }
  }
  
  //std::cout << "\n\nAll done, lets loop through ExternTrigs...." << std::endl;
  //for ( std::vector<raw::ExternalTrigger>::const_iterator TrigIt = ExternTrigs.begin(); TrigIt != ExternTrigs.end(); TrigIt++ )
  //  std::cout << "Looking at index " << std::distance((std::vector<raw::ExternalTrigger>::const_iterator)ExternTrigs.begin(), TrigIt) << " which has indexes " << TrigIt->GetTrigID() << ", " << TrigIt->GetTrigTime() << std::endl;
  
  return ExternTrigs;
}
//=======================================================================================
void DAQToOffline::CollectCounterBits(uint8_t* payload, size_t payload_size, std::vector<std::bitset<TypeSizes::CounterWordSize> > &fCounterBits,
				      std::vector<int> &fCounterTimes, lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t timestamp,
				      std::pair<std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::CounterWordSize> >,
				      std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::TriggerWordSize> > > &PrevTimeStampWords,
				      int const& PTBIgnoreBit
				      ) {
  std::bitset<TypeSizes::CounterWordSize> bits;
  for (size_t ib = 0; ib < payload_size; ib++){
    std::bitset<TypeSizes::CounterWordSize> byte = payload[ib];
    bits ^= (byte << 8*(payload_size-(ib+1)));
  }

  if ( (int)timestamp - (int)PrevTimeStampWords.first.first < PTBIgnoreBit
       && PrevTimeStampWords.first.second == bits ) {
    std::cout << "Throwing away counter as too close and identical to previous counter." << std::endl;
    std::cout << (int)timestamp << " " << (int)PrevTimeStampWords.first.first << " " << PrevTimeStampWords.first.second << " " << bits << std::endl;
    return;
  }

  //std::cout << "HasCounter"<< bits << std::endl;
  fCounterBits.push_back(bits);
  fCounterTimes.push_back(timestamp);

  PrevTimeStampWords.first.first  = timestamp;
  PrevTimeStampWords.first.second = bits;

  return;
}
//=======================================================================================
void DAQToOffline::CollectMuonTrigger(uint8_t* payload, size_t payload_size, std::vector<std::bitset<TypeSizes::TriggerWordSize> > &fMuonTriggerBits,
				      std::vector<int> &fMuonTriggerTimes, lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t timestamp,
				      std::pair<std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::CounterWordSize> >,
				      std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::TriggerWordSize> > > &PrevTimeStampWords,
				      int const& PTBIgnoreBit
				      ) {
  std::bitset<TypeSizes::TriggerWordSize> bits;
  for (size_t ib = 0; ib < payload_size; ib++){
    std::bitset<TypeSizes::TriggerWordSize> byte = payload[ib];
    bits ^= (byte << 8*(payload_size-(ib+1)));
  }

  if ( (int)timestamp - (int)PrevTimeStampWords.second.first < PTBIgnoreBit
       && PrevTimeStampWords.second.second == bits ) {
    std::cout << "Throwing away trigger as too close and identical to previous trigger." << std::endl;
    std::cout << (int)timestamp << " " << (int)PrevTimeStampWords.second.first << " " << PrevTimeStampWords.second.second << " " << bits << std::endl;
    return;
  }

  //std::cout<<"HasTrigger" << bits <<std::endl;
  //Bits collected, now get the trigger type
  std::bitset<TypeSizes::TriggerWordSize> trigger_type_bits; 
  trigger_type_bits ^= (bits >> (TypeSizes::TriggerWordSize - 5));
  int trigger_type = static_cast<int>(trigger_type_bits.to_ulong());
  //If we have a muon trigger, get which trigger it was and increase the hit rate
  if ( trigger_type == 16 ){
    std::bitset<TypeSizes::TriggerWordSize> muon_trigger_bits;
    //Shift the bits left by 5 first to remove the trigger type bits
    muon_trigger_bits ^= (bits << 5);
    //Now we need to shift the entire thing set to the LSB so we can record the number
    //The muon trigger pattern words are 4 bits
    muon_trigger_bits >>= (TypeSizes::TriggerWordSize-4);
    //int muon_trigger = static_cast<int>(muon_trigger_bits.to_ulong());
    //std::cout << "HasTriggerType" << muon_trigger_bits << " " << muon_trigger << std::endl;
    fMuonTriggerBits.push_back(muon_trigger_bits); // Is this right!!???? Probably not wholly, but it works.
    fMuonTriggerTimes.push_back(timestamp);

    PrevTimeStampWords.second.first  = timestamp;
    PrevTimeStampWords.second.second = bits;
  }
  return;
}
//=======================================================================================
void DAQToOffline::BuildPTBChannelMap(std::string PTBMapFile, std::map<int,int>& PTBChannelMap) {
  /// Builds PTB channel map from the map txt file
  PTBChannelMap.clear();

  int onlineChannel;
  int offlineChannel;

  std::string MyDir = "/dune/data/users/warburton/Splitter/MyMap/";
  std::ostringstream LocStream;
  LocStream  << MyDir << PTBMapFile;
  std::string fullname = LocStream.str();
  
  std::string FileLoc;
  cet::search_path sp("FW_SEARCH_PATH");
  if (sp.find_file(PTBMapFile, FileLoc)) fullname = FileLoc;

  //std::string fullname;
  //cet::search_path sp("FW_SEARCH_PATH");
  //sp.find_file(PTBMapFile, fullname);
  
  if (fullname.empty())
    mf::LogWarning("DAQToOffline") << "Input PTB channel map file " << PTBMapFile << " not found.  Using online channel numbers!" << std::endl;

  else {
    mf::LogVerbatim("DAQToOffline") << "Build PTB Online->Offline channel Map from " << fullname;
    std::ifstream infile(fullname);
    if (!infile.good()) std::cout << "Couldn't find the channel map so using online channels instead." << std::endl;
    while (infile.good()) {
      infile >> onlineChannel >> offlineChannel;
      PTBChannelMap.insert(std::make_pair(onlineChannel,offlineChannel));
      mf::LogVerbatim("DAQToOffline") << "Looking at PTB    " << onlineChannel << " -> " << offlineChannel;
    }
  }
}
//=======================================================================================
