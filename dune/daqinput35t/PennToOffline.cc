////////////////////////////////////////////////////////////////////////
// File:   PennToOffline.cc
// Author: Karl Warburton (Oct 2015)
//
// Utility to provide methods for reformatting the raw online PENN data
// into a structure which can be used in the splitter. 
// Heavily uses lbne-artdaq/OnlineMonitoring/DataReformatter.cxx as a 
// template (written by N. Barros, D. Brailsford, M. Wallbank, )
////////////////////////////////////////////////////////////////////////

#include "PennToOffline.h"
#include <map>
#include <iostream>
#include <iomanip>
//=======================================================================================
std::vector<raw::ExternalTrigger> 
DAQToOffline::PennFragmentToExternalTrigger( artdaq::Fragments const& Fragments, std::map<int,int>& channelMap, lbne::PennMicroSlice::Payload_Timestamp *&FirstPTBTimestamp ) {
  
  std::vector<raw::ExternalTrigger> ExternTrigs;

  std::vector<lbne::PennMicroSlice::Payload_Trigger::trigger_type_t> trigger_types = {lbne::PennMicroSlice::Payload_Trigger::TA, // The 'telescope' coincidence...........110
										      lbne::PennMicroSlice::Payload_Trigger::TB, // North upper, South lower coincidence..112
										      lbne::PennMicroSlice::Payload_Trigger::TC, // North lower, South upper coincidence..113
										      lbne::PennMicroSlice::Payload_Trigger::TD};// East lower, West upper coincidence....111
  std::vector<lbne::PennMicroSlice::Payload_Trigger::trigger_type_t> calib_types = {lbne::PennMicroSlice::Payload_Trigger::C1,
										    lbne::PennMicroSlice::Payload_Trigger::C2,
										    lbne::PennMicroSlice::Payload_Trigger::C3,
										    lbne::PennMicroSlice::Payload_Trigger::C4};
  std::vector<lbne::PennMicroSlice::Payload_Timestamp::timestamp_t> fCounterTimes;
  std::vector<lbne::PennMicroSlice::Payload_Timestamp::timestamp_t> fMuonTriggerTimes;
  std::vector<lbne::PennMicroSlice::Payload_Timestamp::timestamp_t> fSSPTriggerTimes;
  //std::vector<lbne::PennMicroSlice::Payload_Timestamp::timestamp_t> fRCETriggerTimes;
  std::vector<lbne::PennMicroSlice::Payload_Timestamp::timestamp_t> fCalibrationTriggerTimes;
  
  std::vector<lbne::PennMicroSlice::Payload_Counter> fCounterWords;
  std::vector<lbne::PennMicroSlice::Payload_Trigger> fMuonTriggers;
  std::vector<lbne::PennMicroSlice::Payload_Trigger> fSSPTriggers;
  //std::vector<lbne::PennMicroSlice::Payload_Trigger> fRCETriggers;
  std::vector<lbne::PennMicroSlice::Payload_Trigger> fCalibrationTriggers;
  
  //std::cout << "There are " << Fragments.size() << " fragments" << std::endl;
  for ( size_t idx = 0; idx < Fragments.size(); ++idx ) {
    const auto& frag(Fragments[idx]);
    
    //Get the PennMilliSliceFragment from the artdaq fragment
    lbne::PennMilliSliceFragment msf(frag);

    lbne::PennMicroSlice::Payload_Header *word_header = nullptr;
    lbne::PennMicroSlice::Payload_Counter *word_p_counter = nullptr;
    lbne::PennMicroSlice::Payload_Trigger *word_p_trigger = nullptr;
    lbne::PennMicroSlice::Payload_Timestamp *previous_timestamp = nullptr;
    lbne::PennMicroSlice::Payload_Header *future_timestamp_header = nullptr;
    lbne::PennMicroSlice::Payload_Timestamp *future_timestamp = nullptr;
    uint8_t* payload_data = nullptr;
    uint32_t payload_index = 0;

    uint16_t counter, trigger, timestamp, payloadCount;
    payloadCount = msf.payloadCount(counter, trigger, timestamp);

    //std::cout << "There are a total of " << counter << " counters, " << trigger << " triggers, " << timestamp << " timestamps, giving a total of " << payloadCount << std::endl;
    
    while (payload_index < uint32_t(payloadCount-1) ) {
      payload_data = msf.get_next_payload(payload_index,word_header);
      //std::cout << "\nGot payload data " << payload_data << " for index " << payload_index << std::endl;
      if (payload_data == nullptr)
	continue;
      switch(word_header->data_packet_type) {

      case lbne::PennMicroSlice::DataTypeCounter:
	word_p_counter = reinterpret_cast<lbne::PennMicroSlice::Payload_Counter*>(payload_data);
	fCounterWords.push_back(*word_p_counter);
	GetTimestamp( msf, word_header, previous_timestamp, future_timestamp, future_timestamp_header, fCounterTimes );
	//std::cout << "Got a counter " << word_p_counter << " with timestamp " << fCounterTimes.back() << std::endl;
	break; // Counter Type

      case lbne::PennMicroSlice::DataTypeTrigger:
	word_p_trigger = reinterpret_cast<lbne::PennMicroSlice::Payload_Trigger*>(payload_data);

	if (word_p_trigger->has_muon_trigger()) { // If a 'Muon Trigger' ie counter coincidence
	  fMuonTriggers.push_back(*word_p_trigger);
	  GetTimestamp( msf, word_header, previous_timestamp, future_timestamp, future_timestamp_header, fMuonTriggerTimes );
	  //std::cout << "It is a Muon Trigger " << word_p_trigger << " with timestamp " << fMuonTriggerTimes.back() << std::endl;
	}
	if (word_p_trigger->has_ssp_trigger()) { // If an 'SSP Trigger'
	  fSSPTriggers.push_back(*word_p_trigger);
          GetTimestamp( msf, word_header, previous_timestamp, future_timestamp, future_timestamp_header, fSSPTriggerTimes );
	  //std::cout << "It is an SSP Trigger " << word_p_trigger << " with timestamp " << fSSPTriggerTimes.back() <<  std::endl;
	}
	/*
	if (word_p_trigger->has_rce_trigger() ) { // If an 'RCE Trigger'
	  fRCETriggers.push_back(*word_p_trigger);
	  GetTimestamp( msf, word_header, previous_timestamp, future_timestamp, future_timestamp_header, fRCETriggerTimes );
	  //std::cout << "It is an RCE Trigger " << word_p_trigger << " with timestamp " << fRCETriggerTimes.back() <<  std::endl;
	}
	*/
	if (word_p_trigger->has_calibration()) { // If a 'Calibration Trigger'
	  fCalibrationTriggers.push_back(*word_p_trigger);
	  GetTimestamp( msf, word_header, previous_timestamp, future_timestamp, future_timestamp_header, fCalibrationTriggerTimes );
	  //std::cout << "It is a Calibration Trigger " << word_p_trigger << " with timestamp " << fCalibrationTriggerTimes.back() << std::endl;
	}
	
	break; // Trigger Type
	
      case lbne::PennMicroSlice::DataTypeTimestamp:
	previous_timestamp = reinterpret_cast<lbne::PennMicroSlice::Payload_Timestamp*>(payload_data);
	if (FirstPTBTimestamp == nullptr) FirstPTBTimestamp = previous_timestamp;
	//std::cout << "Got a timestamp " << previous_timestamp << " " << previous_timestamp->nova_timestamp << ", the first timestamp in this event is " << FirstPTBTimestamp << " " << FirstPTBTimestamp->nova_timestamp << std::endl;
	break; // Timestamp Type
	
      default:
	//std::cout << "It isn't a counter, trigger or timestamp. It is a " << std::bitset<3>(word_header->data_packet_type) << std::endl;
	break; // Default case - do nothing.
      } // switch/case
    }
  } // NumFragments

  // *************** Now I want to loop over all counters and find out when they turned on to make ExternalTriggers ****************
  for (int counter_index=0; counter_index<98; ++counter_index ) {
    bool counter_previously_on = true;
    for (uint32_t pos = 0; pos < fCounterWords.size(); ++pos) {
      bool counter_currently_on = fCounterWords.at(pos).get_counter_status(counter_index);
      //std::cout << "Looking at counter " << counter_index << " which in the map is channel " << channelMap[counter_index]
      //	  << ", position " << pos << ", PrevOn? " << counter_previously_on << ", NowOn? " << counter_currently_on << std::endl;
      bool MakeNewExtCount = MakeNewExtTrig(pos, counter_previously_on, counter_currently_on);
      if (MakeNewExtCount) {
	raw::ExternalTrigger counter( channelMap[counter_index], fCounterTimes.at(pos) );
	ExternTrigs.push_back(counter);
	//std::cout << "Made my new trigger, TrigID " << ExternTrigs.back().GetTrigID() << " (" << counter.GetTrigID() << ") TrigTime "
	//	  << ExternTrigs.back().GetTrigTime() << " (" << counter.GetTrigTime() << ")" << std::endl;
      }
    } // Loop over counter word size   
  } // Loop over channels.

  // *************** Now to loop through 'MUON TRIGGERS' ****************
  bool trigA_previously_on = false, trigB_previously_on = false, trigC_previously_on = false, trigD_previously_on = false;
  for (uint32_t pos = 0; pos < fMuonTriggers.size(); ++pos) {
    bool trigA_currently_on = fMuonTriggers.at(pos).has_muon_TA();
    bool trigB_currently_on = fMuonTriggers.at(pos).has_muon_TB();
    bool trigC_currently_on = fMuonTriggers.at(pos).has_muon_TC();
    bool trigD_currently_on = fMuonTriggers.at(pos).has_muon_TD();
    lbne::PennMicroSlice::Payload_Timestamp::timestamp_t current_trigger_time = fMuonTriggerTimes.at(pos);
    bool MakeNewExtTrigA, MakeNewExtTrigB, MakeNewExtTrigC, MakeNewExtTrigD;
    if (pos == 0) {
      MakeNewExtTrigA = trigA_previously_on = trigA_currently_on;
      MakeNewExtTrigB = trigA_previously_on = trigB_currently_on;
      MakeNewExtTrigC = trigA_previously_on = trigC_currently_on;
      MakeNewExtTrigD = trigA_previously_on = trigD_currently_on;
    } else {
      MakeNewExtTrigA = MakeNewExtTrig( pos, trigA_previously_on, trigA_currently_on );
      MakeNewExtTrigB = MakeNewExtTrig( pos, trigB_previously_on, trigB_currently_on );
      MakeNewExtTrigC = MakeNewExtTrig( pos, trigC_previously_on, trigC_currently_on );
      MakeNewExtTrigD = MakeNewExtTrig( pos, trigD_previously_on, trigD_currently_on );
    }
    //std::cout << "Looking at element " << pos << " of " << fMuonTriggers.size() << "." << std::endl;
    //std::cout << "Trigs prev on (A,B,C,D) (" << trigA_previously_on <<","<< trigB_previously_on <<","<< trigC_previously_on <<","<< trigD_previously_on <<")." << std::endl;
    //std::cout << "And now (" << trigA_currently_on <<","<< trigB_currently_on <<","<< trigC_currently_on <<","<< trigD_currently_on << ")." << std::endl;;
    //std::cout << "So do I make a trigger? (" << MakeNewExtTrigA <<","<< MakeNewExtTrigB <<","<< MakeNewExtTrigC <<","<< MakeNewExtTrigD <<")."<< std::endl;
    if (MakeNewExtTrigA) {
      //std::cout << "Making an external trigger from Trigger A!!" << std::endl;
      raw::ExternalTrigger counter( 110, current_trigger_time );
      ExternTrigs.push_back(counter);
    }
    if (MakeNewExtTrigB) {
      //std::cout << "Making an external trigger from Trigger B!!" << std::endl;
      raw::ExternalTrigger counter( 112, current_trigger_time );
      ExternTrigs.push_back(counter);
    }
    if (MakeNewExtTrigC) {
      //std::cout << "Making an external trigger from Trigger C!!" << std::endl;
      raw::ExternalTrigger counter( 113, current_trigger_time );
      ExternTrigs.push_back(counter);
    }
    if (MakeNewExtTrigD) {
      //std::cout << "Making an external trigger from Trigger D!!" << std::endl;
      raw::ExternalTrigger counter( 111, current_trigger_time );
      ExternTrigs.push_back(counter);
    }
  } // Loop over Muon Triggers

  // *************** Now to loop through 'SSP TRIGGERS' ****************
  for (uint32_t pos = 0; pos < fSSPTriggers.size(); ++pos) {
    //std::cout << "Making an external trigger from SSP trigger!!" << std::endl;
    lbne::PennMicroSlice::Payload_Timestamp::timestamp_t current_trigger_time = fSSPTriggerTimes.at(pos);
    raw::ExternalTrigger counter( 115, current_trigger_time );
  }

  // *************** Now to loop through 'RCE TRIGGERS' ****************
  /*
  for (uint32_t pos = 0; pos < fRCETriggers.size(); ++pos) {
    //std::cout << "Making an external trigger from RCE trigger!!" << std::endl;
    lbne::PennMicroSlice::Payload_Timestamp::timestamp_t current_trigger_time = fRCETriggerTimes.at(pos);
    raw::ExternalTrigger counter( 116, current_trigger_time ); // WHAT CHANNEL DOES AN RCE TRIGGER GO ON!!!????
  }
  */

  //std::cout << "\n\nAll done, lets loop through ExternTrigs...." << std::endl;
  for ( std::vector<raw::ExternalTrigger>::const_iterator TrigIt = ExternTrigs.begin(); TrigIt != ExternTrigs.end(); TrigIt++ ) {
    //std::cout << "Looking at index " << std::distance((std::vector<raw::ExternalTrigger>::const_iterator)ExternTrigs.begin(), TrigIt) << " which has indexes " << TrigIt->GetTrigID() << ", " << TrigIt->GetTrigTime() << std::endl;
  }

  return ExternTrigs;
}
//=======================================================================================
void DAQToOffline::GetTimestamp( lbne::PennMilliSliceFragment msf,
				 lbne::PennMicroSlice::Payload_Header*& word_header,
				 lbne::PennMicroSlice::Payload_Timestamp* const& previous_timestamp,
				 lbne::PennMicroSlice::Payload_Timestamp*& future_timestamp,
				 lbne::PennMicroSlice::Payload_Header*& future_timestamp_header,
				 std::vector<lbne::PennMicroSlice::Payload_Timestamp::timestamp_t> &TimeVector ) {

  if (previous_timestamp != nullptr) {
    TimeVector.push_back(word_header->get_full_timestamp_post(previous_timestamp->nova_timestamp));
  } else {
    if (future_timestamp == nullptr) {
      future_timestamp = reinterpret_cast<lbne::PennMicroSlice::Payload_Timestamp*>(msf.get_next_timestamp(future_timestamp_header));
      if (future_timestamp == nullptr) {
	std::cerr << "CAN'T FIND PTB TIMESTAMP WORDS IN MILLISLICE FRAGMENT!!! Logic will fail." << std::endl;
	return;
      }
    }
    TimeVector.push_back(word_header->get_full_timestamp_pre(future_timestamp->nova_timestamp));
  }
  return;
}
//=======================================================================================
bool DAQToOffline::MakeNewExtTrig( uint32_t pos, bool &PrevOn, bool NowOn ) {
  if (pos == 0 && NowOn )
    return false;
  else if (PrevOn && NowOn )
    return false;
  else if (PrevOn && !NowOn) {
    PrevOn = false;
    return false;
  } else if (!PrevOn && !NowOn )
    return false;
  else if (!PrevOn && NowOn) {
    PrevOn = true;
    return true;
  } else {
    std::cout<<"ERROR IN PTBFORMATTER'S LOGIC IN COUNTER ANALYSIS"<<std::endl; //We should never get here
    return false;
  }
}
//=======================================================================================
void DAQToOffline::BuildPTBChannelMap(std::string MapDir, std::string MapFile, std::map<int,int>& channelMap) {
  /// Builds PTB channel map from the map txt file
  channelMap.clear();

  int onlineChannel;
  int offlineChannel;
  
  std::ostringstream LocStream;
  LocStream  << MapDir << MapFile;
  std::string fullname = LocStream.str();

  std::string FileLoc;
  cet::search_path sp("FW_SEARCH_PATH");
  if (sp.find_file(MapFile, FileLoc)) fullname = FileLoc;
    
  if (fullname.empty())
    mf::LogWarning("DAQToOffline") << "Input PTB channel map file " << MapFile << " not found in FW_SEARCH_PATH, or in location " << LocStream.str() << ".  Using online channel numbers!" << std::endl;

  else {
    mf::LogVerbatim("DAQToOffline") << "Build PTB Online->Offline channel Map from " << fullname;
    std::ifstream infile(fullname);
    if (!infile.good()) std::cout << "Not finding " << fullname << std::endl;
    while (infile.good()) {
      infile >> onlineChannel >> offlineChannel;
      channelMap.insert(std::make_pair(onlineChannel,offlineChannel));
      mf::LogVerbatim("DAQToOffline") << "   " << onlineChannel << " -> " << offlineChannel;
    }
    mf::LogVerbatim("DAQToOffline")<< "channelMap has size " << channelMap.size();
  }
}
//=======================================================================================
void DAQToOffline::MakeCounterPositionMap( std::string CounterDir, std::string CounterFile, std::map< unsigned int, std::pair < TVector3, std::vector< TVector3 > > >& CounterPositionMap, double fExtendCountersX, double fExtendCountersY, double fExtendCountersZ ) {

  if ( fExtendCountersY == 0 ) fExtendCountersY = fExtendCountersX;
  if ( fExtendCountersZ == 0 ) fExtendCountersZ = fExtendCountersX;
  // A function to make a map of counter corners and their centres.
  // Until the LArSoft counters are changed this is the preferred way to get their centres.
  // A text file is loaded in from pardata, though a user defined one can be uploaded instead.
  // The geometries of the LArSoft counters ( which are correct ) are used to determine the size of the corners
  //     of the counters from the loaded centres.
  // NOTE. The counters are assumed to rectangular as this make the maths much simpler.
  //       In truth the short width is 27.06 cm, whilst the long width is 32.5 cm.
  //       This means that there an overlap of counter corner produced though is relatively insignificant.
  // The Map is structured as follows:
  //     Key             - Counter Index, numbered 0 - 92
  //     Value           - a pair of TVector3 and vector of TVector3
  //     Value.first     - The centre position of the counter
  //     Value.second[0] - The Top Left corner of the counter
  //     Value.second[1] - The Top Right corner of the counter
  //     Value.second[2] - The Bottom Left corner of the counter
  //     Value.second[3] - The Bottom Right corner of the counter
  CounterPositionMap.clear();

  auto fGeom = lar::providerFrom<geo::Geometry>();
    
  std::ostringstream CountStream;
  CountStream << CounterDir << CounterFile;
  std::string CountPath = CountStream.str();

  std::string CountLoc;
  cet::search_path sp("FW_SEARCH_PATH");
  if (sp.find_file(CounterFile, CountLoc)) CountPath = CountLoc;

  if (CountPath.empty()) {
    mf::LogWarning("DAQToOffline") << "ERROR::Cannot find the counter position map file " << CounterFile << " in FW_SEARCH_PATH or " << CounterDir << std::endl;
    return;
  } else {
    std::cout << "Loading the counter position map from " << CountPath << std::endl;
    std::ifstream infile(CountPath);
    while (infile.good()) {
      int CountInd;
      double CentreX, CentreY, CentreZ;
      char Type, Side, Orientation;
      infile >> CountInd >> CentreX >> CentreY >> CentreZ >> Type >> Side >> Orientation;
      //std::cout << "Read in new line " << CountInd << " " << CentreX << " " << CentreY << " " << CentreZ << " " << Type << " " << Side << " " << Orientation << std::endl;
      
      geo::AuxDetGeo const& auxDet = fGeom->AuxDet(CountInd);

      // Make my TVector3's
      TVector3 Centre( CentreX, CentreY, CentreZ );
      TVector3 TL, TR, BL, BR;

      // Access the counter dimensions from the geometry
      double HalfLength = 0.5 * auxDet.Length();
      double HalfWidth1  = auxDet.HalfWidth1();
      double HalfWidth2  = auxDet.HalfWidth2();
      
      // Call the counter corner alg.
      MakeCounterCorners( CountInd, HalfLength, HalfWidth1, HalfWidth2, Centre, TL, TR, BL, BR, fExtendCountersX, fExtendCountersY, fExtendCountersZ );

      // Make my vector of TVector3's
      std::vector<TVector3> Corners;
      Corners.push_back( TL );
      Corners.push_back( TR );
      Corners.push_back( BL );
      Corners.push_back( BR );

      // Add this counter to my map.
      CounterPositionMap.insert( std::make_pair( CountInd, std::make_pair( Centre, Corners ) ) );
    } // Loading in CountPath.
  } // If found CountPath
} // MakeCounterPositionMap
//=======================================================================================
void DAQToOffline::MakeCounterCorners( int CountInd, double HalfLength, double HalfWidth1, double HalfWidth2, TVector3 Centre, TVector3& TL, TVector3& TR, TVector3& BL, TVector3& BR, double fExtendCountersX, double fExtendCountersY, double fExtendCountersZ ) {
  // The actual corners are calculated. As above:
  // NOTE: The corners are calculated as if the counters were rectangular I pass both widths incase somebody wishes to correct this...
  //       The telescope corners are not yet calculated...

  if ( fExtendCountersY == 0 ) fExtendCountersY = fExtendCountersX;
  if ( fExtendCountersZ == 0 ) fExtendCountersZ = fExtendCountersX;
  //std::cout << "Extending counters by " << fExtendCountersX << " " << fExtendCountersY << " " << fExtendCountersZ << std::endl;

  if ( CountInd < 44 ) { 
    if ( (CountInd >=6 && CountInd <=15) || (CountInd >=28 && CountInd <=37) ) { // For the East / West counters
      //std::cout << "Looking at an East West counter " << std::endl;
      // Top Left
      TL[0] = Centre[0] - fExtendCountersX - HalfWidth1;
      TL[1] = Centre[1] + fExtendCountersY + HalfLength;
      TL[2] = Centre[2];
      // Top Right
      TR[0] = Centre[0] + fExtendCountersX + HalfWidth1;
      TR[1] = Centre[1] + fExtendCountersY + HalfLength;
      TR[2] = Centre[2];
      // Bottom Left
      BL[0] = Centre[0] - fExtendCountersX - HalfWidth1;
      BL[1] = Centre[1] - fExtendCountersY - HalfLength;
      BL[2] = Centre[2];
      // Bottom Right
      BR[0] = Centre[0] + fExtendCountersX + HalfWidth1;
      BR[1] = Centre[1] - fExtendCountersY - HalfLength;
      BR[2] = Centre[2];
    } else { // For the North / South counters
      //std::cout << "Looking at a North South counter " << std::endl;
      // Top Left
      TL[0] = Centre[0];
      TL[1] = Centre[1] + fExtendCountersY + HalfLength;
      TL[2] = Centre[2] - fExtendCountersZ - HalfWidth1;
      // Top Right
      TR[0] = Centre[0];
      TR[1] = Centre[1] + fExtendCountersY + HalfLength;
      TR[2] = Centre[2] + fExtendCountersZ + HalfWidth1;
      // Bottom Left
      BL[0] = Centre[0];
      BL[1] = Centre[1] - fExtendCountersY - HalfLength;
      BL[2] = Centre[2] - fExtendCountersZ - HalfWidth1;
      // Bottom Right
      BR[0] = Centre[0];
      BR[1] = Centre[1] - fExtendCountersY - HalfLength;
      BR[2] = Centre[2] + fExtendCountersZ + HalfWidth1;
    } 
  } else { // For the telescope counters.
    if ( (CountInd >= 44 && CountInd <= 56) || (CountInd >= 67 && CountInd <= 82) ) { // For one group of Telescope counters
      //std::cout << "Looking at a Telescope 1 counter " << std::endl;
      for (int ww=0; ww<3; ++ww )
	TL[ww] = TR[ww] = BL[ww] = BR[ww] = Centre[ww];
    } else {                                                          // For the other group of Telescope counters
      //std::cout << "Looking at a Telescope 2 counter " << std::endl;
      for (int ww=0; ww<3; ++ww )
	TL[ww] = TR[ww] = BL[ww] = BR[ww] = Centre[ww];
    }
  }
}
