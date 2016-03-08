// tpcFragmentToRawDigits.cc

#include "tpcFragmentToRawDigits.h"
#include "dune/DuneCommon/DuneTimeConverter.h"

std::vector<raw::RawDigit>
DAQToOffline::tpcFragmentToRawDigits(artdaq::Fragments const& rawFragments,
				     std::vector<std::pair <std::pair<unsigned int,unsigned int>,lbne::TpcNanoSlice::Header::nova_timestamp_t> > &DigitsIndexList,
				     lbne::TpcNanoSlice::Header::nova_timestamp_t& firstTimestamp,
				     art::ServiceHandle<lbne::ChannelMapService> const& channelMap, bool useChannelMap,
				     bool debug, raw::Compress_t compression, unsigned int zeroThreshold) {

  DigitsIndexList.clear();
  //Create a map containing (fragmentID, fragIndex) for the event, will be used to check if each channel is present
  unsigned int numFragments = rawFragments.size();
  bool TimestampSet = false;

  std::map < unsigned int, unsigned int > mapFragID;

  for(size_t fragIndex = 0; fragIndex < rawFragments.size(); fragIndex++){

    const artdaq::Fragment &singleFragment = rawFragments[fragIndex];

    unsigned int fragmentID = singleFragment.fragmentID();

    mapFragID.insert(std::pair<unsigned int, unsigned int>(fragmentID,fragIndex));
  }


  if(debug){
    std::cout << numFragments<< " rawFragments" << std::endl;
  }

  std::vector<raw::RawDigit> rawDigitVector;

  //JPD -- first go at unpacking the information
  //    -- seems to make sense to look through channel number,
  //    -- then we'll create a rawDigit object for each channel
  //    -- will need some helper functions to do this for us, so I created a utilites directory

  art::ServiceHandle<geo::Geometry> geometry;
  size_t numChans = geometry->Nchannels();
  for(size_t chan=0;chan < numChans;chan++){

    //Each channel is uniquely identified by (fragmentID, group, sample) in an online event

    unsigned int fragmentID = UnpackFragment::getFragIDForChan(chan);
    unsigned int sample = UnpackFragment::getNanoSliceSampleForChan(chan);

    if (debug) {
      std::cout << "channel: " << chan
                << "\tfragment: " << fragmentID
        //<< "\tgroup: " << group
                << "\tsample: " << sample
                << std::endl;
    }

    //Check that the necessary fragmentID is present in the event
    //i.e. do we have data for this channel?

    if( mapFragID.find(fragmentID) == mapFragID.end() ){

      if (debug) std::cout << "Fragment not found" << std::endl;
      continue;

    }

    unsigned int fragIndex = mapFragID[fragmentID];

    if (debug) std::cout << "fragIndex: " << fragIndex << std::endl;

    std::vector<short> adcvec;


    const artdaq::Fragment &singleFragment = rawFragments[fragIndex];
    lbne::TpcMilliSliceFragment millisliceFragment(singleFragment);

    //Properties of fragment
    auto numMicroSlices = millisliceFragment.microSliceCount();
    bool FirstMicro = false;
    unsigned int FirstGoodIndex = 0;
    unsigned int LastGoodIndex  = 0;
    bool PrevMicroRCE  = false;
    lbne::TpcNanoSlice::Header::nova_timestamp_t ThisTimestamp = 0;
    bool MadeList = DigitsIndexList.size();
    for(unsigned int i_micro=0;i_micro<numMicroSlices;i_micro++){
      std::unique_ptr <const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(i_micro);
      auto numNanoSlices = microSlice->nanoSliceCount();
      
      // Push back adcvec with the ADC values.
      for(uint32_t i_nano=0; i_nano < numNanoSlices; i_nano++){
	uint16_t val = std::numeric_limits<uint16_t>::max();
	bool success = microSlice->nanosliceSampleValue(i_nano, sample, val);
        if(success) adcvec.push_back(short(val));
      }

      if ( chan%128==0 ) {
	//std::cout << "Looking at microslice " << i_micro << ", it has " << numNanoSlices << " nanoslices." << std::endl;
	// Get the First Timestamp for this channel
	if (!FirstMicro && numNanoSlices) {
	  lbne::TpcNanoSlice::Header::nova_timestamp_t Timestamp = microSlice->nanoSlice(0)->nova_timestamp();
	  FirstMicro=true;
	  if (!TimestampSet || Timestamp < firstTimestamp) {
	    //std::cout << "!!!Resetting timestamp to " << Timestamp << " on Chan " << chan << ",Micro " << i_micro << "!!!" << std::endl;
	    firstTimestamp = Timestamp;
	    TimestampSet = true;
	  }
	}
	// Which indexes have RCE information?
	if (!MadeList) {
	  if (numNanoSlices) {
	    LastGoodIndex = LastGoodIndex + numNanoSlices;
	    if (PrevMicroRCE == false)
	      ThisTimestamp = microSlice->nanoSlice(0)->nova_timestamp();
	  }
	  if (PrevMicroRCE == true && ( !numNanoSlices || i_micro==numMicroSlices-1 ) ) {
	    //std::cout << "This is the end of a good set of RCEs, so want to add a pair to my vector...." << std::endl;
	    DigitsIndexList.push_back( std::make_pair( std::make_pair(FirstGoodIndex,LastGoodIndex-1), ThisTimestamp) );
	    FirstGoodIndex = LastGoodIndex;
	  }
	  PrevMicroRCE = (bool)numNanoSlices;
	} // If not made index list
      } // Do stuff once for each RCE...
    }
    // Make my list of good RCEs
    if (DigitsIndexList.size() && !MadeList ) {
      //std::cout << "Finished looking through microslices. DigitsIndexList has size " << DigitsIndexList.size() << ", the useful things were in these microslices." << std::endl;
      //for (size_t qq=0; qq<DigitsIndexList.size(); ++qq)
      //std::cout << "Looking at element " << qq << " start " << DigitsIndexList[qq].first.first << ", end " << DigitsIndexList[qq].first.second << " timestamp " << DigitsIndexList[qq].second << std::endl;
    }
    
    if (debug) std::cout << "adcvec->size(): " << adcvec.size() << std::endl;
    unsigned int numTicks = adcvec.size();
    raw::Compress(adcvec, compression, zeroThreshold);
    int offlineChannel = -1;
    if (useChannelMap) offlineChannel = channelMap->Offline(chan);
    else offlineChannel = chan;
    raw::RawDigit theRawDigit(offlineChannel, numTicks, adcvec, compression);
    rawDigitVector.push_back(theRawDigit);            // add this digit to the collection
  }

  return rawDigitVector;
}

// NOvA time standard is a 56 bit timestamp at an LSB resolution of 15.6 ns (64 MHz)
// and a starting epoch of January 1, 2010 at 00:00:00.
// D. Adams March 2016: Update this to store the sec in the high word and remaing ns in
// the low word. Previous version stored 0 in the high word and sec in the low word.
art::Timestamp DAQToOffline::
make_art_timestamp_from_nova_timestamp(lbne::TpcNanoSlice::Header::nova_timestamp_t novaTime) {
  return DuneTimeConverter:: fromNova(novaTime);
}

art::Timestamp DAQToOffline::
old_make_art_timestamp_from_nova_timestamp(lbne::TpcNanoSlice::Header::nova_timestamp_t this_nova_timestamp) {
  lbne::TpcNanoSlice::Header::nova_timestamp_t seconds_since_nova_epoch = (this_nova_timestamp/nova_time_ticks_per_second);
  TTimeStamp time_of_event(20100101u,
                           0u,
                           0u,
                           true,
                           seconds_since_nova_epoch);
  return art::Timestamp(time_of_event.GetSec());
}
