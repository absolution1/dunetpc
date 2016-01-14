#include "tpcFragmentToRawDigits.h"

#include <map>
#include <iostream>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"

#include "TTimeStamp.h"

// From dunetpc
#include "utilities/UnpackFragment.h"

// From larcore
#include "Geometry/Geometry.h"

// From lardata
#include "RawData/raw.h"

std::vector<raw::RawDigit>
DAQToOffline::tpcFragmentToRawDigits(artdaq::Fragments const& rawFragments, 
				     lbne::TpcNanoSlice::Header::nova_timestamp_t& firstTimestamp,
				     std::map<int,int> const& channelMap,
				     bool debug, raw::Compress_t compression, unsigned int zeroThreshold)
{
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
    //std::cout << "Channel " << chan << " has " << numMicroSlices << " microslices." << std::endl;
    for(unsigned int i_micro=0;i_micro<numMicroSlices;i_micro++){

      std::unique_ptr <const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(i_micro);
      auto numNanoSlices = microSlice->nanoSliceCount();
      
      if (numNanoSlices) {
	lbne::TpcNanoSlice::Header::nova_timestamp_t Timestamp = microSlice->nanoSlice(0)->nova_timestamp();
	if (!FirstMicro && chan%128==0) {
	  std::cout << "Channel " << chan << ", microslice " << i_micro << ", nanoslice 0 has timestamp " << Timestamp
		    << ". nanoslice 1 has timestamp " << microSlice->nanoSlice(1)->nova_timestamp() << std::endl;
	  FirstMicro=true;
	}
	if (!TimestampSet || Timestamp < firstTimestamp) {
	  std::cout << "!!!Resetting timestamp from " << firstTimestamp << " to " << Timestamp << " on Chan " << chan << ",Micro " << i_micro << "!!!" << std::endl;
	  firstTimestamp = Timestamp;
	  TimestampSet = true;
	}
	/*
	std::unique_ptr<const lbne::TpcNanoSlice> nanoSlice0 = microSlice->nanoSlice(0);
	std::unique_ptr<const lbne::TpcNanoSlice> nanoSlice1 = microSlice->nanoSlice(1);
	std::unique_ptr<const lbne::TpcNanoSlice> nanoSlice2 = microSlice->nanoSlice(2);
	std::unique_ptr<const lbne::TpcNanoSlice> nanoSliceX = microSlice->nanoSlice(numNanoSlices-1);
	if ( chan == 129 || chan == 257 )
	  std::cout << "LOOKING AT " << chan << " " << i_micro << " " << nanoSlice0->nova_timestamp() << " " << nanoSlice1->nova_timestamp()<< " " << nanoSlice2->nova_timestamp() << " " << nanoSliceX->nova_timestamp() << std::endl;
	*/
      }

      for(uint32_t i_nano=0; i_nano < numNanoSlices; i_nano++){

	uint16_t val = std::numeric_limits<uint16_t>::max();
	bool success = microSlice->nanosliceSampleValue(i_nano, sample, val);
        if(success) adcvec.push_back(short(val));
      }
    }

    if (debug) std::cout << "adcvec->size(): " << adcvec.size() << std::endl;
    unsigned int numTicks = adcvec.size();
    raw::Compress(adcvec, compression, zeroThreshold);
    int offlineChannel = -1;
    if (channelMap.size() == 0) offlineChannel = chan;
    else offlineChannel = channelMap.at(chan);
    raw::RawDigit theRawDigit(offlineChannel, numTicks, adcvec, compression);
    rawDigitVector.push_back(theRawDigit);            // add this digit to the collection
  }

  return rawDigitVector;
}

void DAQToOffline::BuildTPCChannelMap(std::string channelMapFile, std::map<int,int>& channelMap) {

  /// Builds TPC channel map from the map txt file

  channelMap.clear();

  int onlineChannel;
  int offlineChannel;
    
  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(channelMapFile, fullname);
    
  if (fullname.empty())
    mf::LogWarning("DAQToOffline") << "Input TPC channel map file " << channelMapFile << " not found in FW_SEARCH_PATH.  Using online channel numbers!" << std::endl;

  else {
    mf::LogVerbatim("DAQToOffline") << "Build TPC Online->Offline channel Map from " << fullname;
    std::ifstream infile(fullname);
    while (infile.good()) {
      infile >> onlineChannel >> offlineChannel;
      channelMap.insert(std::make_pair(onlineChannel,offlineChannel));
      mf::LogVerbatim("DAQToOffline") << "   " << onlineChannel << " -> " << offlineChannel;
    }
    std::cout << "channelMap has size " << channelMap.size() << ". If this is 2048, then it's fine even if the above lines skipped a 'few' channels..." << std::endl;
  }
    
}

art::Timestamp DAQToOffline::make_art_timestamp_from_nova_timestamp(lbne::TpcNanoSlice::Header::nova_timestamp_t this_nova_timestamp){

/*

"NOvA time standard"
which is a 56 bit timestamp at an LSB resolution of 15.6 ns (64 MHz) and a starting epoch of
January 1, 2010 at 00:00:00.

*/
  lbne::TpcNanoSlice::Header::nova_timestamp_t seconds_since_nova_epoch = (this_nova_timestamp/nova_time_ticks_per_second);
  TTimeStamp time_of_event(20100101u,
                           0u,
                           0u,
                           true,
                           seconds_since_nova_epoch);

  return art::Timestamp(time_of_event.GetSec());
}
