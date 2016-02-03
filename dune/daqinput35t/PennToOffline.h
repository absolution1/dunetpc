////////////////////////////////////////////////////////////////////////
// File:   PennToOffline.cc
// Author: Karl Warburton (Oct 2015)
//
// Utility to provide methods for reformatting the raw online PENN data
// into a structure which can be used in the splitter. 
// Heavily uses lbne-artdaq/OnlineMonitoring/DataReformatter.cxx as a base
////////////////////////////////////////////////////////////////////////

#ifndef PennToOffline_h
#define PennToOffline_h

#include "art/Framework/Principal/Handle.h"
#include "lbne-raw-data/Overlays/PennMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "utilities/UnpackFragment.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Geometry/Geometry.h"
#include "RawData/ExternalTrigger.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <bitset>
#include <utility>
#include <iterator>
#include "TMath.h"

namespace DAQToOffline {

    namespace TypeSizes {
    static int const CounterWordSize = 104;
    static int const TriggerWordSize = 32;
  }
    void CollectCounterBits(lbne::PennMicroSlice::Payload_Header *header,lbne::PennMicroSlice::Payload_Counter *trigger);
    void CollectTrigger(lbne::PennMicroSlice::Payload_Header *header,lbne::PennMicroSlice::Payload_Trigger *trigger);

    typedef std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::TriggerWordSize> > PTBTrigger;

    // Unpack the given artdaq::Fragment objects, and create a vector of raw::ExternalTrigger objects. The
    // Fragments are expected to be carrying Penn Trigger board data; this is not checked.
    std::vector<raw::ExternalTrigger> PennFragmentToExternalTrigger( artdaq::Fragments const& Fragments );

    // Function to get the full timestamp of a given payload
    void GetTimestamp( lbne::PennMilliSliceFragment msf,
		       lbne::PennMicroSlice::Payload_Header*& word_header,
		       lbne::PennMicroSlice::Payload_Timestamp* const& previous_timestamp,
		       lbne::PennMicroSlice::Payload_Timestamp*& future_timestamp,
		       lbne::PennMicroSlice::Payload_Header*& future_timestamp_header,
		       std::vector<lbne::PennMicroSlice::Payload_Timestamp::timestamp_t> &TimeVector );

    // Function to decide whether to make a new External Trigger
    bool MakeNewExtTrig( uint32_t pos, bool &PrevOn, bool NowOn );
}
#endif
