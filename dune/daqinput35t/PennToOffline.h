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
#include "artdaq-core/Data/Fragment.hh"
#include "utilities/UnpackFragment.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <bitset>
#include <utility>
#include <iterator>
#include "TMath.h"
#include "TVector3.h"

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
    std::vector<raw::ExternalTrigger> PennFragmentToExternalTrigger( artdaq::Fragments const& Fragments, std::map<int,int>& channelMap, lbne::PennMicroSlice::Payload_Timestamp *&FirstPTBTimestamp );

    // Function to get the full timestamp of a given payload
    void GetTimestamp( lbne::PennMilliSliceFragment msf,
		       lbne::PennMicroSlice::Payload_Header*& word_header,
		       lbne::PennMicroSlice::Payload_Timestamp* const& previous_timestamp,
		       lbne::PennMicroSlice::Payload_Timestamp*& future_timestamp,
		       lbne::PennMicroSlice::Payload_Header*& future_timestamp_header,
		       std::vector<lbne::PennMicroSlice::Payload_Timestamp::timestamp_t> &TimeVector );

    // Function to decide whether to make a new External Trigger
    bool MakeNewExtTrig( uint32_t pos, bool &PrevOn, bool NowOn );

    void BuildPTBChannelMap(std::string MapDir, std::string MapFile, std::map<int,int>& channelMap);

    void MakeCounterPositionMap( std::string CounterDir, std::string CounterFile,
				 std::map< unsigned int, std::pair < TVector3, std::vector< TVector3 > > >& CounterPositionMap,
				 double fExtendCountersX=0, double fExtendCountersY=0, double fExtendCountersZ=0 );

    void MakeCounterCorners( int CountInd, double HalfLength, double HalfWidth1, double HalfWidth2, TVector3 Centre,
			     TVector3& TL, TVector3& TR, TVector3& BL, TVector3& BR,
			     double fExtendCountersX=0, double fExtendCountersY=0, double fExtendCountersZ=0 );
}
#endif
