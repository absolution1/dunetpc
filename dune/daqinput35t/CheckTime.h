////////////////////////////////////////////////////////////////////////
// File:   PennToOffline.cc
// Author: Karl Warburton (Oct 2015)
//
// Utility to provide methods for reformatting the raw online PENN data
// into a structure which can be used in the splitter. 
// Heavily uses lbne-artdaq/OnlineMonitoring/DataReformatter.cxx as a base
////////////////////////////////////////////////////////////////////////

#ifndef CheckTime_h
#define CheckTime_h

//lbne-artdaq includes
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "lbne-raw-data/Overlays/SSPFragment.hh"
#include "lbne-raw-data/Overlays/anlTypes.hh"
#include "lbne-raw-data/Overlays/PennMilliSlice.hh"
#include "lbne-raw-data/Overlays/PennMicroSlice.hh"
#include "lbne-raw-data/Overlays/PennMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "lardataobj/RawData/RawDigit.h"
#include "utilities/UnpackFragment.h"

#include <memory>
#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"

namespace DAQToOffline {

  void GetRCEFirstTimestamp( artdaq::Fragments const& Fragments, int &ConsistRCE, int &NumADCs, long long &RCETime );
  void GetSSPFirstTimestamp( artdaq::Fragments const& Fragments, int &nSSPPayloads, long long &SSPTime );
  void GetPTBFirstTimestamp( artdaq::Fragments const& PTBrawFragments, int &nPTBPayloads, long long &PTBTime );
}
#endif
