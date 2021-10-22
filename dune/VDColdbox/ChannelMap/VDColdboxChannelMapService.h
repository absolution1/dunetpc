///////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       VDColdboxChannelMapService
// Module type: service
// File:        VDColdboxChannelMapService.h
// Author:      Tom Junk, October 2021
//
// Implementation of hardware-offline channel mapping reading from a file.  
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef VDColdboxChannelMapService_H
#define VDColdboxChannelMapService_H

#include <map>
#include <vector>
#include <iostream>
#include <limits>
#include <fstream>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace dune {
  class VDColdboxChannelMapService;
}

class dune::VDColdboxChannelMapService {

 public:

  VDColdboxChannelMapService(fhicl::ParameterSet const& pset);
  VDColdboxChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  struct VDCBChanInfo {
    int offlchan;        // in gdml and channel sorting convention
    int wib;             // slot number +1:    1, 2, 3, or 4
    int wibconnector;    // which board on the WIB  1, 2, 3, or 4
    int cebchan;         // cold electronics channel on board:  0 to 127
    int femb;            // FEMB ID:  1 to 14
    int asic;            // ASIC:   1 to 8
    int asicchan;        // ASIC channel:  0 to 15
    int connector;       // detector connector
    std::string stripid;  // strip name, with plane and number.  e.g. U79
    bool valid;          // true if valid, false if not
  };

  // Access methods

  // The function below gets cold electronics info from offline channel number.  Sets valid to be false if 
  // there is no corresponding cold electronics channel

  VDCBChanInfo getChanInfoFromOfflChan(int offlchan);

  // The function below uses conventions from Nitish's spreadsheet.  WIB: 1-3, wibconnector: 1-4, cechan: 0-127
  // It returns an offline channel number of -1 if the wib, wibconnector, and cechan are not in the map

  int getOfflChanFromWIBConnectorInfo(int wib, int wibconnector, int cechan);

  // For convenience, the function below  uses conventions from the DAQ WIB header, with two FEMBs per fiber
  // on FELIX readout:  slot: 0-2, fiber=1 or 2, cehcan: 0-255

  int getOfflChanFromSlotFiberChan(int slot, int fiber, int chan);


 private:

  //  map so we can look up channel info by offline channel number

  std::unordered_map<int,VDCBChanInfo> chantoinfomap;

  // map so we can look up offline channel number by channel info  (three keys), wib, wibconnector, and cechan

  std::unordered_map<int,std::unordered_map<int,std::unordered_map<int,int> > > infotochanmap;

};

DECLARE_ART_SERVICE(dune::VDColdboxChannelMapService, LEGACY)

#endif
