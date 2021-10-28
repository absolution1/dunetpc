///////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       VDColdboxChannelMapService
// Module type: service
// File:        VDColdboxChannelMapService.h
// Author:      Tom Junk, October 2021
//
// Implementation of hardware-offline channel mapping reading from a file.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "VDColdboxChannelMapService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


dune::VDColdboxChannelMapService::VDColdboxChannelMapService(fhicl::ParameterSet const& pset) {

  std::string channelMapFile = pset.get<std::string>("FileName");
  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(channelMapFile, fullname);

  if (fullname.empty()) {
    std::cout << "Input file " << channelMapFile << " not found" << std::endl;
    throw cet::exception("File not found");
  }
  else
    std::cout << "VD Coldbox Channel Map: Building TPC wiremap from file " << channelMapFile << std::endl;

  std::ifstream inFile(fullname, std::ios::in);
  std::string line;
  while (std::getline(inFile,line)) {
    VDCBChanInfo chinfo;
    std::stringstream linestream(line);
    linestream >>
      chinfo.offlchan >>
      chinfo.wib >>
      chinfo.wibconnector >>
      chinfo.cebchan >>
      chinfo.femb >>
      chinfo.asic >>
      chinfo.asicchan >>
      chinfo.connector >>
      chinfo.stripid;
    chinfo.valid = true;
    chantoinfomap[chinfo.offlchan] = chinfo;
    infotochanmap[chinfo.wib][chinfo.wibconnector][chinfo.cebchan] = chinfo.offlchan; 
  }
  inFile.close();
}

dune::VDColdboxChannelMapService::VDColdboxChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&) : VDColdboxChannelMapService(pset) {
}

// The function below gets cold electronics info from offline channel number.  Sets valid to be false if 
// there is no corresponding cold electronics channel
// if not found in the map, return a chan info struct filled with -1's and set the valid flag to false.

dune::VDColdboxChannelMapService::VDCBChanInfo dune::VDColdboxChannelMapService::getChanInfoFromOfflChan(int offlchan)
{
  VDCBChanInfo r;
  auto fm = chantoinfomap.find(offlchan);
  if (fm == chantoinfomap.end())
    {
      r.offlchan = -1;
      r.wib = -1;
      r.wibconnector = -1;
      r.cebchan = -1;
      r.femb = -1;
      r.asic = -1;
      r.asicchan = -1;
      r.connector = -1;
      r.stripid = "INVALID";
      r.valid = false;
    }
  else
    {
      r = fm->second;
    }
  return r;
}


  // The function below uses conventions from Nitish's spreadsheet.  WIB: 1-3, wibconnector: 1-4, cechan: 0-127
  // It returns an offline channel number of -1 if the wib, wibconnector, and cechan are not in the map

int dune::VDColdboxChannelMapService::getOfflChanFromWIBConnectorInfo(int wib, int wibconnector, int cechan)
{
  int r = -1;
  auto fm1 = infotochanmap.find(wib);
  if (fm1 == infotochanmap.end()) return r;
  auto m1 = fm1->second;
  auto fm2 = m1.find(wibconnector);
  if (fm2 == m1.end()) return r;
  auto m2 = fm2->second;
  auto fm3 = m2.find(cechan);
  if (fm3 == m2.end()) return r;
  r = fm3->second;  
  return r;
}

  // For convenience, the function below  uses conventions from the DAQ WIB header, with two FEMBs per fiber
  // on FELIX readout:  slot: 0-2, fiber=1 or 2, cehcan: 0-255
  // it returns a channel number of -1 if the slot, fiber, chan combination is not in the map

int dune::VDColdboxChannelMapService::getOfflChanFromSlotFiberChan(int slot, int fiber, int chan)
{
  int wc = fiber*2 - 1;
  if (chan>127)
    {
      chan -= 128;
      wc++;
    }
  return getOfflChanFromWIBConnectorInfo(slot+1,wc,chan);
}


DEFINE_ART_SERVICE(dune::VDColdboxChannelMapService)
