// WiredAdcChannelDataMap.h
//
// David Adams
// September 2016
//
// WiredAdcChannelDataMap is a simple struct that hold a map of
// AdcChannelDataMaps and a map of larsoft wire containers.
// Both are indexed by name.
//
// It is used to hold intermediate states in data preparation.

#ifndef WiredAdcChannelDataMap_H
#define WiredAdcChannelDataMap_H

#include <string>
#include <vector>
#include "dune/DuneInterface/AdcChannelData.h"
#include "lardataobj/RecoBase/Wire.h"

class WiredAdcChannelDataMap {

public:

  typedef std::string Name;
  typedef std::vector<Name> NameVector;
  typedef std::vector<recob::Wire> WireContainer;

  // Data and wire maps indexed by state name.
  std::map<Name, AdcChannelDataMap> dataMaps;
  std::map<Name, WireContainer*> wires;

  WiredAdcChannelDataMap() =default;

  // Ctor from state name and maximum number of channels.
  WiredAdcChannelDataMap(const NameVector names, AdcChannel nchmax) {
    for ( Name name : names ) {
      dataMaps[name];
      WireContainer* pwires = new WireContainer;
      wires[name] = pwires;
      pwires->reserve(nchmax);
    }
  }

  // Return if this object holds ADC channel data for the given stae name.
  bool hasData(Name sname) {
    return dataMaps.find(sname) != dataMaps.end();
  }

  // Return if this object holds wires for the given state name.
  bool hasWires(Name sname) {
    return wires.find(sname) != wires.end();
  }

};

#endif
