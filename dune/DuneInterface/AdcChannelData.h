// AdcChannelData.h
//
// David Adams
// June 2016
//
// Struct to hold the prepared time samples for a single TPC channel.
//
//          run - Run number
//       subRun - Sub-run number
//        event - Event number
//      channel - Offline channel number
//       fembID - FEMB ID
//  fembChannel - Channel number in FEMB (0, 1,..., 127).
// triggerClock - Time counter for the trigger
//     pedestal - Pedestal subtracted from the raw count
//  pedestalRms - Pedestal RMS or sigma
//          raw - Uncompressed array holding the raw ADC count for each tick
//      samples - Array holding the prepared signal value for each tick
//  sampleNoise - Noise level (e.g. RMS for samples)
//        flags - Array holding the status flag for each tick
//       signal - Array holding bools indicating which ticks have signals
//         rois - Array of ROIs indicating ticks of interest (e.g. have signals)
//      dftmags - Array of magnitudes for the DFT of the samples.
//    dftphases - Array of phases for the DFT of the samples.
//        digit - Corresponding raw digit
//         wire - Corresponding wire
//   digitIndex - Index for the digit in the event digit container
//    wireIndex - Index for the wire in the event wire container
//   sampleUnit - Unit for samples array (ADC counts, fC, ke, ...)
//    dftphases - Phases for the DFT of samples
//      dftmags - Magnitudes for the DFT of samples
//     metadata - Extra attributes
//
// User can compare values against the defaults below to know if a value has been set.
// For arrays, check if the size in nonzero.
//
// If filled, the DFT fields should have lengths (nsam+2)/2 for the magnitudes and
// (nsam+1)/2 for the phases with the first phase zero or pi.
//

#ifndef AdcChannelData_H
#define AdcChannelData_H

#include <string>
#include <map>
#include "dune/DuneInterface/AdcTypes.h"

namespace raw {
  class RawDigit;
}
namespace recob {
  class Wire;
}

class AdcChannelData {

public:

  using FloatMap = std::map<std::string, float>;

  static const AdcIndex badIndex =-1;
  static const AdcChannel badChannel =-1;
  static const size_t badSignal =-99999;

  AdcIndex run =badIndex;
  AdcIndex subRun =badIndex;
  AdcIndex event =badIndex;
  AdcChannel channel =badIndex;
  AdcIndex fembID =badIndex;
  AdcIndex fembChannel =badIndex;
  AdcLongIndex triggerClock =0;
  AdcSignal pedestal =badSignal;
  AdcSignal pedestalRms =0.0;
  AdcCountVector raw;
  AdcSignalVector samples;
  AdcSignal sampleNoise =0.0;
  AdcFlagVector flags;
  AdcFilterVector signal;
  AdcRoiVector rois;
  std::string sampleUnit;
  const raw::RawDigit* digit =nullptr;
  const recob::Wire* wire =nullptr;
  AdcIndex digitIndex =badIndex;
  AdcIndex wireIndex =badIndex;
  AdcSignalVector dftmags;
  AdcSignalVector dftphases;
  FloatMap metadata;

  // Hide copy and assignment but allow move.
  // Root dictionary (6.08/06) requires we keep copy.
  AdcChannelData() =default;
  AdcChannelData(AdcChannelData&&) =default;
  AdcChannelData& operator=(AdcChannelData&&) =default;
  AdcChannelData& operator=(const AdcChannelData&) =delete;
#if 0
  AdcChannelData(const AdcChannelData&) =delete;
#else
  AdcChannelData(const AdcChannelData&) =default;
#endif

  // Check if a metadata field is defined.
  bool hasMetadata(std::string mname) const {
    return metadata.find(mname) != metadata.end();
  }

  // Fetch metadata.
  float getMetadata(std::string mname, float def =0.0) const {
    FloatMap::const_iterator imtd = metadata.find(mname);
    if ( imtd == metadata.end() ) return def;
    return imtd->second;
  }

  // Fetch any property including metadata.
  float getAttribute(std::string mname, float def =0.0) const {
    // For basic types, return the value.
    if ( mname == "run" ) return run;
    if ( mname == "subRun" ) return subRun;
    if ( mname == "event" ) return event;
    if ( mname == "fembID" ) return fembID;
    if ( mname == "fembChannel" ) return fembChannel;
    if ( mname == "triggerClock" ) return triggerClock;  // lose precision here
    if ( mname == "pedestal" ) return pedestal;
    if ( mname == "pedestalRms" ) return pedestalRms;
    if ( mname == "sampleNoise" ) return sampleNoise;
    if ( mname == "digitIndex" ) return digitIndex;
    if ( mname == "wireIndex" ) return wireIndex;
    // For vectors, return the size.
    if ( mname == "raw" ) return raw.size();
    if ( mname == "samples" ) return samples.size();
    if ( mname == "flags" ) return flags.size();
    if ( mname == "signal" ) return signal.size();
    if ( mname == "rois" ) return rois.size();
    if ( mname == "dftmags" ) return dftmags.size();
    if ( mname == "dftphases" ) return dftphases.size();
    if ( mname == "digit" ) return digit != nullptr;
    // For pointer, return false for null.
    if ( mname == "wire" ) return wire != nullptr;
    if ( mname == "metadata" ) return metadata.size();
    // Otherwise return the metatdata field.
    return getMetadata(mname, def);
  }

  // Clear the data.
  void clear();

  // Fill rois from signal.
  void roisFromSignal();

};

//**********************************************************************

typedef std::map<AdcChannel, AdcChannelData> AdcChannelDataMap;

#ifdef __ROOTCLING__
#pragma link C++ class AdcChannelDataMap;
#endif

//**********************************************************************

inline void AdcChannelData::clear() {
  run = badIndex;
  subRun = badIndex;
  event = badIndex;
  channel = badIndex;
  fembID = badIndex;
  fembChannel = badIndex;
  triggerClock = 0;
  pedestal = badSignal;
  pedestalRms = 0.0;
  raw.clear();
  samples.clear();
  sampleNoise = 0.0;
  flags.clear();
  signal.clear();
  rois.clear();
  dftmags.clear();
  dftphases.clear();
  digit = nullptr;
  wire = nullptr;
  digitIndex = badIndex;
  wireIndex = badIndex;
  dftmags.clear();
  dftphases.clear();
  metadata.clear();
}

//**********************************************************************

inline void AdcChannelData::roisFromSignal() {
  rois.clear();
  bool inRoi = false;
  for ( unsigned int isig=0; isig<signal.size(); ++isig ) {
    if ( inRoi ) {
      AdcRoi& roi = rois.back();
      if ( signal[isig] ) roi.second = isig;
      else inRoi = false;
    } else {
      if ( signal[isig] ) {
        rois.push_back(AdcRoi(isig, isig));
        inRoi = true;
      }
    }
  }
}

//**********************************************************************
#endif
