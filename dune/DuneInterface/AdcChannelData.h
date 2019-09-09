// AdcChannelData.h
//
// David Adams
// June 2016
//
// Struct to hold the prepared time samples for a single TPC channel.
//
//           run - Run number
//        subRun - Sub-run number
//         event - Event number
//          time - Unix time (sec)
//       timerem - Time remainder (ns)
//       trigger - Index indicating which trigger(s) fired.
//  triggerClock - Time counter for the trigger
//
//       channel - Offline channel number
// channelStatus - Channel status (0=ok, 1=bad, 2=noisy)
//        fembID - FEMB ID
//   fembChannel - Channel number in FEMB (0, 1,..., 127).
//
//      pedestal - Pedestal subtracted from the raw count
//   pedestalRms - Pedestal RMS or sigma
//         tick0 - Tick offset between event/trigger time and raw/samples data.
//           raw - Uncompressed array holding the raw ADC count for each tick
//       samples - Array holding the prepared signal value for each tick
//    sampleUnit - Unit for samples array (ADC counts, fC, ke, ...)
//   sampleNoise - Noise level (e.g. RMS for samples)
//         flags - Array holding the status flag for each tick
//        signal - Array holding bools indicating which ticks have signals
//          rois - Array of ROIs indicating ticks of interest (e.g. have signals)
//       dftmags - Array of magnitudes for the DFT of the samples.
//     dftphases - Array of phases for the DFT of the samples.
//      metadata - Extra attributes
//
//         digit - Corresponding raw digit
//          wire - Corresponding wire
//    digitIndex - Index for the digit in the event digit container
//     wireIndex - Index for the wire in the event wire container
//
//         views - Map of vectors of AdcChannelData providing alternative representations of this data.
//                 E.g. non-contiguous time samples, processing snapshot or alternative, ....
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
#include <vector>
#include <map>
#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/Data/DuneEventInfo.h"

namespace raw {
  class RawDigit;
}
namespace recob {
  class Wire;
}

class AdcChannelData {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using FloatMap = std::map<Name, float>;
  using View = std::vector<AdcChannelData>;
  using ViewMap = std::map<Name, View>;

  static const AdcIndex badIndex =-1;
  static const AdcChannel badChannel =-1;
  static const size_t badSignal =-99999;

  // Event, trigger and channel data.
  // Everything here should be included in copy and assignment.
  AdcIndex run =badIndex;
  AdcIndex subRun =badIndex;
  AdcIndex event =badIndex;
  time_t time =0;
  int timerem =0;
  AdcIndex trigger =badIndex;
  AdcLongIndex triggerClock =0;
  AdcChannel channel =badIndex;
  AdcIndex channelStatus =badIndex;
  AdcIndex fembID =badIndex;
  AdcIndex fembChannel =badIndex;

  // ADC and derived data.
  AdcSignal pedestal =badSignal;
  AdcSignal pedestalRms =0.0;
  AdcInt tick0 = 0;
  AdcCountVector raw;
  AdcSignalVector samples;
  Name sampleUnit;
  AdcSignal sampleNoise =0.0;
  AdcFlagVector flags;
  AdcFilterVector signal;
  AdcRoiVector rois;
  AdcSignalVector dftmags;
  AdcSignalVector dftphases;
  FloatMap metadata;

  // Connections to persistent data.
  const raw::RawDigit* digit =nullptr;
  const recob::Wire* wire =nullptr;
  AdcIndex digitIndex =badIndex;
  AdcIndex wireIndex =badIndex;

private:

  // Constituents and altenate views.
  ViewMap m_views;

public:

  // Return event info.
  DuneEventInfo getEventInfo() const {
    DuneEventInfo devt;
    devt.event = event;
    devt.run = run;
    devt.subRun = subRun;
    devt.triggerClock = triggerClock;
    return devt;
  }

  // Copy ctor.
  // Only copy event and channel data and not samples and
  // derived data.
  // Root dictionary (6.08/06) and AdcChannelDataMap (STL) require we keep copy.
  AdcChannelData(const AdcChannelData& rhs) {
    run           = rhs.run;
    subRun        = rhs.subRun;
    event         = rhs.event;
    time          = rhs.time;
    timerem       = rhs.timerem;
    trigger       = rhs.trigger;
    triggerClock  = rhs.triggerClock;
    channel       = rhs.channel;
    channelStatus = rhs.channelStatus;
    fembID        = rhs.fembID;
  };

  // Delete assignment.
  AdcChannelData& operator=(const AdcChannelData&) =delete;

  // Keel default ctor and move operations.
  AdcChannelData() =default;
  AdcChannelData(AdcChannelData&&) =default;
  AdcChannelData& operator=(AdcChannelData&&) =default;

  // Check if a metadata field is defined.
  bool hasMetadata(Name mname) const {
    return metadata.find(mname) != metadata.end();
  }

  // Fetch metadata.
  float getMetadata(Name mname, float def =0.0) const {
    FloatMap::const_iterator imtd = metadata.find(mname);
    if ( imtd == metadata.end() ) return def;
    return imtd->second;
  }

  // Fetch any property including metadata.
  float getAttribute(Name mname, float def =0.0) const {
    // For basic types, return the value.
    if ( mname == "run" ) return run;
    if ( mname == "subRun" ) return subRun;
    if ( mname == "event" ) return event;
    if ( mname == "trigger" ) return trigger;
    if ( mname == "triggerClock" ) return triggerClock;  // lose precision here
    if ( mname == "channel" ) return channel;
    if ( mname == "fembID" ) return fembID;
    if ( mname == "fembChannel" ) return fembChannel;
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
    // For pointer, return false for null.
    if ( mname == "digit" ) return digit != nullptr;
    if ( mname == "wire" ) return wire != nullptr;
    if ( mname == "metadata" ) return metadata.size();
    // Otherwise return the metatdata field.
    return getMetadata(mname, def);
  }

  // Clear the data.
  void clear();

  // Fill rois from signal.
  void roisFromSignal();

  // Return normalization specifier for the DFT held here.
  // See DuneCommon/RealDftData.h.
  static AdcIndex dftNormalization() { return 22; }

  // Return the number of views.
  size_t viewSize() const { return m_views.size(); }

    NameVector vnams;
  // Return the vector of available view names.
  NameVector viewNames() const {
    NameVector vnams;
    for ( ViewMap::value_type ivie : m_views ) {
      vnams.push_back(ivie.first);
    }
    return vnams;
  }

  // Return if a view is defined.
  bool hasView(Name vnam) const {
    return m_views.find(vnam) != m_views.end();
  }

  // Read a view. Empty if undefined.
  const View& view(Name vnam) const {
    static const View empty;
    ViewMap::const_iterator ivie = m_views.find(vnam);
    if ( ivie == m_views.end() ) return empty;
    return ivie->second;
  }
 
  // Obtain a mutable view. Entries may be added, removed or modified.
  // View is added if not already existing.
  View& updateView(Name vnam) {
    return m_views[vnam];
  }

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
  time = 0;
  timerem = 0;
  trigger = badIndex;
  triggerClock = 0;
  channel = badIndex;
  channelStatus = badIndex;
  fembID = badIndex;
  fembChannel = badIndex;
  pedestal = badSignal;
  pedestalRms = 0.0;
  tick0 = 0;
  raw.clear();
  samples.clear();
  sampleUnit = "";
  sampleNoise = 0.0;
  flags.clear();
  signal.clear();
  rois.clear();
  dftmags.clear();
  dftphases.clear();
  digit = nullptr;
  wire = nullptr;
  dftmags.clear();
  dftphases.clear();
  digitIndex = badIndex;
  wireIndex = badIndex;
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
