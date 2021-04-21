// AdcChannelData.h
//
// David Adams
// June 2016
//
// Struct to hold the prepared time samples for a single TPC channel.
//
// Held in event info pointer:
//           run - Run number
//        subRun - Sub-run number
//         event - Event number
//          time - Unix time (sec)
//       timerem - Time remainder (ns)
//       trigger - Index indicating which trigger(s) fired.
//  triggerClock - Time counter for the trigger
//  triggerTick0 - Tick counter for the event
//
// Held in channel info pointer:
//       channel - Offline channel number
//        fembID - FEMB ID
//   fembChannel - Channel number in FEMB (0, 1,..., 127).
// channelStatus - Channel status (0=ok, 1=bad, 2=noisy)
//
// Held directly:
//  channelClock - Time counter for the channel data
//         tick0 - Tick offset between triggerTick0 and first tick in raw/samples data.
//      pedestal - Pedestal subtracted from the raw count
//   pedestalRms - Pedestal RMS or sigma
//
//           raw - Uncompressed array holding the raw ADC count for each tick
//       samples - Array holding the prepared signal value for each tick
//    binSamples - Array holding binned samples: binSamples[ibin][itick]
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
//         views - Map of data views indexed by name.
//                 Each object implicitly has a self view with name "" holding only itself.
//    viewParent - Pointer to the data object holding this as a view. Creator of the view
//                 object is responsible for filling this field.
//
// User can compare values against the defaults below to know if a value has been set.
// For arrays, check if the size in nonzero.
//
// There separate clocks for the trigger which corresponds to (time, timerem) and
// the TPC readout (channelClock). For protoDUNE, these run at 50 MHz.
//
// If filled, the DFT fields should have lengths (nsam+2)/2 for the magnitudes and
// (nsam+1)/2 for the phases with the first phase zero or pi.
//
// A data view is a vector of AdcChannelData objects providing an alternative
// representation of its parent.
// E.g. non-contiguous time samples, processing snapshot or alternative, ....
//
// A view path is a vector of names {vn0, vn1, ...} with string representation
// "vn0/vn1/..." that specifies the path to a collection of AdcChannelData objects
// at any depth in the view hiearchy.
// Example view access:
//   string myview = "constituents/trun1000";
//   AdcIndex nvie = mydata.viewCount(myview);
//   for ( AdcIndex ivie=0; ivie<nvie; ++ivie ) {
//     AdcChannelData* pdat = mydata.viewEntry(myview, ivie);
//     cout << "T0 for view " << ivie " << ": " << pdata->tick0 << endl;
//   }

#ifndef AdcChannelData_H
#define AdcChannelData_H

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <iostream>
#include "dune/DuneInterface/Data/AdcTypes.h"
#include "dune/DuneInterface/Data/DuneEventInfo.h"
#include "dune/DuneInterface/Data/DuneChannelInfo.h"

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
  using EventInfo = DuneEventInfo;
  using EventInfoPtr = std::shared_ptr<const EventInfo>;
  using ChannelInfo = DuneChannelInfo;
  using ChannelInfoPtr = std::shared_ptr<const ChannelInfo>;
  using Index = EventInfo::Index;
  using LongIndex = EventInfo::LongIndex;
  using Channel = ChannelInfo::Channel;

  static Index badIndex() { return EventInfo::badIndex(); }
  static Index badLongIndex() { return EventInfo::badLongIndex(); }
  static Index badChannel() { return EventInfo::badIndex(); }
  static Index badSignal() { return -99999; }

  static const EventInfo& badEventInfo() {
    return EventInfo::badEventInfo();
  }

  static const ChannelInfo& badChannelInfo() {
    return ChannelInfo::badChannelInfo();
  }

  // ADC and derived data.
  AdcLongIndex channelClock =0;
  AdcInt tick0 = 0;
  AdcSignal pedestal =badSignal();
  AdcSignal pedestalRms =0.0;
  AdcCountVector raw;
  AdcSignalVector samples;
  AdcSignalVectorVector binSamples;
  Name sampleUnit;
  AdcSignal sampleNoise =0.0;
  AdcFlagVector flags;
  AdcFilterVector signal;
  AdcRoiVector rois;
  AdcSignalVector dftmags;
  AdcSignalVector dftphases;
  FloatMap metadata;
  AdcChannelData* viewParent =nullptr;

  // Connections to persistent data.
  const raw::RawDigit* digit =nullptr;
  const recob::Wire* wire =nullptr;
  AdcIndex digitIndex =badIndex();
  AdcIndex wireIndex =badIndex();

private:

  EventInfoPtr m_peventInfo;
  ChannelInfoPtr m_pchanInfo;

  // Constituents and alternate views.
  ViewMap m_views;

public:

  // Set event info.
  void setEventInfo(EventInfoPtr pevi) { m_peventInfo = pevi; }
  void setEventInfo(const EventInfo* pevi) { m_peventInfo.reset(pevi); }
  void setEventInfo(Index a_run, Index a_event, Index a_subRun =badIndex(),
                    time_t a_time =0, int a_timerem =0,
                    Index a_trigger =badIndex(), LongIndex a_triggerClock =badLongIndex()) {
    setEventInfo(new EventInfo(a_run, a_event, a_subRun, a_time, a_timerem, a_trigger, a_triggerClock));
  }

  // Return event info.
  bool hasEventInfo() const { return bool(m_peventInfo); }
  EventInfoPtr getEventInfoPtr() const { return m_peventInfo; }
  const DuneEventInfo& getEventInfo() const {
    return hasEventInfo() ? *m_peventInfo : badEventInfo();
  }
  AdcIndex run()              const { return getEventInfo().run; }
  AdcIndex subRun()           const { return getEventInfo().subRun; }
  AdcIndex event()            const { return getEventInfo().event; }
  time_t time()               const { return getEventInfo().time; }
  int timerem()               const { return getEventInfo().timerem; }
  AdcIndex trigger()          const { return getEventInfo().trigger; }
  AdcLongIndex triggerClock() const { return getEventInfo().triggerClock; }
  AdcLongIndex triggerTick0() const { return getEventInfo().triggerTick0; }

  // Set channel info.
  void setChannelInfo(ChannelInfoPtr pchi) { m_pchanInfo = pchi; }
  void setChannelInfo(const ChannelInfo* pchi) { m_pchanInfo.reset(pchi); }
  void setChannelInfo(Channel a_channel,
                      Index a_fembID =badIndex(),
                      Index a_fembChannel =badIndex(),
                      Index a_channelStatus =badIndex()) {
    setChannelInfo(new ChannelInfo(a_channel, a_fembID, a_fembChannel, a_channelStatus));
  }

  // Return channel info.
  bool hasChannelInfo() const { return bool(m_pchanInfo); }
  ChannelInfoPtr getChannelInfoPtr() const { return m_pchanInfo; }
  const DuneChannelInfo& getChannelInfo() const {
    return hasChannelInfo() ? *m_pchanInfo : badChannelInfo();
  }
  Channel channel()     const { return getChannelInfo().channel; }
  Index fembID()        const { return getChannelInfo().fembID; }
  Index fembChannel()   const { return getChannelInfo().fembChannel; }
  Index channelStatus() const { return getChannelInfo().channelStatus; }

  // Return the tick offset for the raw and samples data in this channel.
  // This should allow us to align data from multiple channels and events.
  AdcLongIndex tickOffset() const {
    AdcLongIndex tickoff = triggerTick0();
    if ( tickoff == EventInfo::badLongIndex() ) tickoff = 0;
    if ( tick0 < 0 ) return tickoff - AdcLongIndex(-tick0);
    return tickoff + AdcLongIndex(tick0);
  }

  // Copy ctor.
  // Only copy event and channel data and not samples and
  // derived data.
  // Root dictionary (6.08/06) and AdcChannelDataMap (STL) require we keep copy.
  AdcChannelData(const AdcChannelData& rhs) {
    channelClock  = rhs.channelClock;
    m_peventInfo  = rhs.m_peventInfo;
    m_pchanInfo  = rhs.m_pchanInfo;
  };

  // Delete assignment.
  AdcChannelData& operator=(const AdcChannelData&) =delete;

  // Keel default ctor and move operations.
  AdcChannelData() =default;
  AdcChannelData(AdcChannelData&&) =default;
  AdcChannelData& operator=(AdcChannelData&&) =default;

  // Check if a metadata field is defined.
  bool hasMetadata(Name mname) const {
    if ( mname.substr(0,3) == "../" ) {
      if ( viewParent == nullptr ) return false;
      return viewParent->hasMetadata(mname.substr(3));
    }
    return metadata.find(mname) != metadata.end();
  }

  // Fetch metadata.
  float getMetadata(Name mname, float def =0.0) const {
    if ( mname.substr(0,3) == "../" ) {
      if ( viewParent == nullptr ) return def;
      return viewParent->getMetadata(mname.substr(3));
    }
    FloatMap::const_iterator imtd = metadata.find(mname);
    if ( imtd == metadata.end() ) return def;
    return imtd->second;
  }

  // Set metadata.
  void setMetadata(Name mname, float val) {
    metadata[mname] = val;
  }

  // Check if an attribute is defined.
  bool hasAttribute(Name mname, float def =0.0) const {
    if ( mname.substr(0,3) == "../" ) {
      if ( viewParent == nullptr ) return false;
      return viewParent->hasAttribute(mname.substr(3));
    }
    if ( mname == "run" ) return true;
    if ( mname == "subRun" ) return true;
    if ( mname == "event" ) return true;;
    if ( mname == "trigger" ) return true;;
    if ( mname == "triggerClock" ) return true;
    if ( mname == "channelClock" ) return true;
    if ( mname == "channelClockOffset" ) return true;
    if ( mname == "channel" ) return true;
    if ( mname == "fembID" ) return true;
    if ( mname == "fembChannel" ) return true;
    if ( mname == "pedestal" ) return true;
    if ( mname == "pedestalRms" ) return true;
    if ( mname == "sampleNoise" ) return true;
    if ( mname == "digitIndex" ) return true;
    if ( mname == "wireIndex" ) return true;
    if ( mname == "raw" ) return true;
    if ( mname == "samples" ) return true;
    if ( mname == "flags" ) return true;
    if ( mname == "signal" ) return true;
    if ( mname == "rois" ) return true;
    if ( mname == "dftmags" ) return true;
    if ( mname == "dftphases" ) return true;
    if ( mname == "digit" ) return true;
    if ( mname == "wire" ) return true;
    if ( mname == "metadata" ) return metadata.size();
    return hasMetadata(mname);
  }

  // Fetch any property including metadata.
  float getAttribute(Name mname, float def =0.0) const {
    if ( mname.substr(0,3) == "../" ) {
      if ( viewParent == nullptr ) return def;
      return viewParent->getAttribute(mname.substr(3));
    }
    // For basic types, return the value.
    if ( mname == "run" ) return run();
    if ( mname == "subRun" ) return subRun();
    if ( mname == "event" ) return event();
    if ( mname == "trigger" ) return trigger();
    if ( mname == "triggerClock" ) return triggerClock();  // lose precision here
    if ( mname == "channelClock" ) return channelClock;  // lose precision here
    if ( mname == "channelClockOffset" ) return channelClock - triggerClock();
    if ( mname == "channel" ) return channel();
    if ( mname == "fembID" ) return fembID();
    if ( mname == "fembChannel" ) return fembChannel();
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
  // See RealDftNormalization.h
  static AdcIndex dftNormalization() { return 22; }

  // Return the number of views.
  // This does not include the self view.
  size_t viewSize() const { return m_views.size(); }

  // Return the vector of available view names.
  // This does not include the self view.
  NameVector viewNames() const;

  // Return if a view or view path exists.
  // This includes the self view.
  bool hasView(Name vnam) const;

  // Read a view. Empty if undefined.
  // This does not include the self view.
  const View& view(Name vnam) const;
 
  // Obtain a mutable view. Entries may be added, removed or modified.
  // View is added if not already existing.
  View& updateView(Name vnam);

  // Return the number of entries for a specified view path.
  // This includes the self view.
  AdcIndex viewSize(Name vpnam) const;

  // Return the const channel data for an entry in a view path.
  // This includes the self view.
  // Returns null if the path does not exist or does not have the
  // requested entry.
  const AdcChannelData* viewEntry(Name vpnam, AdcIndex ient) const;

  // Return the mutable channel data for an entry in a view path.
  // This includes the self view.
  // Returns null if the path does not exist or does not have the
  // requested entry.
  AdcChannelData* mutableViewEntry(Name vpnam, AdcIndex ient);
  
};

//**********************************************************************

typedef std::map<AdcChannel, AdcChannelData> AdcChannelDataMap;

#ifdef __ROOTCLING__
#pragma link C++ class AdcChannelDataMap;
#endif

//**********************************************************************

inline void AdcChannelData::clear() {
  channelClock = 0;
  pedestal = badSignal();
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
  digitIndex = badIndex();
  wireIndex = badIndex();
  metadata.clear();
  m_peventInfo.reset();
  m_pchanInfo.reset();
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

AdcChannelData::NameVector AdcChannelData::viewNames() const {
  NameVector vnams;
  for ( ViewMap::value_type ivie : m_views ) {
    vnams.push_back(ivie.first);
  }
  return vnams;
}

//**********************************************************************

bool AdcChannelData::hasView(Name vpnam) const {
  if ( vpnam.size() == 0 ) return true;
  Name::size_type ipos = vpnam.find("/");
  if ( ipos == Name::npos ) return m_views.count(vpnam);
  Name vnam = vpnam.substr(0, ipos);
  if ( ! hasView(vnam) ) return false;
  Name vsnam = vpnam.substr(ipos + 1);
  for ( const AdcChannelData& acd : view(vnam) ) {
    if ( acd.hasView(vsnam) ) return true;
  }
  return false;
}

//**********************************************************************

const AdcChannelData::View& AdcChannelData::view(Name vnam) const {
  static const View empty;
  ViewMap::const_iterator ivie = m_views.find(vnam);
  if ( ivie == m_views.end() ) return empty;
  return ivie->second;
}
 
//**********************************************************************

AdcChannelData::View& AdcChannelData::updateView(Name vnam) {
  const Name myname = "AdcChannelData::updateView: ";
  if ( vnam.size() == 0 ) {
    std::cout << myname << "WARNING: View name may not be blank." << std::endl;
    static View emptyView;
    if ( emptyView.size() ) {
      std::cout << myname << "ERROR: Empty view has been modified." << std::endl;
    }
    return emptyView;
  }
  return m_views[vnam];
}

//**********************************************************************

AdcIndex AdcChannelData::viewSize(Name vpnam) const {
  if ( vpnam.size() == 0 ) return 1;
  Name::size_type ipos = vpnam.find("/");
  if ( ipos == Name::npos ) return view(vpnam).size();
  Name vnam = vpnam.substr(0, ipos);
  Name vsnam = vpnam.substr(ipos + 1);
  AdcIndex nvie = 0;
  for ( const AdcChannelData& dat : view(vnam) ) nvie += dat.viewSize(vsnam);
  return nvie;
}

//**********************************************************************

const AdcChannelData* AdcChannelData::viewEntry(Name vpnam, AdcIndex ient) const {
  if ( vpnam.size() == 0 ) {
    if ( ient == 0 ) return this;
    else return nullptr;
  }
  Name::size_type ipos = vpnam.find("/");
  if ( ipos == Name::npos ) {
    const View& myview = view(vpnam);
    if ( ient < myview.size() ) return &myview[ient];
    return nullptr;
  }
  Name vnam = vpnam.substr(0, ipos);
  Name vsnam = vpnam.substr(ipos + 1);
  AdcIndex ientRem = ient;
  for ( const AdcChannelData& dat : view(vnam) ) {
    AdcIndex nvie = dat.viewSize(vsnam);
    if ( ientRem < nvie ) return dat.viewEntry(vsnam, ientRem);
    ientRem -= nvie;
  }
  return nullptr;
}
  
//**********************************************************************

AdcChannelData* AdcChannelData::mutableViewEntry(Name vpnam, AdcIndex ient) {
  if ( vpnam.size() == 0 ) {
    if ( ient == 0 ) return this;
    else return nullptr;
  }
  Name::size_type ipos = vpnam.find("/");
  if ( ipos == Name::npos ) {
    View& myview = updateView(vpnam);
    if ( ient < myview.size() ) return &myview[ient];
    return nullptr;
  }
  Name vnam = vpnam.substr(0, ipos);
  Name vsnam = vpnam.substr(ipos + 1);
  AdcIndex ientRem = ient;
  for ( AdcChannelData& dat : updateView(vnam) ) {
    AdcIndex nvie = dat.viewSize(vsnam);
    if ( ientRem < nvie ) return dat.mutableViewEntry(vsnam, ientRem);
    ientRem -= nvie;
  }
  return nullptr;
}
  
//**********************************************************************

#endif
