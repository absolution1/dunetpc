// DuneChannelInfo.h

// David Adams
// June 2019
//
// Channel info for dataprep in DUNE.

#ifndef DuneChannelInfo_H
#define DuneChannelInfo_H

#include "dune/DuneInterface/Data/AdcTypes.h"
#include <sstream>
#include <memory>

class DuneChannelInfo {

public:

  using Channel = AdcChannel;
  using Index = AdcIndex;

  Channel channel =badChannel();
  Index fembID =badIndex();
  Index fembChannel =badIndex();
  Index channelStatus =badIndex();

  static Index badChannel() { return -1u; }
  static Index badIndex() { return -1u; }
  static const DuneChannelInfo& badChannelInfo() {
    static DuneChannelInfo bci;
    return bci;
  }

  // Reset to no channel.
  void clear() {
    channel = badChannel();
    fembID = badIndex();
    fembChannel = badIndex();
    channelStatus = badIndex();
  }

  // Default ctor.
  DuneChannelInfo() = default;

  // Ctor from components.
  DuneChannelInfo(Index a_channel, Index a_fembID, Index a_fembChannel, Index a_channelStatus)
  : channel(a_channel), fembID(a_fembID), fembChannel(a_fembChannel),
    channelStatus(a_channelStatus) { }

  // Is this valid info?
  bool isValid() const { return channel != badIndex(); }

  // Order two channel info objects.
  bool operator<(const DuneChannelInfo& rhs) {
    if ( channel < rhs.channel ) return true;
    return false;
  }

  // Compare two channel info objects.
  bool operator==(const DuneChannelInfo& rhs) const {
    if ( channel       != rhs.channel ) return false;
    if ( fembID        != rhs.fembID ) return false;
    if ( fembChannel   != rhs.fembChannel ) return false;
    if ( channelStatus != rhs.channelStatus ) return false;
    return true;
  }

  bool operator!=(const DuneChannelInfo& rhs) const {
    return ! (*this == rhs);
  }

  // Return the channel number as a string.
  std::string channelString() const {
    std::ostringstream ssout;
    ssout << channel;
    return ssout.str();
  }

  // Return femb as ID-chan.
  //   opt = 0:  Return RRR.
  //         1:  Return RRR-SSS
  //         2:  if SSS==0, return  RRR, otherwise RRR-SSS.
  std::string fembString(Index wid =2, Index wch =3) const {
    std::string sid = std::to_string(fembID);
    while ( sid.size() < wid ) sid = "0" + sid;
    std::string sch = std::to_string(fembChannel);
    while ( sch.size() < wch ) sch = "0" + sch;
    return sid + "-" + sch;
  }

};

#endif
