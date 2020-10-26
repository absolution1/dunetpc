// ProtoduneChannelHelper.h
//
// David Adams
// October 2020
//
// Tool to intepret PDSP channel numbers.

#ifndef ProtoduneChannelHelper_H
#define ProtoduneChannelHelper_H

#include <string>

class ProtoduneChannelHelper {

public:

  using Index = unsigned int;
  using Name = std::string;

  static Index badIndex() { return -1; }

  // Helpers to extract info from a channel index.
  // Set isOff true if the index is offline.
  static Index      tpcSet(Index chan, bool isOff); // 0-5
  static Index         apa(Index chan, bool isOff); // 1-6 (PDSP-1 convention)
  static Index        femb(Index chan, bool isOff); // 1-20
  static Index        asic(Index chan, bool isOff); // 1-8
  static Index asicChannel(Index chan, bool isOff); // 0-15

  // Convert offline to online channel index.
  static Index onlineChannel(Index chanOff, Index dbg =0);

  // Ctor specifying whether channel indices are online or offline.
  ProtoduneChannelHelper(bool isOff);

  bool isOffline() const { return m_isOff; }

  // Return channel component indices.
  Index      tpcSet(Index chan) const { return      tpcSet(chan, isOffline()); }
  Index         apa(Index chan) const { return         apa(chan, isOffline()); }
  Index        femb(Index chan) const { return        femb(chan, isOffline()); }
  Index        asic(Index chan) const { return        asic(chan, isOffline()); }
  Index asicChannel(Index chan) const { return asicChannel(chan, isOffline()); }

  // Return the ASIC location: AFF-CXX
  // A = APA, FF=FEMB, C=ascic (1-8), XX=asic channel (0-15)
  Name asicChannelName(Index icha) const;

private:

  bool m_isOff;

};

#endif
