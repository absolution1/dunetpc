// ProtoduneChannelHelper.h
//
// David Adams
// October 2020
//
// Tool to intepret PDSP channel numbers.

#ifndef ProtoduneChannelHelper_H
#define ProtoduneChannelHelper_H

class ProtoduneChannelHelper {

public:

  using Index = unsigned int;

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

  Index      tpcSet(Index chan) const { return      tpcSet(chan, m_isOff); }
  Index         apa(Index chan) const { return         apa(chan, m_isOff); }
  Index        femb(Index chan) const { return        femb(chan, m_isOff); }
  Index        asic(Index chan) const { return        asic(chan, m_isOff); }
  Index asicChannel(Index chan) const { return asicChannel(chan, m_isOff); }

private:

  bool m_isOff;

};

#endif
