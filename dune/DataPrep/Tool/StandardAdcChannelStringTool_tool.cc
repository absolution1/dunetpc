// StandardAdcChannelStringTool_tool.cc

#include "StandardAdcChannelStringTool.h"
#include "dune/DuneCommon/StringManipulator.h"
#include <iostream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;

//**********************************************************************
// Class methods.
//**********************************************************************

StandardAdcChannelStringTool::
StandardAdcChannelStringTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_RunWidth(ps.get<Index>("RunWidth")),
  m_SubRunWidth(ps.get<Index>("SubRunWidth")),
  m_EventWidth(ps.get<Index>("EventWidth")),
  m_ChannelWidth(ps.get<Index>("ChannelWidth")),
  m_CountWidth(ps.get<Index>("CountWidth")) {
  const string myname = "StandardAdcChannelStringTool::ctor: ";
  m_reps[0] = "RUN";
  m_reps[1] = "SUBRUN";
  m_reps[2] = "EVENT";
  m_reps[3] = "CHAN";
  m_reps[4] = "COUNT";
  m_reps[5] = "CHAN1";
  m_reps[6] = "CHAN2";
  m_wids[0] = m_RunWidth;
  m_wids[1] = m_SubRunWidth;
  m_wids[2] = m_EventWidth;
  m_wids[3] = m_ChannelWidth;
  m_wids[4] = m_CountWidth;
  m_wids[5] = m_ChannelWidth;
  m_wids[6] = m_ChannelWidth;
  m_bads[0] = "RunNotFound";
  m_bads[1] = "SubRunNotFound";
  m_bads[2] = "EventNotFound";
  m_bads[3] = "ChannelNotFound";
  m_bads[4] = "CountNotFound";
  m_bads[5] = "Channel1NotFound";
  m_bads[6] = "Channel2NotFound";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "      RunWidth: " << m_RunWidth << endl;
    cout << myname << "   SubRunWidth: " << m_SubRunWidth << endl;
    cout << myname << "    EventWidth: " << m_EventWidth << endl;
    cout << myname << "  ChannelWidth: " << m_ChannelWidth << endl;
    cout << myname << "    CountWidth: " << m_CountWidth << endl;
  }
}

//**********************************************************************

string StandardAdcChannelStringTool::
build(const AdcChannelData& acd, const DataMap& dm, string spat) const {
  const string myname = "StandardAdcChannelStringTool::build: ";
  Index vals[m_nrep] = {acd.run, acd.subRun, acd.event, acd.channel,
                        Index(dm.getInt("count")),
                        Index(dm.getInt("chan1")),
                        Index(dm.getInt("chan2"))};
  bool isBad[m_nrep] = {
    acd.run     == AdcChannelData::badIndex,
    acd.subRun  == AdcChannelData::badIndex,
    acd.event   == AdcChannelData::badIndex,
    acd.channel == AdcChannelData::badChannel,
    !dm.haveInt("count"),
    !dm.haveInt("chan1"),
    !dm.haveInt("chan2")
  };
  string sout = spat;
  StringManipulator sman(sout);
  for ( Index irep=0; irep<m_nrep; ++irep ) {
    string srep = "%" + m_reps[irep] + "%";
    Index w = m_wids[irep];
    Index val = vals[irep];
    if ( isBad[irep] ) {
      sman.replace(srep, m_bads[irep]);
    } else if ( w == 0 ) {
      sman.replace(srep, val);
    } else {
      sman.replaceFixedWidth(srep, val, w);
    }
  }
  string sunit = acd.sampleUnit;
  string sunitSpaced = sunit.size() ? " " + sunit : "";
  string sunitWrapped = sunit.size() ? "(" + sunit + ")" : "";
  string sunitSpWrapped = sunit.size() ? " " + sunitWrapped : "";
  string sunitBarred = sunit.size() ? "[" + sunit + "]" : "";
  string sunitSpBarred = sunit.size() ? " " + sunitBarred : "";
  sman.replace("%SUNIT%", sunit);
  sman.replace("% SUNIT%", sunitSpaced);
  sman.replace("%(SUNIT)%", sunitWrapped);
  sman.replace("% (SUNIT)%", sunitSpWrapped);
  sman.replace("%[SUNIT]%", sunitBarred);
  sman.replace("% [SUNIT]%", sunitSpBarred);
  if ( m_LogLevel >= 2 ) cout << myname << spat << " --> " << sout << endl;
  return sout;
}

//**********************************************************************
