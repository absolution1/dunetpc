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
  m_wids[0] = m_RunWidth;
  m_wids[1] = m_SubRunWidth;
  m_wids[2] = m_EventWidth;
  m_wids[3] = m_ChannelWidth;
  m_wids[4] = m_CountWidth;
  m_bads[0] = "RunNotFound";
  m_bads[1] = "SubRunNotFound";
  m_bads[2] = "EventNotFound";
  m_bads[3] = "ChannelNotFound";
  m_bads[4] = "CountNotFound";
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
build(const AdcChannelData& acd, string spat, Index count) const {
  const string myname = "StandardAdcChannelStringTool::build: ";
  Index vals[m_nrep] = {acd.run, acd.subRun, acd.event, acd.channel, count};
  bool isBad[m_nrep] = {
    acd.run     == AdcChannelData::badIndex,
    acd.subRun  == AdcChannelData::badIndex,
    acd.event   == AdcChannelData::badIndex,
    acd.channel == AdcChannelData::badChannel,
    count       == AdcChannelData::badIndex
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
  if ( m_LogLevel >= 2 ) cout << myname << spat << " --> " << sout << endl;
  return sout;
}

//**********************************************************************
