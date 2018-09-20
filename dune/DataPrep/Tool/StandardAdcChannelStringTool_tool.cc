// StandardAdcChannelStringTool_tool.cc

#include "StandardAdcChannelStringTool.h"
#include "dune/DuneCommon/StringManipulator.h"
#include <iostream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
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
  m_CountWidth(ps.get<Index>("CountWidth")),
  m_FembWidth(ps.get<Index>("FembWidth")) {
  const string myname = "StandardAdcChannelStringTool::ctor: ";
  m_reps[0] = "RUN";
  m_reps[1] = "SUBRUN";
  m_reps[2] = "EVENT";
  m_reps[3] = "CHAN";
  m_reps[4] = "COUNT";
  m_reps[5] = "CHAN1";
  m_reps[6] = "CHAN2";
  m_reps[7] = "FEMB";
  m_wids[0] = m_RunWidth;
  m_wids[1] = m_SubRunWidth;
  m_wids[2] = m_EventWidth;
  m_wids[3] = m_ChannelWidth;
  m_wids[4] = m_CountWidth;
  m_wids[5] = m_ChannelWidth;
  m_wids[6] = m_ChannelWidth;
  m_wids[7] = m_FembWidth;
  m_bads[0] = "RunNotFound";
  m_bads[1] = "SubRunNotFound";
  m_bads[2] = "EventNotFound";
  m_bads[3] = "ChannelNotFound";
  m_bads[4] = "CountNotFound";
  m_bads[5] = "Channel1NotFound";
  m_bads[6] = "Channel2NotFound";
  m_bads[7] = "FembNotFound";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "      RunWidth: " << m_RunWidth << endl;
    cout << myname << "   SubRunWidth: " << m_SubRunWidth << endl;
    cout << myname << "    EventWidth: " << m_EventWidth << endl;
    cout << myname << "  ChannelWidth: " << m_ChannelWidth << endl;
    cout << myname << "    CountWidth: " << m_CountWidth << endl;
    cout << myname << "     FembWidth: " << m_CountWidth << endl;
  }
}

//**********************************************************************

string StandardAdcChannelStringTool::
build(const AdcChannelData& acd, const DataMap& dm, string spat) const {
  const string myname = "StandardAdcChannelStringTool::build: ";
  // First replace the indices (run, event, ....)
  Index vals[m_nrep] = {acd.run, acd.subRun, acd.event, acd.channel,
                        Index(dm.getInt("count")),
                        Index(dm.getInt("chan1")),
                        Index(dm.getInt("chan2")),
                        acd.fembID};
  bool isBad[m_nrep] = {
    acd.run     == AdcChannelData::badIndex,
    acd.subRun  == AdcChannelData::badIndex,
    acd.event   == AdcChannelData::badIndex,
    acd.channel == AdcChannelData::badChannel,
    !dm.haveInt("count"),
    !dm.haveInt("chan1"),
    !dm.haveInt("chan2"),
    acd.fembID  == AdcChannelData::badIndex,
  };
  string sout = spat;
  for ( Index irep=0; irep<m_nrep; ++irep ) {
    string smat = m_reps[irep] + "%";
    string::size_type ipos = 1;  // Current index in string
    Index wNew = 0;  // Width of the replacement field.
    while ( (ipos = sout.find(smat, ipos)) != string::npos ) {
      // Check for replacement with natural width.
      string::size_type iposRep = ipos;  // The position where we make the replacement.
      string::size_type wOld = 0;  // Width of the field to be replaced
      if ( sout[--iposRep] == '%' ) {
        wOld = smat.size() + 1;
      // Check for replacement with padded width.
      } else if ( iposRep > 0 && sout[--iposRep] == '%' ) {
        istringstream sswid(sout.substr(iposRep+1, 1));
        wNew = 999;
        sswid >> wNew;
        if ( wNew == 0 ) wNew = m_wids[irep];
        if ( wNew == 999 ) wNew = 0;
        wOld = smat.size() + 2;
      } else {
        // If we get here we found "XXX%" but not the preceding "%'.
        // Continue to the later part of the string.
        ++ipos;
      }
      if ( wOld > 0 ) {
        string sval;  // replacement field
        if ( isBad[irep] ) {
          sval = m_bads[irep];
        } else {
          Index ival = vals[irep];
          ostringstream ssval;
          ssval << ival;
          sval = ssval.str();
          while ( wNew > sval.size() ) sval = "0" + sval;
        }
        sout.replace(iposRep, wOld, sval);
        ipos = iposRep + sval.size();
      }
    }
  }
  // Next replace the units strings.
  StringManipulator sman(sout);
  string sunit = acd.sampleUnit;
  string sunitSpaced = sunit.size() ? " " + sunit : "";
  string sunitWrapped = sunit.size() ? "(" + sunit + ")" : "";
  string sunitSpWrapped = sunit.size() ? " " + sunitWrapped : "";
  string sunitOptWrapped = sunit.find(" ") != string::npos ? sunitWrapped : "";
  string sunitSpOptWrapped = sunitOptWrapped.size() ? " " + sunitOptWrapped : "";
  string sunitBarred = sunit.size() ? "[" + sunit + "]" : "";
  string sunitSpBarred = sunit.size() ? " " + sunitBarred : "";
  sman.replace("%SUNIT%", sunit);
  sman.replace("% SUNIT%", sunitSpaced);
  sman.replace("%(SUNIT)%", sunitWrapped);
  sman.replace("% (SUNIT)%", sunitSpWrapped);
  sman.replace("%((SUNIT))%", sunitOptWrapped);
  sman.replace("% ((SUNIT))%", sunitSpOptWrapped);
  sman.replace("%[SUNIT]%", sunitBarred);
  sman.replace("% [SUNIT]%", sunitSpBarred);
  if ( m_LogLevel >= 2 ) cout << myname << spat << " --> " << sout << endl;
  return sout;
}

//**********************************************************************
