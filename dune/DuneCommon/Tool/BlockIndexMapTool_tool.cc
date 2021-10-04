// BlockIndexMapTool_tool.cc

#include "BlockIndexMapTool.h"
#include "fhiclcpp/make_ParameterSet.h"

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setw;
using Name = std::string;
using Index = BlockIndexMapTool::Index;

//**********************************************************************

BlockIndexMapTool::BlockIndexMapTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_Starts(ps.get<IndexVector>("Starts")),
  m_Values(ps.get<IndexVector>("Values")),
  m_Period(ps.get<Index>("Period")) {
  const Name myname = "BlockIndexMapTool::ctor: ";
  // Build teh map.
  m_map[0] = 0;
  for ( Index ival=0; ival<m_Starts.size(); ++ival ) {
     Index idx = m_Starts[ival];
     Index val = ival < m_Values.size() ? m_Values[ival] : 0;
     m_map[idx] = val;
  }
  // Display the parameters.
  if ( m_LogLevel >= 1 ) {
    // Use largest value to set the display width for indices and values.
    Index imax = 0;
    for ( Index val : m_Starts ) if ( val > imax ) imax = val;
    for ( Index val : m_Values ) if ( val > imax ) imax = val;
    Index wid = 1;
    Index chk = 1;
    while ( chk < imax ) {
      chk *=10;
      ++wid;
    }
    cout << myname << "   LogLevel: " << m_LogLevel << endl;
    cout << myname << "     Starts: [";
    bool first = true;
    for ( Index val : m_Starts ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << setw(wid) << val;
    }
    cout << "]" << endl;
    cout << myname << "     Values: [";
    first = true;
    for ( Index val : m_Values ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << setw(wid) << val;
    }
    cout << "]" << endl;
    cout << myname << "    Period: " << m_Period << endl;
  }
}

//**********************************************************************

Index BlockIndexMapTool::get(Index icha) const {
  const Name myname = "BlockIndexMapTool::get: ";
  Index kcha = m_Period ? icha%m_Period : icha;
  return (--m_map.upper_bound(kcha))->second;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(BlockIndexMapTool)
