// FclIndexRangeTool_tool.cc

#include "FclIndexRangeTool.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::setw;
using NameVector = IndexRange::NameVector;

//**********************************************************************

FclIndexRangeTool::FclIndexRangeTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")) {
  const Name myname = "FclIndexRangeTool::ctor: ";
  Index namSize = 0;
  Index labSize = 0;
  for ( string rnam : ps.get_pset_names() ) {
    if ( rnam == "LogLevel" ) continue;
    cout << myname << "  " << rnam << endl;
    fhicl::ParameterSet psr = ps.get<fhicl::ParameterSet>(rnam);
    if ( m_Ranges.find(rnam) != m_Ranges.end() ) {
      cout << myname << "Overwriting duplicate range: " << rnam << endl;
    }
    IndexRange& ran = m_Ranges[rnam];
    ran.name = rnam;
    ran.labels = psr.get<NameVector>("labels");
    ran.begin = psr.get<Index>("begin");
    ran.end = psr.get<Index>("end");
    if ( ran.name.size() > namSize ) namSize = ran.name.size();
    if ( ran.label().size() > labSize ) labSize = ran.label().size();
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "     LogLevel: " << m_LogLevel << endl;
    cout << myname << "  Range count: " << m_Ranges.size() << endl;
    if ( m_LogLevel >= 2 ) {
      for ( IndexRangeMap::value_type iran : m_Ranges ) {
        const IndexRange& ran = iran.second;
        cout << myname << setw(namSize+2) << ran.name
             << setw(6) << ran.begin
             << setw(6) << ran.end
             << setw(labSize+1) << ran.label() << endl;
      }
    }
  }
}

//**********************************************************************

IndexRange FclIndexRangeTool::get(Name nam) const {
  const Name myname = "FclIndexRangeTool::runData: ";
  IndexRangeMap::const_iterator iran = m_Ranges.find(nam);
  if ( iran == m_Ranges.end() ) return IndexRange();
  return iran->second;
}

//**********************************************************************
