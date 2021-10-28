// FclIndexVectorMapTool_tool.cc

#include "FclIndexVectorMapTool.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::setw;
using Name = std::string;

//**********************************************************************

FclIndexVectorMapTool::FclIndexVectorMapTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_Entries(ps.get<IndexVectorVector>("Entries")) {
  const Name myname = "FclIndexVectorMapTool::ctor: ";
  for ( const IndexVector ent : m_Entries ) {
    if ( ent.size() == 0 ) continue;
    Index ient = ent.front();
    IndexVector& vals = m_entryMap[ient];
    vals.insert(vals.end(), ++ent.begin(), ent.end());
  }
  Index vecCount = m_entryMap.size();
  Index valCount = 0;
  for ( const IndexVectorMap::value_type& ent : m_entryMap ) valCount += ent.second.size();
  if ( m_LogLevel >= 1 ) {
    cout << myname << "  LogLevel: " << m_LogLevel << endl;
    cout << myname << "   Entries: " << vecCount << " vectors with "
         << valCount << " values";
    if ( m_LogLevel > 1 ) {
      cout << ":" << endl;
      for ( const IndexVectorMap::value_type& ent : m_entryMap ) {
        cout << myname << "  " << ent.first << ": [";
        bool first = true;
        for ( Index val : ent.second ) {
          if ( first ) first = false;
          else cout << ", ";
          cout << val;
        }
        cout << "]" << endl;
      }
    } else {
      cout << "." << endl;
    }
  }
}

//**********************************************************************

IndexVectorMapTool::IndexVector FclIndexVectorMapTool::get(Index idx) const {
  const Name myname = "FclIndexVectorMapTool::get: ";
  IndexVectorMap::const_iterator ient = m_entryMap.find(idx);
  static const IndexVector empty;
  if ( ient == m_entryMap.end() ) return empty;
  return ient->second;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(FclIndexVectorMapTool)
