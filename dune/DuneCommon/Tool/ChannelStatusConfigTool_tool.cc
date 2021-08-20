// ChannelStatusConfigTool_tool.cc

#include "ChannelStatusConfigTool.h"
#include "dune/DuneInterface/Data/AdcTypes.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setw;
using Name = std::string;
using Index = ChannelStatusConfigTool::Index;

//**********************************************************************

ChannelStatusConfigTool::ChannelStatusConfigTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_DefaultIndex(ps.get<Index>("DefaultIndex")),
  m_IndexVectors(ps.get<IndexVectorVector>("IndexVectors")) {
  const Name myname = "ChannelStatusConfigTool::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "     LogLevel: " << m_LogLevel << endl;
    Index nidx = m_IndexVectors.size();
    Index nsho = 0;
    for ( Index iidx=0; iidx<nidx; ++iidx ) {
      const IndexVector& vals = m_IndexVectors[iidx];
      Index nval = vals.size();
      if ( nval == 0 ) continue;
      if ( nsho == 0 ) {
        cout << myname << "  Index  Count" << endl;
      }
      ++nsho;
      cout << myname << setw(7) << iidx << setw(7) << nval << endl;
    }
  }
  // Build the map.
  Index idxSize = 0;
  for ( const IndexVector& vec : m_IndexVectors ) {
    for ( Index idx : vec ) {
      if ( idx >= idxSize ) idxSize = idx + 1;
    }
  }
  m_vals.resize(idxSize, m_DefaultIndex);
  Index ivec = m_IndexVectors.size();
  IndexVector counts(idxSize, 0);
  while ( ivec ) {
    const IndexVector& vec = m_IndexVectors[--ivec];
    for ( Index idx : vec ) {
      m_vals[idx] = ivec;
      ++counts[idx];
    }
  }
  if ( m_LogLevel == 2 ) {
    Index ndup = 0;
    for ( Index idx=0; idx<idxSize; ++idx ) {
      if ( counts[idx] > 1 ) {
        if ( ndup >= 20 ) {
          cout << ", ...";
          break;
        }
        if ( ndup == 0 ) {
          cout << myname << "Found indices with duplicate entries: ";
        } else {
          cout << ", ";
        }
        cout << idx;
        ++ndup;
      }
    }
    if ( ndup == 0 ) {
      cout << myname << "No duplicate entries.";
    }
    cout << endl;
  }
}

//**********************************************************************

Index ChannelStatusConfigTool::get(Index icha) const {
  const Name myname = "ChannelStatusConfigTool::get: ";
  if ( icha < m_vals.size() ) return m_vals[icha];
  return m_DefaultIndex;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(ChannelStatusConfigTool)
