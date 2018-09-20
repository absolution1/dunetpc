// FixedTimeOffsetTool_tool.cc

#include "FixedTimeOffsetTool.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

using Offset = FixedTimeOffsetTool::Offset;

//**********************************************************************

FixedTimeOffsetTool::FixedTimeOffsetTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_Value(ps.get<long>("Value")),
  m_Rem(ps.get<double>("Rem")),
  m_Unit(ps.get<Name>("Unit")) {
  const Name myname = "FixedTimeOffsetTool::ctor: ";
  if ( m_LogLevel ) {
    cout << "Configuration:" << endl;
    cout << "  LogLevel: " << m_LogLevel << endl;
    cout << "     Value: " << m_Value << endl;
    cout << "       Rem: " << m_Rem << endl;
    cout << "      Unit: " << m_Unit << endl;
  }
}

//**********************************************************************

Offset FixedTimeOffsetTool::offset(const Data& dat) const {
  const Name myname = "FixedTimeOffsetTool::offset: ";
  Name ifname = "artdaqTimestamp-Run" + std::to_string(dat.run) + 
                "-Event" + std::to_string(dat.event) + ".dat";
  Offset res;
  res.value = m_Value;
  res.rem = m_Rem;
  res.unit = m_Unit;
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Offset is " << res.value;
    if ( res.rem ) cout << " + " << res.rem;
    cout << " " << res.unit << endl;
  }
  return res;
}

//**********************************************************************
