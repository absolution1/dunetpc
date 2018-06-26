// TimingRawDecoderOffsetTool_tool.cc

#include "TimingRawDecoderOffsetTool.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

using Offset = TimingRawDecoderOffsetTool::Offset;

//**********************************************************************

TimingRawDecoderOffsetTool::TimingRawDecoderOffsetTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_TpcTickPhase(ps.get<Index>("TpcTickPhase")),
  m_Unit(ps.get<Name>("Unit")) {
  const Name myname = "TimingRawDecoderOffsetTool::ctor: ";
  if ( m_LogLevel ) {
    cout << "Configuration:" << endl;
    cout << "      LogLevel: " << m_LogLevel << endl;
    cout << "  TpcTickPhase: " << m_TpcTickPhase << endl;
    cout << "          Unit: " << m_Unit << endl;
  }
}

//**********************************************************************

Offset TimingRawDecoderOffsetTool::offset(const Data& dat) const {
  const Name myname = "TimingRawDecoderOffsetTool::offset: ";
  Name ifname = "artdaqTimestamp-Run" + std::to_string(dat.run) + 
                "-Event" + std::to_string(dat.event) + ".dat";
  Offset res;
  ifstream fin(ifname.c_str());
  if ( ! fin ) {
    cout << myname << "Unable to find time offset file: " << ifname << endl;
    return res.setStatus(1);
  }
  long daqVal;
  fin >> daqVal;
  if ( m_Unit == "daq" ) {
    res.value = daqVal;
  } else if ( m_Unit == "ns" ) {
    res.value = 20*daqVal;
  } else if ( m_Unit == "tick" ) {
    res.value = (daqVal + m_TpcTickPhase)/25;
  } else {
    cout << myname << "Invalid unit: " << m_Unit << ifname << endl;
    return res.setStatus(2);
  }
  res.unit = m_Unit;
  if ( m_LogLevel >= 2 ) cout << myname << "Offset is " << res.value
                              << " " << res.unit << endl;
  return res;
}

//**********************************************************************
