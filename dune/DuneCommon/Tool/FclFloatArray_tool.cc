// FclFloatArray_tool.cc

#include "FclFloatArray.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::setw;

using Name = std::string;

//**********************************************************************

FclFloatArray::FclFloatArray(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_DefaultValue(ps.get<float>("DefaultValue")),
  m_Offset(ps.get<Index>("Offset")),
  m_Label(ps.get<Name>("Label")),
  m_Unit(ps.get<Name>("Unit")),
  m_Values(ps.get<FloatVector>("Values")) {
  const Name myname = "FclFloatArray::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "        Offset: " << m_Offset << endl;
    cout << myname << "  DefaultValue: " << m_DefaultValue << endl;
    cout << myname << "         Label: " << m_Label << endl;
    cout << myname << "          Unit: " << m_Unit << endl;
    cout << myname << "          Size: " << size() << endl;
    cout << myname << "        Values: [";
    if ( m_LogLevel == 1 || size() <= 10 ) {
      Index nshow = size() <= 10 ? size() : 5;
      for ( Index ival=0; ival<nshow; ++ival ) {
        if ( ival ) cout << ", ";
        cout << values()[ival];
      }
      if ( nshow < size() ) cout << ", ...";
      cout << "]" << endl;
    } else {
      Index iwrd = 0;
      string sgap = "            ";
      for ( Index ival=0; ival<size(); ++ival ) {
        if ( ival ) cout << ",";
        if ( ++iwrd == 10 ) {
          cout << "\n" << myname << sgap << " ";
          iwrd = 0;
        }
        cout << " " << setw(8) << values()[ival];
      }
      cout << "\n" << myname << sgap << "]" << endl;
    }
  }
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(FclFloatArray)
