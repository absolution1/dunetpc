// TestTool_tool.cc

#include "TestTool.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setw;
using std::fixed;
using std::setprecision;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

TestTool::TestTool(fhicl::ParameterSet const& ps)
: m_Label(ps.get<string>("Label")),
  m_count(0) {
  const string myname = "TestTool::ctor: ";
  cout << myname << "Creating \"" << label() << "\" at @" << this << endl;
}

//**********************************************************************

TestTool::~TestTool() {
  const string myname = "TestTool::dtor: ";
  cout << myname << "Deleting \"" << label() << "\" at @" << this << endl;
}

//**********************************************************************
