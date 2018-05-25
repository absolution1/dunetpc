// test_AdcResultDumper.cxx
//
// David Adams
// April 2018
//
// Test AdcResultDumper.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;

//**********************************************************************

int test_AdcResultDumper(bool useExistingFcl =false) {
  const string myname = "test_AdcResultDumper: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcResultDumper.fcl";
  string gname = "protodune_geo";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool1: {" << endl;
    fout << "    tool_type: AdcResultDumper" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "    Tool: \"NoSuchTool\"" << endl;
    fout << "  }" << endl;
    fout << "  mytool2: {" << endl;
    fout << "    tool_type: AdcResultDumper" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "    Tool: \"mytool1\"" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetching tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance(fclfile);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  tm.print();
  assert( tm.toolNames().size() == 2 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool 1." << endl;
  auto pard1 = tm.getPrivate<AdcChannelTool>("mytool1");
  assert( pard1 != nullptr );

  cout << myname << line << endl;
  cout << myname << "Fetching tool 2." << endl;
  auto pard2 = tm.getPrivate<AdcChannelTool>("mytool2");
  assert( pard2 != nullptr );

  cout << myname << line << endl;
  cout << myname << "Calling tool 1." << endl;
  AdcChannelData acd;
  DataMap dm1 = pard1->view(acd);
  dm1.print();
  assert( dm1.status() == 101 );

  cout << myname << line << endl;
  cout << myname << "Calling tool 2." << endl;
  DataMap dm2 = pard2->view(acd);
  dm2.print();
  assert( dm2.status() == 101 );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_AdcResultDumper(useExistingFcl);
}

//**********************************************************************
