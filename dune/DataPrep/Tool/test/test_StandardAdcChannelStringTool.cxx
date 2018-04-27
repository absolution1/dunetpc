// test_StandardAdcChannelStringTool.cxx
//
// David Adams
// April 2017
//
// Test StandardAdcChannelStringTool.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using fhicl::ParameterSet;

//**********************************************************************

int test_StandardAdcChannelStringTool(bool useExistingFcl =false) {
  const string myname = "test_StandardAdcChannelStringTool: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_StandardAdcChannelStringTool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: StandardAdcChannelStringTool" << endl;
    fout << "    LogLevel:     2" << endl;
    fout << "    RunWidth:     4" << endl;
    fout << "    SubRunWidth:  3" << endl;
    fout << "    EventWidth:   6" << endl;
    fout << "    ChannelWidth: 5" << endl;
    fout << "    CountWidth:   6" << endl;
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
  assert( tm.toolNames().size() == 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto past = tm.getPrivate<AdcChannelStringTool>("mytool");
  assert( past != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData acd;
  acd.run = 123;
  acd.subRun = 45;
  acd.event = 246;
  acd.channel = 1357;

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  string rawName = "data_run%RUN%-%SUBRUN%_ev%EVENT%_ch%CHAN%_%COUNT%.dat";
  string expName = "data_run0123-045_ev000246_ch01357_000023.dat";
  string outName = past->build(acd, rawName, 23);
  cout << myname << "  Raw name: " << rawName << endl;
  cout << myname << "  Out name: " << outName << endl;
  assert( outName == expName );

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
  return test_StandardAdcChannelStringTool(useExistingFcl);
}

//**********************************************************************
