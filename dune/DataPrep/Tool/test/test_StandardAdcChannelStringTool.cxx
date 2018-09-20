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
using std::vector;
using fhicl::ParameterSet;
using Index = unsigned int;

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
    fout << "    FembWidth:    3" << endl;
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
  acd.sampleUnit = "my units";
  acd.fembID = 24;

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  vector<string> rawNames = {
    "run%RUN%",
    "srun%SUBRUN%",
    "run%RUN%_srun%XXBRUN%",
    "run%RUN%_srun%SUBRUN%",
    "data_run%RUN%-%SUBRUN%_ev%EVENT%_ch%CHAN%_%COUNT%.dat",
    "data_run%5RUN%-%5SUBRUN%_ev%5EVENT%_ch%5CHAN%_%5COUNT%.dat",
    "data_run%0RUN%-%0SUBRUN%_ev%0EVENT%_ch%0CHAN%_%0COUNT%.dat",
    "data_run%0RUN%_run%8RUN%.dat",
    "Units are %SUNIT%",
    "Units are %(SUNIT)%",
    "Units are %((SUNIT))%",
    "Units are% ((SUNIT))%",
    "FEMB is %FEMB%",
    "femb%4FEMB%",
    "femb%0FEMB%"
  };
  vector<string> expNames = {
    "run123",
    "srun45",
    "run123_srun%XXBRUN%",
    "run123_srun45",
    "data_run123-45_ev246_ch1357_23.dat",
    "data_run00123-00045_ev00246_ch01357_00023.dat",
    "data_run0123-045_ev000246_ch01357_000023.dat",
    "data_run0123_run00000123.dat",
    "Units are my units",
    "Units are (my units)",
    "Units are (my units)",
    "Units are (my units)",
    "FEMB is 24",
    "femb0024",
    "femb024"
  };
  DataMap dm;
  dm.setInt("count", 23);
  assert( dm.haveInt("count") );

  for ( Index inam=0; inam<rawNames.size(); ++inam ) {
    string rawName = rawNames[inam];
    string expName = expNames[inam];
    string outName = past->build(acd, dm, rawName);
    cout << myname << "  Raw name: " << rawName << endl;
    cout << myname << "  Exp name: " << expName << endl;
    cout << myname << "  Out name: " << outName << endl;
    assert( outName == expName );
  }

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
