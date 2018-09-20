// test_TimingRawDecoderOffsetTool.cxx
//
// David Adams
// May 2018
//
// Test TimingRawDecoderOffsetTool.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/TimeOffsetTool.h"
#include "TH1F.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using fhicl::ParameterSet;
using Index = unsigned int;

//**********************************************************************

int test_TimingRawDecoderOffsetTool(bool useExistingFcl =false) {
  const string myname = "test_TimingRawDecoderOffsetTool: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_TimingRawDecoderOffsetTool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: TimingRawDecoderOffsetTool" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "    TpcTickPhase: 0" << endl;
    fout << "    Unit: \"tick\"" << endl;
    fout << "    RunDataTool: \"\"" << endl;
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
  assert( tm.toolNames().size() >= 1 );

  cout << myname << line << endl;
  cout << myname << "Create time offset file." << endl;
  Index run = 123;
  Index evt = 2468;
  Index daqVal = 369258;

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto tot = tm.getPrivate<TimeOffsetTool>("mytool");
  assert( tot != nullptr );

  cout << myname << line << endl;
  cout << "Fetch the time offset in ticks." << endl;
  TimeOffsetTool::Data dat;
  dat.run = run;
  dat.event = evt;
  dat.triggerClock = daqVal;
  TimeOffsetTool::Offset off = tot->offset(dat);
  assert( off.isValid() );
  assert( off.value == daqVal/25 );
  assert( off.unit == "tick" );

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
      cout << "Usage: " << argv[0] << " [keepFCL] [RUN]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      cout << "  If RUN is nonzero, the data for that run are displayed." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_TimingRawDecoderOffsetTool(useExistingFcl);
}

//**********************************************************************
