// test_FclRunDataTool.cxx
//
// David Adams
// April 2017
//
// Test FclRunDataTool.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "TH1F.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using fhicl::ParameterSet;

//**********************************************************************

int test_FclRunDataTool(bool useExistingFcl =false) {
  const string myname = "test_FclRunDataTool: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_FclRunDataTool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: FclRunDataTool" << endl;
    fout << "    LogLevel: 3" << endl;
    fout << "    FileNames: [\"rdtest/rundata.fcl\", \"rdtest/rundata/run000123.fcl\"]" << endl;
    fout << "    RunPrefix: \"\"" << endl;
    fout << "    SubPrefix: \"\"" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Create rundata fcl" << endl;
  if ( ! useExistingFcl ) {
    system("rm -rf rdtest");
    system("mkdir rdtest");
    ofstream fout("rdtest/rundata.fcl");
    fout << "cryostat: protodune\n";
    fout << "apas: [1, 2, 3, 4, 5, 6]\n";
    fout << "gain: 14.0\n";
    fout << "shaping: 2.0\n";
    fout << "leakage: 500\n";
    fout << "hvfrac: 1.0\n";
    fout << "pulserAmplitude: 0\n";
    fout.close();
  }
  if ( ! useExistingFcl ) {
    system("mkdir rdtest/rundata");
    ofstream fout("rdtest/rundata/run000123.fcl");
    fout << "run: 123\n";
    fout << "gain: 25.0\n";
    fout << "shaping: 3.0\n";
    fout << "hvfrac: 0.5\n";
    fout << "pulserAmplitude: 8\n";
    fout.close();
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
  auto rdt = tm.getPrivate<RunDataTool>("mytool");
  assert( rdt != nullptr );

  cout << myname << line << endl;
  cout << "Fetch run data." << endl;
  RunData rdat = rdt->runData(123);
  cout << rdat << endl;
  assert( rdat.run() == 123 );
  assert( rdat.cryostat() == "protodune" );
  assert( rdat.apas().size() == 6 );
  assert( rdat.gain() == 25.0 );
  assert( rdat.shaping() == 3.0 );
  assert( rdat.leakage() == 500 );
  assert( rdat.hvfrac() == 0.5 );
  assert( rdat.pulserAmplitude() == 8 );

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
  return test_FclRunDataTool(useExistingFcl);
}

//**********************************************************************
