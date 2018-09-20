// test_FclRunDataTool.cxx
//
// David Adams
// April 2017
//
// Test FclRunDataTool.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
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

int test_FclRunDataTool(bool useExistingFcl =false, Index runin =0) {
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
    if ( runin == 0 ) {
      fout << "tools: {" << endl;
      fout << "  mytool: {" << endl;
      fout << "    tool_type: FclRunDataTool" << endl;
      fout << "    LogLevel: 3" << endl;
      fout << "    FileNames: [\"rdtest/rundata.fcl\", \"rdtest/rundata/run000123.fcl\"]" << endl;
      fout << "  }" << endl;
      fout << "}" << endl;
    } else {
      fout << "#include \"dunecommon_tools.fcl\"" << endl;
      fout << "tools.mytool: @local::tools.protoduneRunDataTool" << endl;
    }
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
  assert( tm.toolNames().size() >= 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto rdt = tm.getPrivate<RunDataTool>("mytool");
  assert( rdt != nullptr );

  Index run = runin == 0 ? 123 : runin;

  cout << myname << line << endl;
  cout << "Fetch run data." << endl;
  RunData rdat = rdt->runData(run);
  cout << rdat << endl;
  assert( rdat.run() == run );
  if ( run == 123 ) {
    assert( rdat.cryostat() == "protodune" );
    assert( rdat.apas().size() == 6 );
    assert( rdat.gain() == 25.0 );
    assert( rdat.shaping() == 3.0 );
    assert( rdat.leakage() == 500 );
    assert( rdat.hvfrac() == 0.5 );
    assert( rdat.pulserAmplitude() == 8 );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  Index run = 0;
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
  if ( argc > 2 ) {
    string sarg(argv[2]);
    istringstream ssarg(argv[2]);
    ssarg >> run;
  }
  return test_FclRunDataTool(useExistingFcl, run);
}

//**********************************************************************
