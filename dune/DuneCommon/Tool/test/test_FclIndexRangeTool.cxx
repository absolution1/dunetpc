// test_FclIndexRangeTool.cxx
//
// David Adams
// April 2017
//
// Test FclIndexRangeTool.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
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

int test_FclIndexRangeTool(bool useExistingFcl =false, Index runin =0) {
  const string myname = "test_FclIndexRangeTool: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_FclIndexRangeTool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    if ( runin == 0 ) {
      fout << "tools: {" << endl;
      fout << "  mytool: {" << endl;
      fout << "    tool_type: FclIndexRangeTool" << endl;
      fout << "    LogLevel: 2" << endl;
      fout << "    range1: { labels:[\"Range 1\"] begin:10 end:20}" << endl;
      fout << "    range2: { labels:[\"Range 2\"] begin:20 end:30}" << endl;
      fout << "  }" << endl;
      fout << "}" << endl;
    } else {
      fout << "#include \"dunecommon_tools.fcl\"" << endl;
      fout << "tools.mytool: @local::tools.protoduneIndexRangeTool" << endl;
    }
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
  cout << myname << "Fetching tool." << endl;
  auto irt = tm.getPrivate<IndexRangeTool>("mytool");
  assert( irt != nullptr );

  cout << myname << line << endl;
  cout << "Fetch range1." << endl;
  IndexRange ir1 = irt->get("range1");
  cout << ir1.rangeString() << endl;
  assert( ir1.isValid() );
  assert( ir1.name == "range1" );
  assert( ir1.label() == "Range 1" );
  assert( ir1.begin == 10 );
  assert( ir1.end == 20 );

  cout << myname << line << endl;
  cout << "Fetch range2." << endl;
  IndexRange ir2 = irt->get("range2");
  cout << ir2.rangeString() << endl;
  assert( ir2.isValid() );
  assert( ir2.name == "range2" );
  assert( ir2.label() == "Range 2" );
  assert( ir2.begin == 20 );
  assert( ir2.end == 30 );

  cout << myname << line << endl;
  cout << "Fetch bad range" << endl;
  IndexRange irb = irt->get("rangebad");
  cout << irb.rangeString() << endl;
  assert( ! irb.isValid() );

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
  return test_FclIndexRangeTool(useExistingFcl, run);
}

//**********************************************************************
