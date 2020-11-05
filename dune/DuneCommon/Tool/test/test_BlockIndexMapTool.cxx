// test_BlockIndexMapTool.cxx
//
// David Adams
// April 2020
//
// Test BlockIndexMapTool.

#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "TH1F.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::setw;
using std::map;
using fhicl::ParameterSet;
using Index = unsigned int;

//**********************************************************************

int test_BlockIndexMapTool(bool useExistingFcl =false) {
  const string myname = "test_BlockIndexMapTool: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_BlockIndexMapTool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: BlockIndexMapTool" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    Starts: [20, 10, 30]" << endl;
    fout << "    Values: [ 2,  1,  3]" << endl;
    fout << "    Period: 100" << endl;
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
  cout << myname << "Fetching tool." << endl;
  auto ptoo = tm.getPrivate<IndexMapTool>("mytool");
  assert( ptoo != nullptr );

  cout << myname << line << endl;
  map<Index, Index> chk;
  for ( Index idx=  0; idx< 10; ++idx ) chk[idx] = 0;
  for ( Index idx= 10; idx< 20; ++idx ) chk[idx] = 1;
  for ( Index idx= 20; idx< 30; ++idx ) chk[idx] = 2;
  for ( Index idx= 30; idx<100; ++idx ) chk[idx] = 3;
  for ( Index idx=100; idx<110; ++idx ) chk[idx] = 0;
  for ( Index idx=110; idx<120; ++idx ) chk[idx] = 1;
  for ( Index idx=120; idx<130; ++idx ) chk[idx] = 2;
  for ( Index idx=130; idx<200; ++idx ) chk[idx] = 3;
  for ( auto kchk : chk ) {
    Index idx = kchk.first;
    Index chkval = kchk.second;
    Index val = ptoo->get(idx);
    cout << myname << setw(10) << idx << ": " << setw(5) << val << endl;
    assert( val == chkval );
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
      cout << "Usage: " << argv[0] << " [keepFCL]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_BlockIndexMapTool(useExistingFcl);
}

//**********************************************************************
