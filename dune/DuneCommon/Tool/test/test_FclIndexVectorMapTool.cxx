// test_FclIndexVectorMapTool.cxx
//
// David Adams
// April 2017
//
// Test FclIndexVectorMapTool.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexVectorMapTool.h"
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
using IndexVector = std::vector<Index>;

//**********************************************************************

int test_FclIndexVectorMapTool(bool useExistingFcl =false, Index runin =0) {
  const string myname = "test_FclIndexVectorMapTool: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_FclIndexVectorMapTool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    if ( runin == 0 ) {
      fout << "tools: {" << endl;
      fout << "  mytool: {" << endl;
      fout << "    tool_type: FclIndexVectorMapTool" << endl;
      fout << "    LogLevel: 3" << endl;
      fout << "    Entries: [ [100, 2, 3], [101, 12] ]" << endl;
      fout << "  }" << endl;
      fout << "}" << endl;
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
  auto pvm = ptm->getPrivate<IndexVectorMapTool>("mytool");
  assert( pvm != nullptr );

  cout << myname << line << endl;
  cout << myname << "Check get." << endl;
  const IndexVector& vec0 = pvm->get(100);
  assert( vec0.size() == 2 );
  const IndexVector& vec1 = pvm->get(101);
  assert( vec1.size() == 1 );
  const IndexVector& vec2 = pvm->get(102);
  assert( vec2.size() == 0 );
  
  cout << myname << line << endl;
  cout << myname << "Check contains." << endl;
  assert( pvm->contains(100, 2) );
  assert( pvm->contains(100, 3) );
  assert( ! pvm->contains(100, 4) );
  assert( pvm->contains(101, 12) );
  assert( ! pvm->contains(101, 4) );
  assert( ! pvm->contains(102, 4) );

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
  return test_FclIndexVectorMapTool(useExistingFcl, run);
}

//**********************************************************************
