// test_FclFloatArray.cxx
//
// David Adams
// January 2019
//
// Test FclFloatArray.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/FloatArrayTool.h"
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

bool isequal(float lhs, float rhs ) {
  string myname = "isequal: ";
  if ( fabs(lhs - rhs) < 1.e-10 ) return true;
  cout << myname << lhs << " != " << rhs << endl;
  return false;
}

int test_FclFloatArray(bool useExistingFcl =false) {
  const string myname = "test_FclFloatArray: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_FclFloatArray.fcl";
  Index ioff = 10;
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "   tool_type: FclFloatArray" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "      Offset: " << ioff << endl;
    fout << "       Label: myvals" << endl;
    fout << "      Values: [   0, 1.01, 2.02, 3.03, 4.04, 5.05, 6.06, 7.07, 8.08, 9.09," << endl;
    fout << "                10.10, 11.11, 12.12, 13.13, 14.14, 15.15 ]" << endl;
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
  auto pto = tm.getPrivate<FloatArrayTool>("mytool");
  assert( pto != nullptr );

  cout << myname << line << endl;
  cout << "Check tool params." << endl;
  Index nval = 16;
  assert( pto->offset() == ioff );
  assert( pto->label() == "myvals" );
  assert( pto->size() == nval );
  assert( pto->offset() == ioff );

  cout << myname << line << endl;
  cout << "Check fill." << endl;
  FloatArrayTool::FloatVector fvals;
  float fval = 1234.5;
  assert( pto->fill(fvals, fval) == 0 );
  assert( fvals.size() == nval + ioff );
  for ( Index ival=0; ival<ioff; ++ival ) {
    assert( fvals[ival] == fval );
  }
  
  cout << myname << line << endl;
  cout << "Check values." << endl;
  for ( Index ival=0; ival<nval; ++ival ) {
    assert( pto->values()[ival] == fvals[ioff+ival] );
    assert( isequal(pto->values()[ival], 1.01*ival) );
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
      cout << "Usage: " << argv[0] << " [keepFCL] [RUN]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      cout << "  If RUN is nonzero, the data for that run are displayed." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_FclFloatArray(useExistingFcl);
}

//**********************************************************************
