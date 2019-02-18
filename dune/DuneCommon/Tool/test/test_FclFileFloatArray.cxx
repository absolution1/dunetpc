// test_FclFileFloatArray.cxx
//
// David Adams
// April 2017
//
// Test FclFileFloatArray.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
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
using std::setw;
using fhicl::ParameterSet;
using Index = unsigned int;

//**********************************************************************

int test_FclFileFloatArray(bool useExistingFcl =false) {
  const string myname = "test_FclFileFloatArray: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  if ( ! useExistingFcl ) {
    cout << myname << "Creating float array FCL." << endl;
    {
      ofstream fout("test_FclFileFloatArray_arr1.fcl");
      fout << "       tool_type: FclFloatArray" << endl;
      fout << "        LogLevel: 2" << endl;
      fout << "    DefaultValue: -1" << endl;
      fout << "          Offset: 10" << endl;
      fout << "           Label: arr1" << endl;
      fout << "            Unit: cm" << endl;
      fout << "          Values: [   110, 111, 112, 113, 114, 115, 116, 117, 118, 119 ]" << endl;
      fout.close();
    }
    {
      ofstream fout("test_FclFileFloatArray_arr3.fcl");
      fout << "       tool_type: FclFloatArray" << endl;
      fout << "        LogLevel: 2" << endl;
      fout << "    DefaultValue: -1" << endl;
      fout << "          Offset: 30" << endl;
      fout << "           Label: arr3" << endl;
      fout << "            Unit: cm" << endl;
      fout << "          Values: [   130, 131, 132, 133, 134, 135, 136, 137, 138, 139 ]" << endl;
      fout.close();
    }
  }
  string fclfile = "test_FclFileFloatArray.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: FclFileFloatArray" << endl;
    fout << "    LogLevel: 3" << endl;
    fout << "    FileNames: [\"test_FclFileFloatArray_arr1.fcl\", \"test_FclFileFloatArray_arr3.fcl\"]" << endl;
    fout << "    DefaultValue: 999" << endl;
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
  auto ptoo = tm.getPrivate<FloatArrayTool>("mytool");
  assert( ptoo != nullptr );

  cout << myname << line << endl;
  cout << "size =  " << ptoo->size() << endl;
  assert( ptoo->size() == 40 );

  cout << myname << line << endl;
  cout << "Unit =  " << ptoo->unit() << endl;
  assert( ptoo->unit() == "cm" );

  cout << myname << line << endl;
  cout << myname << "Check values." << endl;
  FloatArrayTool::FloatVector exp(40, 999);
  for ( int ival=10; ival<20; ++ival ) exp[ival] = 100 + ival;
  for ( int ival=30; ival<40; ++ival ) exp[ival] = 100 + ival;
  for ( Index ival=0; ival<ptoo->size(); ++ival ) {
    cout << setw(4) << ival << ": " << ptoo->value(ival) << endl;
    assert( fabs(ptoo->value(ival) - exp[ival]) < 1.e-6 );
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
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_FclFileFloatArray(useExistingFcl);
}

//**********************************************************************
