// test_Result.cxx
//
// David Adams
// September 2017
//
// Test Result.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Data/Result.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

int test_Result(bool useExistingFcl =false) {
  const string myname = "test_Result: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Creating success result." << endl;
  Result res;
  assert( ! res );

  cout << myname << line << endl;
  cout << myname << "Creating failed result." << endl;
  Result resf(3);
  assert( resf );
  assert( int(resf) != 0 );
  assert( resf == 3 );

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
  return test_Result(useExistingFcl);
}

//**********************************************************************
