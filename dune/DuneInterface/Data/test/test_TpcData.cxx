// test_TpcData.cxx
//
// David Adams
// January 2021
//
// Test TpcData.

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "dune/DuneInterface/Data/TpcData.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::vector;

//**********************************************************************

int test_TpcData() {
  const string myname = "test_TpcData: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = myname + "-----------------------------";

  cout << "Checking empty data." << endl;
  TpcData tpd;
  assert( tpd.getTpcData("") == &tpd );
  assert( tpd.getTpcData(".") == &tpd );
  assert( tpd.getTpcData("sub1") == nullptr );

  cout << "Add sub1" << endl;
  TpcData* pdat1 = tpd.addTpcData("sub1");
  assert( pdat1 != nullptr );
  assert( tpd.addTpcData("sub1") == nullptr );
  assert( tpd.getTpcData("sub1") == pdat1 );

  cout << "Add sub11" << endl;
  assert( tpd.getTpcData("sub1/sub11") == nullptr );
  TpcData* pdat11 = tpd.addTpcData("sub1/sub11");
  assert( pdat11 != nullptr );
  assert( tpd.addTpcData("sub1") == nullptr );
  assert( tpd.addTpcData("sub1/sub11") == nullptr );
  assert( tpd.getTpcData("sub1") == pdat1 );
  assert( tpd.getTpcData("sub1/sub11") == pdat11 );

  cout << line << endl;
  cout << myname << "All tests passed." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  if ( argc > 1 ) {
    cout << "Usage: " << argv[0] << " [ARG]" << endl;
    return 0;
  }
  return test_TpcData();
}

//**********************************************************************
