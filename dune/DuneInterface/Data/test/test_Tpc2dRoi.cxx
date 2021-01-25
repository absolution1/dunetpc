// test_Tpc2dRoi.cxx
//
// David Adams
// January 2021
//
// Test Tpc2dRoi.

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "dune/DuneInterface/Data/Tpc2dRoi.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::vector;

using Index = Tpc2dRoi::Index;

//**********************************************************************

int test_Tpc2dRoi() {
  const string myname = "test_Tpc2dRoi: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = myname + "-----------------------------";

  cout << line << endl;
  cout << "Checking default roi." << endl;
  Tpc2dRoi roi0;
  assert( roi0.data().size() == 0 );
  assert( roi0.sampleSize() == 0 );
  assert( roi0.channelSize() == 0 );
  assert( roi0.sampleOffset() == 0 );
  assert( roi0.channelOffset() == 0 );

  cout << line << endl;
  cout << "Checking finite roi." << endl;
  Index nsam = 100;
  Index ncha = 20;
  Index ndat = ncha*nsam;
  Index isam0 = 10;
  Index icha0 = 500;
  Tpc2dRoi roi1(20, 100, 500, 10);
  assert( roi1.data().size() == ndat );
  assert( roi1.sampleSize() == nsam );
  assert( roi1.channelSize() == ncha );
  assert( roi1.sampleOffset() == isam0 );
  assert( roi1.channelOffset() == icha0 );

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
  return test_Tpc2dRoi();
}

//**********************************************************************
