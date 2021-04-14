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
  string myline = myname + "-----------------------------";

  cout << myline << endl;
  cout << myname << "Checking default roi." << endl;
  Tpc2dRoi roi0;
  assert( roi0.data().size() == 0 );
  assert( roi0.sampleSize() == 0 );
  assert( roi0.channelSize() == 0 );
  assert( roi0.sampleOffset() == 0 );
  assert( roi0.channelOffset() == 0 );

  cout << myline << endl;
  cout << myname << "Checking finite roi." << endl;
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

  cout << myline << endl;
  cout << myname << "Set values." << endl;
  Float2dData::IndexArray idxs;
  Index& kcha = idxs[0];
  Index& ksam = idxs[1];
  Index ichk = 0;
  Index ntst = 0;
  for ( kcha=0; kcha<ncha; ++kcha ) {
    for ( ksam=0; ksam<nsam; ++ksam ) {
      float val = nsam*kcha + ksam + 1;
      roi1.data().setValue(idxs, val, &ichk);
      assert( ichk == 0 );
      float chkval = roi1.data().value(idxs, &ichk);
      assert( chkval == val );
      ++ntst;
    }
  }
  cout << myname << "Test count: " << ntst << endl;

  cout << myline << endl;
  cout << myname << "Retrieve values." << endl;
  ichk = 0;
  ntst = 0;
  for ( kcha=0; kcha<ncha; ++kcha ) {
    for ( ksam=0; ksam<nsam; ++ksam ) {
      float expval = nsam*kcha + ksam + 1;
      float chkval1 = roi1.data().value(idxs, &ichk);
      assert( ichk == 0 );
      assert( chkval1 == expval );
      Index icha = icha0 + kcha;
      Index isam = isam0 + ksam;
      float chkval2 = roi1.value(icha, isam);
      assert( chkval2 == expval );
      ++ntst;
    }
  }
  cout << myname << "Test count: " << ntst << endl;

  cout << myline << endl;
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
