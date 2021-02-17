// test_FftwReal2dDftData.cxx
//
// David Adams
// February 2021
//
// Test FftwReal2dDftData.

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "dune/DuneInterface/Data/FftwReal2dDftData.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::vector;

using DftData = FftwDouble2dDftData;
using Index = DftData::Index;
using IndexArray = DftData::IndexArray;
using IndexArrayVector = DftData::IndexArrayVector;
using Complex = DftData::Complex;
using Float = DftData::Float;

Index realval(Index idat) { return 2*idat + 1; }
Index imagval(Index idat) { return 2*idat + 2; }

//**********************************************************************

int test_FftwReal2dDftData(Index nsam0, Index nsam1) {
  const string myname = "test_FftwReal2dDftData: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = myname + "-----------------------------";

  cout << line << endl;
  cout << myname << "Sample sizes: " << nsam0 << ", " << nsam1 << endl;

  cout << line << endl;
  cout << myname << "Checking index array." << endl;
  IndexArray nsams({nsam0, nsam1});
  assert( nsams[0] == nsam0 );
  assert( nsams[1] == nsam1 );
  Index nComplexExp = nsam0*(nsam1/2 + 1);
  Index ndftComplex = DftData::dftComplexDataSize(nsams);
  cout << myname << "  Complex size: " << ndftComplex << endl;
  Index ndftFloat = DftData::dftFloatDataSize(nsams);
  cout << myname << "    Float size: " << ndftFloat << endl;
  assert( ndftComplex == nComplexExp );
  assert( ndftFloat == 2*nComplexExp );

  cout << line << endl;
  cout << myname << "Construct DFT object." << endl;
  DftData::Norm norm(22);
  DftData dft(norm, nsams);
  assert( dft.rank() == 2 );
  assert( dft.normalization() == norm );
  assert( dft.size() == ndftComplex );
  assert( dft.nSamples(0) == nsam0 );
  assert( dft.size(0) == nsam0 );
  assert( dft.nSamples(1) == nsam1 );
  assert( dft.size(1) == nsam1 );

  cout << myname << "Check indexing." << endl;
  Index cntij = 0;
  for ( Index idat=0; idat<nComplexExp; ++idat ) {
    IndexArrayVector idxas = dft.indexArrays(idat);
    cout << myname << setw(3) << idat << ": ";
    Index nerr = idxas.size() == 0;
    bool first = true;
    for ( IndexArray idxa : idxas ) {
      Index jdat = dft.globalIndex(idxa);
      if ( first ) first = false;
      else cout << ", ";
      cout << "[" << setw(3) << idxa[0] << "," << setw(3) << idxa[1] << "]";
      if ( jdat != idat ) cout << " --> " << jdat << " (!)";
      if ( jdat != idat ) ++nerr;
      ++cntij;
    }
    cout << endl;
    assert( nerr == 0 );
  }
  assert( cntij == nsam0*nsam1 );

  cout << line << endl;
  cout << myname << "Create test data." << endl;
  vector<Complex> datc(nComplexExp);
  vector<Float> datf(2*nComplexExp);
  Index iflt = 0;
  Index isam0 = 0;
  Index isam1 = 0;
  IndexArray isams({isam0, isam1});
  Index idat;
  for ( idat=0; idat<nComplexExp; ++idat ) {
    IndexArrayVector arrs = dft.indexArrays(idat);
    Float re = realval(idat);
    Float im = imagval(idat);
    datf[iflt++] = re;
    datf[iflt++] = im;
    Complex cval(re, im);
    datc[idat] = cval;
  }
  assert( idat == datc.size() );
  assert( iflt == datf.size() );

  cout << line << endl;
  cout << myname << "Copy data in." << endl;
  assert( dft.copyIn(nsams, datc) == 0 );
  for ( Index idat=0; idat<nComplexExp; ++idat ) {
    IndexArrayVector idxas = dft.indexArrays(idat);
    assert( idxas.size() );
    Complex cval = dft.value(idat);
    Index ire = cval.real() + 0.1;
    Index iim = cval.imag() + 0.1;
    cout << myname << setw(3) << idat << ": ["
         << setw(3) << idxas[0][0] << "," << setw(3) << idxas[0][1] << "]: "
         << "(" << setw(3) << ire << "," << setw(3) << iim << ")"
         << endl;
  }

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
  Index stat = 0;
  stat += test_FftwReal2dDftData(3, 4);
  stat += test_FftwReal2dDftData(4, 3);
  return stat;
}

//**********************************************************************
