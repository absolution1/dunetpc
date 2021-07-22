// test_Fw2dFFT.cxx
//
// David Adams
// April 2019
//
// Test Fw2dFFT.

#include "dune/DuneCommon/Utility/Fw2dFFT.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;
using std::vector;
using std::setw;
using std::fixed;
using std::setprecision;

using Index = unsigned int;
using DataFloat = Fw2dFFT::DataFloat;
using DataVector = std::vector<DataFloat>;
using Data = Fw2dFFT::Data;
using DftData = Fw2dFFT::DFT;
using IndexArray = DftData::IndexArray;
using Complex = Fw2dFFT::Complex;

template<typename T>
void printValue(T val) {
  Index prec = 3;
  Index wid = 10;
  cout << setw(wid) << fixed << setprecision(prec) << val;
}

template<>
void printValue(Complex val) {
  printValue(std::real(val));
  cout << ", ";
  printValue(std::imag(val));
}

template<class Data>
void printData(const Data& dat) {
  string myname = "printData: ";
  Index nrow = dat.nSamples()[0];
  Index ncol = dat.nSamples()[1];
  IndexArray ifrqs;
  Index& irow = ifrqs[0];
  Index& icol = ifrqs[1];
  for ( irow=0; irow<nrow; ++irow ) {
    for ( icol=0; icol<ncol; ++icol ) {
       cout << myname;
       cout << setw(4) << irow;
       cout << setw(4) << icol;
       cout << ": ";
       printValue(dat.value(ifrqs));
       cout << endl;
    }
  }
}

template<class Data>
bool printData(const Data& dat1, const Data& dat2) {
  string myname = "printData: ";
  assert( dat1.nSamples() == dat2.nSamples() );
  Index nrow = dat1.nSamples()[0];
  Index ncol = dat1.nSamples()[1];
  IndexArray ifrqs;
  Index& irow = ifrqs[0];
  Index& icol = ifrqs[1];
  int nbad = 0;
  for ( irow=0; irow<nrow; ++irow ) {
    for ( icol=0; icol<ncol; ++icol ) {
       cout << myname;
       cout << setw(4) << irow;
       cout << setw(4) << icol;
       float val1 = dat1.value(ifrqs);
       float val2 = dat2.value(ifrqs);
       cout << ": ";
       printValue(val1);
       cout << "  ?= ";
       printValue(val2);
       if ( fabs(val2 - val1) > 0.001 ) {
         cout << "  !";
         ++nbad;
       }
       cout << endl;
    }
  }
   return nbad == 0;
}

//**********************************************************************

int test_Fw2dFFT(Index ignorm, Index itnorm, int loglev) {
  const string myname = "test_Fw2dFFT: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  RealDftNormalization norm(ignorm, itnorm);
  cout << line << endl;
  cout << myname << " Global norm: " << ignorm << " " << norm.globalName() << endl;
  cout << myname << "   Term norm: " << itnorm << " " << norm.termName() << endl;

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  Index n0 = 3;
  Index n1 = 20;
  Index nsam = n0*n1;
  Index ndft = 2*n0*(n1/2 + 1);
  DataVector sams = {
     2.0,  3.0,   2.0,   1.0,    0.0, -1.0,  -2.0,  -3.0, -2.0, -1.0,
     0.0,  1.0,   2.0,   3.0,    2.0,  1.0,   0.0,  -1.0, -2.0, -3.0,
     3.0,  6.0,  16.1,  28.6,  30.2,  27.7,  16.3,   9.6,  4.2, -1.0,
    -2.3, -4.2,  -9.2, -18.6, -21.9, -29.0, -24.3, -14.2, -5.0, -3.0,
    -2.0, -3.0,  -2.0,  -1.0,    0.0,  1.0,   2.0,   3.0,  2.0,  1.0,
     0.0, -1.0,  -2.0,  -3.0,   -2.0, -1.0,   0.0,   1.0,  2.0,  3.0
  };
  assert( sams.size() == nsam );
  cout << myname << "Sample length: " << nsam << endl;
  float samsum = 0.0;
  for ( float sam : sams ) samsum += sam;
  cout << myname << "Sample mean: " << samsum/nsam << endl;
  cout << myname << "               n0: " << n0 << endl;
  cout << myname << "               n1: " << n1 << endl;
  cout << myname << "        # samples: " << nsam << endl;
  cout << myname << "            # DFT: " << ndft << endl;

  cout << myname << line << endl;
  cout << myname << "Create the sample counts." << endl;
  IndexArray nsams;
  nsams[0] = n0;
  nsams[1] = n1;

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  Data dat(nsams, sams);
  cout << myname << "      Rank: " << dat.rank() << endl;
  cout << myname << "      Size: " << dat.size() << endl;
  assert( dat.isValid() );
  printData(dat);
  assert( dat.size() == nsam );

  cout << myname << line << endl;
  cout << myname << "Create an empty DFT." << endl;
  DftData dft(norm, nsams);
  assert( dft.dataSize(nsams) == nsam );
  assert( dft.dftFloatDataSize(nsams) == ndft );
  printData(dft);

  cout << myname << line << endl;
  cout << myname << "Build transform." << endl;
  Fw2dFFT xf(ndft, 0);

  cout << myname << line << endl;
  cout << myname << "Build plans." << endl;
  assert( ! xf.haveForwardPlan(nsams) );
  assert( ! xf.haveBackwardPlan(nsams) );
  xf.forwardPlan(nsams);
  xf.backwardPlan(nsams);
  assert( xf.haveForwardPlan(nsams) );
  assert( xf.haveBackwardPlan(nsams) );

  cout << myname << line << endl;
  cout << myname << "Transform forward." << endl;
  int fstat = xf.fftForward(dat, dft, loglev);
  assert( fstat == 0 );
  printData(dft);

  cout << myname << line << endl;
  cout << myname << "Transform backward." << endl;
  Data dat2;
  int bstat = xf.fftBackward(dft, dat2, loglev);
  assert( bstat == 0 );
  assert( printData(dat, dat2) );

/*
  cout << myname << line << endl;
  cout << myname << "Check power." << endl;
  float pwr1 = 0.0;
  for ( float sam : sams ) pwr1 += sam*sam;
  float pwr2 = dft.power();
  float pwr3 = 0.0;
  float pwr4 = 0.0;
  for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
    float amp = dft.amplitude(ifrq);
    float xre = dft.real(ifrq);
    float xim = dft.imag(ifrq);
    if ( !dft.normalization().isPower() || ifrq < dft.nCompact() ) {
      pwr3 += amp*amp;
      pwr4 += xre*xre + xim*xim;
    }
  }
  double pfac = 1.0;
  if ( dft.normalization().isStandard() ) pfac = 1.0/nsam;
  if ( dft.normalization().isBin() ) pfac = nsam;
  pwr3 *= pfac;
  pwr4 *= pfac;
  cout << myname << "Tick power: " << pwr1 << endl;
  cout << myname << " DFT power: " << pwr2 << endl;
  cout << myname << " Amp power: " << pwr3 << endl;
  cout << myname << " ReI power: " << pwr4 << endl;
  assert( fabs(pwr2 - pwr1) < 1.e-5*pwr1 );

  cout << myname << line << endl;
  cout << myname << "Call inverse." << endl;
  FloatVector sams2, xres2, xims2;
  int rstat = xf.fftInverse(dft, sams2, loglev);
  if ( rstat != 0 ) {
    cout << myname << "FFT invert returned " << rstat << endl;
    assert( false );
  }

  assert( sams2.size() == nsam );
  for ( Index isam=0; isam<nsam; ++isam ) {
    cout << setw(4) << isam << ":" << setw(10) << fixed << sams[isam]
         << " ?= " << setw(10) << fixed << sams2[isam] << endl;
    assert( fabs(sams2[isam] - sams[isam]) < 1.e-4 );
  }

*/
  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  Index ignorm = 1;
  Index itnorm = 1;
  int loglev = 0;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ignorm [itnorm [loglev]]]" << endl;
      return 0;
    }
    ignorm = std::stoi(sarg);
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    itnorm = std::stoi(sarg);
  }
  if ( argc > 3 ) {
    string sarg(argv[3]);
    loglev = std::stoi(sarg);
  }
  return test_Fw2dFFT(ignorm, itnorm, loglev);
}

//**********************************************************************
