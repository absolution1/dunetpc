// DuneFFT.cxx
#include "DuneFFT.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "TVirtualFFT.h"
#include "TComplex.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::setw;
using std::fixed;

//**********************************************************************
// Class methods.
//**********************************************************************

int DuneFFT::
fftForward(Index nsam, const float* psam, DFT& dft, Index logLevel) {
  const string myname = "DuneFFT::fftForward: ";
  dft.reset(nsam);
  if ( ! dft.isValid() ) return 1;
  if ( nsam == 0 ) return 0;
  vector<double> sams(nsam);
  for ( Index isam=0; isam<nsam; ++isam ) sams[isam] = psam[isam];
  int nsamInt = nsam;
  TVirtualFFT* pfft = TVirtualFFT::FFT(1, &nsamInt, "R2C");
  pfft->SetPoints(&sams[0]);
  pfft->Transform();
  double xre = 0.0;
  double xim = 0.0;
  Index namp = dft.nAmplitude();
  Index npha = dft.nPhase();
  FloatVector amps(namp);
  FloatVector phas(npha);
  // Loop over the compact samples.
  float nfac = 1.0;
  if ( dft.normalization().isConsistent() ) nfac = 1.0/sqrt(nsam);
  if ( dft.normalization().isBin() ) nfac = 1.0/nsam;
  for ( Index ifrq=0; ifrq<namp; ++ifrq ) {
    pfft->GetPointComplex(ifrq, xre, xim);
    // For an even # samples (namp = npha + 1), the Nyquist term is real
    // and we store the sign with the amplitude.
    double xam = ifrq < npha ? sqrt(xre*xre + xim*xim) : xre;
    double xph = atan2(xim, xre);
    float fac = nfac;
    if ( dft.normalization().isPower() && dft.isAliased(ifrq) ) fac *= sqrt(2.0);
    xam *= fac;
    dft.setAmplitude(ifrq, xam);
    if ( ifrq < npha ) {
      dft.setPhase(ifrq, xph);
    }
    if ( logLevel >= 3 ) {
      cout << myname << setw(4) << ifrq << ": ("
           << setw(10) << fixed << xre << ", "
           << setw(10) << fixed << xim << "): "
           << setw(10) << fixed << xam;
      if ( ifrq < npha ) cout << " @ " << setw(10) << fixed << xph;
      cout << endl;
    }
  }
  return 0;
}

//**********************************************************************

int DuneFFT::
fftForward(const FloatVector& sams, DFT& dft, Index logLevel) {
  return fftForward(sams.size(), &sams[0], dft, logLevel);
}

//**********************************************************************

int DuneFFT::
fftInverse(const DFT& dft, FloatVector& sams, Index logLevel) {
  const string myname = "DuneFFT::fftInverse: ";
  if ( ! dft.isValid() ) return 1;
  Index namp = dft.nAmplitude();
  Index npha = dft.nPhase();
  if ( namp == 0 || npha == 0 ) return 1;
  if ( namp < npha ) return 2;
  if ( namp - npha > 1 ) return 3;
  Index nsam = namp + npha - 1;
  int nsamInt = nsam;
  TVirtualFFT* pfft = TVirtualFFT::FFT(1, &nsamInt, "C2R");
  vector<double> xdres(nsam, 0.0);
  vector<double> xdims(nsam, 0.0);
  for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
    double amp = dft.convAmplitude(ifrq);
    double pha = dft.phase(ifrq);
    double xre = amp*cos(pha);
    double xim = amp*sin(pha);
    xdres[ifrq] = xre;
    xdims[ifrq] = xim;
  }
  pfft->SetPointsComplex(&xdres[0], &xdims[0]);
  if ( logLevel >= 3 ) {
    cout << myname << "Inverting" << endl;
    cout << myname << "Real/imag components:" << endl;
    for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
      float xre = xdres[ifrq];
      float xim = xdims[ifrq];
      cout << myname << setw(4) << ifrq << ": (" << setw(10) << fixed << xre
           << ", " << setw(10) << fixed << xim << ")" << endl;
    }
  }
  pfft->Transform();
  sams.resize(nsam);
  float nfac = 1.0/nsam;
  for ( Index isam=0; isam<nsam; ++isam ) sams[isam] = nfac*pfft->GetPointReal(isam);
  return 0;
}

//**********************************************************************
