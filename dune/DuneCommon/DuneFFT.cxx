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
fftForward(Index normOpt, Index nsam, const float* psam,
           FloatVector& xres, FloatVector& xims,  FloatVector& mags, FloatVector& phases,
           Index logLevel) {
  const string myname = "DuneFFT::fftForward: ";
  Index normOptGlobal = getNormOptGlobal(normOpt);
  bool usePowerNormalization = isPower(normOpt);
  vector<double> sams(nsam);
  for ( Index isam=0; isam<nsam; ++isam ) sams[isam] = psam[isam];
  int nsamInt = nsam;
  TVirtualFFT* pfft = TVirtualFFT::FFT(1, &nsamInt, "R2C");
  pfft->SetPoints(&sams[0]);
  pfft->Transform();
  double xre = 0.0;
  double xim = 0.0;
  Index nmag = (nsam+2)/2;
  Index npha = (nsam+1)/2;
  xres.resize(nsam);
  xims.resize(nsam);
  mags.resize(nmag);
  phases.resize(npha);
  // Loop over the compact samples.
  for ( Index iptc=0; iptc<nmag; ++iptc ) {
    pfft->GetPointComplex(iptc, xre, xim);
    // For an even # samples (nmag = npha + 1), the Nyquist term is real
    // and we store the sign with the magnitude.
    double xmg = iptc < npha ? sqrt(xre*xre + xim*xim) : xre;
    double xph = atan2(xim, xre);
    xres[iptc] = xre;
    xims[iptc] = xim;
    mags[iptc] = xmg;
    if ( iptc < npha ) phases[iptc] = xph;
    if ( logLevel >= 3 ) {
      cout << myname << setw(4) << iptc << ": ("
           << setw(10) << fixed << xre << ", "
           << setw(10) << fixed << xim << "): "
           << setw(10) << fixed << xmg;
      if ( iptc < npha ) cout << " @ " << setw(10) << fixed << xph;
      cout << endl;
    }
  }
  // Loop over the remaining complex samples.
  for ( Index iptc=nmag; iptc<nsam; ++iptc ) {
    xres[iptc] =  xres[nsam-iptc];
    xims[iptc] = -xims[nsam-iptc];
  }
  float nfac = 1.0;
  if ( normOptGlobal == 1 ) nfac = 1.0/sqrt(nsam);
  if ( normOptGlobal == 2 ) nfac = 1.0/nsam;
  for ( Index ixsm=0; ixsm<nsam; ++ixsm ) {
    xres[ixsm] *= nfac;
    xims[ixsm] *= nfac;
  }
  for ( Index ixsm=0; ixsm<nmag; ++ixsm ) {
    float fac = nfac;
    if ( usePowerNormalization && ixsm > 0 && ixsm < npha ) fac *= sqrt(2.0);
    mags[ixsm] *= fac;
  }
  return 0;
}

//**********************************************************************

int DuneFFT::
fftForward(Index normOpt, const FloatVector& sams,
           FloatVector& xres, FloatVector& xims,  FloatVector& mags, FloatVector& phases,
           Index logLevel) {
  return fftForward(normOpt, sams.size(), &sams[0], xres, xims, mags, phases, logLevel);
}

//**********************************************************************

int DuneFFT::
fftInverse(Index normOpt, const FloatVector& mags, const FloatVector& phases,
           FloatVector& xres, FloatVector& xims, FloatVector& sams,
           Index logLevel) {
  const string myname = "DuneFFT::fftInverse: ";
  Index normOptGlobal = getNormOptGlobal(normOpt);
  bool usePowerNormalization = isPower(normOpt);
  Index nmag = mags.size();
  Index npha = phases.size();
  if ( nmag == 0 || npha == 0 ) return 1;
  if ( nmag < npha ) return 2;
  if ( nmag - npha > 1 ) return 3;
  Index nsam = nmag + npha - 1;
  int nsamInt = nsam;
  TVirtualFFT* pfft = TVirtualFFT::FFT(1, &nsamInt, "C2R");
  xres.clear();
  xims.clear();
  xres.resize(nsam, 0.0);
  xims.resize(nsam, 0.0);
  for ( Index imag=0; imag<nmag; ++imag ) {
    double mag = mags[imag];
    if ( usePowerNormalization && imag > 0 && imag != npha ) mag /= sqrt(2.0);
    double pha = imag<npha ? phases[imag] : 0.0;
    double xre = mag*cos(pha);
    double xim = mag*sin(pha);
    Index ifrq = imag;
    xres[ifrq] = xre;
    xims[ifrq] = xim;
    if ( ifrq > 0 ) {
      Index ifrq2 = nsam  - ifrq;
      if ( ifrq2 > ifrq ) {
        xres[ifrq2] = xre;
        xims[ifrq2] = -xim;
      }
    }
  }
  vector<double> xdres(nsam, 0.0);
  vector<double> xdims(nsam, 0.0);
  for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
    xdres[ifrq] = xres[ifrq];
    xdims[ifrq] = xims[ifrq];
  }
  pfft->SetPointsComplex(&xdres[0], &xdims[0]);
  if ( logLevel >= 3 ) {
    cout << myname << "Inverting" << endl;
    cout << myname << "Frequency components:" << endl;
    for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
      //double xre, xim;
      //pfft->GetPointComplex(ifrq, xre, xim, true);
      double xre = xres[ifrq];
      double xim = xims[ifrq];
      cout << myname << setw(4) << ifrq << ": (" << setw(10) << fixed << xre
           << ", " << setw(10) << fixed << xim << ")" << endl;
    }
  }
  pfft->Transform();
  sams.resize(nsam);
  float nfac = 1.0;
  if ( normOptGlobal == 1 ) nfac = 1.0/sqrt(nsam);
  if ( normOptGlobal == 0 ) nfac = 1.0/nsam;
  for ( Index isam=0; isam<nsam; ++isam ) sams[isam] = nfac*pfft->GetPointReal(isam);
  return 0;
}

//**********************************************************************
