// AdcChannelFFT_tool.cc

#include "AdcChannelFFT.h"
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

AdcChannelFFT::AdcChannelFFT(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_FirstTick(ps.get<Index>("FirstTick")),
  m_NTick(ps.get<Index>("NTick")),
  m_NormOpt(ps.get<Index>("NormOpt")),
  m_Action(ps.get<Index>("Action")),
  m_ReturnOpt(ps.get<Index>("ReturnOpt")) {
  const string myname = "AdcChannelFFT::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "            LogLevel: " << m_LogLevel << endl;
    cout << myname << "           FirstTick: " << m_FirstTick << endl;
    cout << myname << "               NTick: " << m_NTick << endl;
    cout << myname << "             NormOpt: " << m_NormOpt << endl;
    cout << myname << "              Action: " << m_Action << endl;
    cout << myname << "           ReturnOpt: " << m_ReturnOpt << endl;
  }
}

//**********************************************************************

DataMap AdcChannelFFT::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelFFT::view: ";
  DataMap ret;
  DataMap::FloatVector sams;
  DataMap::FloatVector mags;
  DataMap::FloatVector phas;
  internalView(acd, sams, mags, phas, ret);
  return ret;
}

//**********************************************************************

DataMap AdcChannelFFT::update(AdcChannelData& acd) const {
  const string myname = "AdcChannelFFT::update: ";
  DataMap ret;
  DataMap::FloatVector sams;
  DataMap::FloatVector mags;
  DataMap::FloatVector phas;
  internalView(acd, sams, mags, phas, ret);
  if ( ret ) return ret;
  if ( m_Action == 3 || m_Action == 4 ) {
    if ( mags.size() ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Saving DFT." << endl;
      acd.dftmags = mags;
      acd.dftphases = phas;
    }
  }
  if ( m_Action == 13 || m_Action == 14 ) {
    if ( sams.size() ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Saving samples." << endl;
      acd.samples = sams;
    }
  }
  return ret;
}

//**********************************************************************

void AdcChannelFFT::
internalView(const AdcChannelData& acd, FloatVector& sams, FloatVector& xmgs, FloatVector& xphs, DataMap& ret) const {
  const string myname = "AdcChannelFFT::internalView: ";
  bool doForward = false;
  bool doInverse = false;
  if ( m_Action <= 4 ) {
    if ( m_Action > 0 ) {
      if ( m_Action == 2 || m_Action == 4 ) {
        doForward = acd.dftmags.size() == 0;
      } else {
        doForward = true;
      }
    }
  } else if ( m_Action >= 10  && m_Action <= 14 ) {
    if ( m_Action != 10 ) {
      if ( m_Action == 12 || m_Action == 14 ) {
        doInverse = acd.samples.size() == 0;
      } else {
        doInverse = true;
      }
    }
  } else {
    cout << myname << "ERROR: Invalid action: " << m_Action << endl;
    ret.setStatus(1);
    return;
  }
  DataMap::FloatVector xres;
  DataMap::FloatVector xims;
  Index isam0 = 0;
  Index nsam = 0;
  if ( doForward ) {
    isam0 = m_FirstTick;
    if ( isam0 >= acd.samples.size() ) {
      cout << myname << "No data in range." << endl;
      ret.setStatus(11);
      return;
    }
    nsam = acd.samples.size() - isam0;
    if ( m_NTick > 0 && m_NTick < nsam ) nsam = m_NTick;
    int rstat = fftForward(m_NormOpt, nsam, &acd.samples[isam0], xres, xims, xmgs, xphs);
    if ( rstat ) {
      ret.setStatus(10+rstat);
      return;
    }
  } else if ( doInverse ) {
    int rstat = fftInverse(m_NormOpt, acd.dftmags, acd.dftphases, xres, xims, sams);
    xmgs = acd.dftmags;
    xphs = acd.dftphases;
    if ( rstat ) {
      ret.setStatus(20+rstat);
      return;
    }
  }
  if ( m_ReturnOpt >= 1 ) {
    ret.setInt("fftTick0", isam0);
    ret.setInt("fftNTick", nsam);
  }
  Index dftRet = m_ReturnOpt % 10;
  if ( dftRet >= 1 ) {
    ret.setInt("fftNMag",   xmgs.size());
    ret.setInt("fftNPhase", xphs.size());
  }
  if ( dftRet >= 2 ) {
    ret.setFloatVector("fftMags",   xmgs);
    ret.setFloatVector("fftPhases", xphs);
  }
  if ( dftRet >= 3 ) {
    ret.setFloatVector("fftReals",  xres);
    ret.setFloatVector("fftImags",  xims);
  }
  if ( m_ReturnOpt >= 10 ) {
    if ( sams.size() ) ret.setFloatVector("fftSamples", sams);
    else ret.setFloatVector("fftSamples", acd.samples);
  }
  return;
}

//**********************************************************************

int AdcChannelFFT::
fftForward(Index normOpt, Index nsam, const float* psam,
           FloatVector& xres, FloatVector& xims,  FloatVector& mags, FloatVector& phases) const {
  const string myname = "AdcChannelFFT::fftForward: ";
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
  // Loop over the complex samples.
  vector<double> xres2(nsam, 9.99);
  vector<double> xims2(nsam, 9.99);
  pfft->GetPointsComplex(&xres2[0], &xims2[0]);
  for ( Index iptc=0; iptc<nsam; ++iptc ) {
    pfft->GetPointComplex(iptc, xre, xim);
    double xmg = sqrt(xre*xre + xim*xim);
    double xph = atan2(xim, xre);
    xres[iptc] = xre;
    xims[iptc] = xim;
    if ( iptc < npha ) {
      mags[iptc] = xmg;
      phases[iptc] = xph;
    // For an even # samples (nmag = npha + 1), the Nyquist term is real.
    } else if ( iptc < nmag ) {
      mags[iptc] = xre;
    }
    if ( m_LogLevel >= 2 ) {
      cout << myname << setw(4) << iptc << ": ("
           << setw(10) << fixed << xre << ", "
           << setw(10) << fixed << xim << "): "
           << setw(10) << fixed << xmg << " @ "
           << setw(10) << fixed << xph << endl;
      cout << myname << "      ("
           << setw(10) << fixed << xres2[iptc] << ", "
           << setw(10) << fixed << xims2[iptc] << ")" << endl;
    }
  }
  float nfac = 1.0;
  if ( m_NormOpt == 1 ) nfac = 1.0/sqrt(nsam);
  if ( m_NormOpt == 2 ) nfac = 1.0/nsam;
  for ( Index ixsm=0; ixsm<nsam; ++ixsm ) {
    xres[ixsm] *= nfac;
    xims[ixsm] *= nfac;
  }
  for ( Index ixsm=0; ixsm<nmag; ++ixsm ) {
    float fac = nfac;
    if ( ixsm > 0 && ixsm < npha ) fac *= sqrt(2.0);
    mags[ixsm] *= fac;
  }
  return 0;
}

//**********************************************************************

int AdcChannelFFT::
fftInverse(Index normOpt, const FloatVector& mags, const FloatVector& phases,
           FloatVector& xres, FloatVector& xims, FloatVector& sams) const {
  const string myname = "AdcChannelFFT::fftInverse: ";
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
    if ( imag > 0 && imag != npha ) mag /= sqrt(2.0);
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
  if ( m_LogLevel >= 3 ) {
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
  if ( m_NormOpt == 1 ) nfac = 1.0/sqrt(nsam);
  if ( m_NormOpt == 0 ) nfac = 1.0/nsam;
  for ( Index isam=0; isam<nsam; ++isam ) sams[isam] = nfac*pfft->GetPointReal(isam);
  return 0;
}

//**********************************************************************
