// AdcChannelFFT_tool.cc

#include "AdcChannelFFT.h"
#include <iostream>
#include <sstream>
#include <vector>
#include "TVirtualFFT.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;

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
    cout << myname << "           ReturnOpt: " << m_ReturnOpt << endl;
  }
}

//**********************************************************************

DataMap AdcChannelFFT::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelFFT::view: ";
  DataMap ret;
  bool doForward = false;
  bool doInverse = false;
  bool returnForward = false;
  bool returnInverse = false;
  if ( m_Action <= 4 ) {
    returnForward = true;
    if ( m_Action > 0 ) {
      if ( m_Action == 2 || m_Action == 4 ) {
        doForward = acd.dftmags.size() == 0;
      } else {
        doForward = true;
      }
    }
  } else if ( m_Action >= 10  && m_Action <= 14 ) {
    returnInverse = true;
    if ( m_Action != 10 ) {
      if ( m_Action == 2 || m_Action == 4 ) {
        doForward = acd.dftmags.size() == 0;
      } else {
        doForward = true;
      }
    }
  } else {
    cout << myname << "ERROR: Invalid action: " << m_Action << endl;
    return ret.setStatus(1);
  }
  DataMap::FloatVector sams;
  DataMap::FloatVector xres;
  DataMap::FloatVector xims;
  DataMap::FloatVector xmgs;
  DataMap::FloatVector xphs;
  Index isam0 = 0;
  Index nsam = 0;
  if ( doForward ) {
    isam0 = m_FirstTick;
    if ( isam0 >= acd.samples.size() ) {
      cout << myname << "No data in range." << endl;
      return ret.setStatus(11);
    }
    nsam = acd.samples.size() - isam0;
    if ( m_NTick > 0 && m_NTick < nsam ) nsam = m_NTick;
    int rstat = fftForward(m_NormOpt, nsam, &acd.samples[isam0], xres, xims, xmgs, xphs);
    if ( rstat ) return ret.setStatus(ret);
  } else if ( doInverse ) {
    cout << myname << "Inverse transform is not yet supported." << endl;
    return ret.setStatus(2);
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
  return ret;
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
  for ( Index ixsm=0; ixsm<nsam; ++ixsm ) {
    pfft->GetPointComplex(ixsm, xre, xim);
    double xmg = sqrt(xre*xre + xim*xim);
    double xph = atan2(xim, xre);
    xres[ixsm] = xre;
    xims[ixsm] = xim;
    if ( ixsm < nmag ) mags[ixsm] = xmg;
    if ( ixsm < npha ) phases[ixsm] = xph;
    if ( m_LogLevel >= 2 ) cout << myname << ixsm << ": (" << xre << ", " << xim << "): " << xmg << " @ " << xph << endl;
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
fftInverse(Index normOpt, const FloatVector& mags, const FloatVector& phases, FloatVector& sams) const {
  return 1;
}

//**********************************************************************
