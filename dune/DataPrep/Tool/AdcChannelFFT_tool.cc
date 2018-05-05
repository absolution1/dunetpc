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
  m_NormOpt(ps.get<Index>("NormOpt")) {
  const string myname = "AdcChannelFFT::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "            LogLevel: " << m_LogLevel << endl;
    cout << myname << "           FirstTick: " << m_FirstTick << endl;
    cout << myname << "               NTick: " << m_NTick << endl;
    cout << myname << "             NormOpt: " << m_NormOpt << endl;
  }
}

//**********************************************************************

DataMap AdcChannelFFT::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelFFT::view: ";
  DataMap ret;
  Index isam0 = m_FirstTick;
  if ( isam0 >= acd.samples.size() ) {
    cout << myname << "No data in range." << endl;
    return ret.setStatus(1);
  }
  Index nsam = acd.samples.size() - isam0;
  if ( m_NTick > 0 && m_NTick < nsam ) nsam = m_NTick;
  int nsamInt = nsam;
  vector<double> sams(nsam);
  for ( Index isam=0; isam<nsam; ++isam ) sams[isam] = acd.samples[isam0 + isam];
  TVirtualFFT* pfft = TVirtualFFT::FFT(1, &nsamInt, "R2C");
  pfft->SetPoints(&sams[0]);
  pfft->Transform();
  double xre = 0.0;
  double xim = 0.0;
  Index nxsm = nsam/2 + 1;
  DataMap::FloatVector xres(nsam);
  DataMap::FloatVector xims(nsam);
  DataMap::FloatVector xmgs(nxsm);
  DataMap::FloatVector xphs(nxsm);
  for ( Index ixsm=0; ixsm<nsam; ++ixsm ) {
    pfft->GetPointComplex(ixsm, xre, xim);
    double xmg = sqrt(xre*xre + xim*xim);
    double xph = atan2(xim, xre);
    xres[ixsm] = xre;
    xims[ixsm] = xim;
    if ( ixsm < nxsm ) {
      xmgs[ixsm] = xmg;
      xphs[ixsm] = xph;
    }
    if ( m_LogLevel >= 2 ) cout << myname << ixsm << ": (" << xre << ", " << xim << "): " << xmg << " @ " << xph << endl;
  }
  if ( m_NormOpt == 1 || m_NormOpt == 2 ) {
    float nfac = m_NormOpt == 1 ? 1.0/sqrt(nsam) : 1.0/nsam;
    for ( Index ixsm=0; ixsm<nsam; ++ixsm ) {
      xres[ixsm] *= nfac;
      xims[ixsm] *= nfac;
    }
    bool nsamIsEven = 2*(nsam/2) == nsam;
    for ( Index ixsm=0; ixsm<nxsm; ++ixsm ) {
      float fac = sqrt(2.0)*nfac;
      if ( (ixsm == 0) || (nsamIsEven && (ixsm+1 == nxsm)) ) fac = nfac;
      xmgs[ixsm] *= fac;
    }
  }
  ret.setInt("fftTick0", isam0);
  ret.setInt("fftNTick", nsam);
  ret.setInt("fftNFreq", nxsm);
  ret.setFloatVector("fftReals",  xres);
  ret.setFloatVector("fftImags",  xims);
  ret.setFloatVector("fftMags",   xmgs);
  ret.setFloatVector("fftPhases", xphs);
  return ret;
}

//**********************************************************************
