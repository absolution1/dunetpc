// FwFFT.cxx
#include "FwFFT.h"
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

FwFFT::FwFFT(Index nsamMax, Index opt)
: m_nsamMax(nsamMax),
  m_flag(opt==2 ? FFTW_PATIENT : opt==1 ? FFTW_MEASURE : FFTW_ESTIMATE),
  m_inData((Float*) fftw_malloc(m_nsamMax*sizeof(Float))),
  m_outData((Complex*) fftw_malloc(m_nsamMax*sizeof(Complex))) { }

//**********************************************************************

FwFFT::~FwFFT() {
  for ( auto& iplan : m_forwardPlans ) fftw_destroy_plan(iplan.second);
  for ( auto& iplan : m_backwardPlans ) fftw_destroy_plan(iplan.second);
  fftw_free(m_inData);
  fftw_free(m_outData);
}

//**********************************************************************

FwFFT::Plan& FwFFT::forwardPlan(Index nsam) {
  const string myname = "FwFFT::forwardPlan";
  static Plan badplan;
  if ( nsam > m_nsamMax ) {
    cout << myname << "Sample count is too large. Maximum is " << m_nsamMax << endl;
    return badplan;
  }
  if ( m_forwardPlans.count(nsam) == 0 ) {
    m_forwardPlans[nsam] = fftw_plan_dft_r2c_1d(nsam, m_inData, m_outData, m_flag);
  }
  return m_forwardPlans[nsam];
}

//**********************************************************************

FwFFT::Plan& FwFFT::backwardPlan(Index nsam) {
  const string myname = "FwFFT::backwardPlan";
  static Plan badplan;
  if ( nsam > m_nsamMax ) {
    cout << myname << "Sample count is too large. Maximum is " << m_nsamMax << endl;
    return badplan;
  }
  if ( m_backwardPlans.count(nsam) == 0 ) {
    m_backwardPlans[nsam] = fftw_plan_dft_c2r_1d(nsam, m_outData, m_inData, m_flag);
  }
  return m_backwardPlans[nsam];
}

//**********************************************************************

int FwFFT::
fftForward(Index nsam, const float* psam, DFT& dft, Index logLevel) {
  const string myname = "FwFFT::fftForward: ";
  if ( nsam > m_nsamMax ) {
    cout << myname << "Sample count is too large. Maximum is " << m_nsamMax << endl;
    return 2;
  }
  dft.reset(nsam);
  if ( ! dft.isValid() ) return 1;
  if ( nsam == 0 ) return 0;
  for ( Index isam=0; isam<nsam; ++isam ) m_inData[isam] = psam[isam];
  Plan& plan = forwardPlan(nsam);
  fftw_execute(plan);
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
    xre = m_outData[ifrq][0];
    xim = m_outData[ifrq][1];
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

int FwFFT::
fftForward(const FloatVector& sams, DFT& dft, Index logLevel) {
  return fftForward(sams.size(), &sams[0], dft, logLevel);
}

//**********************************************************************

int FwFFT::
fftInverse(const DFT& dft, FloatVector& sams, Index logLevel) {
  const string myname = "FwFFT::fftInverse: ";
  if ( ! dft.isValid() ) return 1;
  Index namp = dft.nAmplitude();
  Index npha = dft.nPhase();
  if ( namp == 0 || npha == 0 ) return 1;
  if ( namp < npha ) return 2;
  if ( namp - npha > 1 ) return 3;
  Index nsam = namp + npha - 1;
  if ( nsam > m_nsamMax ) {
    cout << myname << "Sample count is too large. Maximum is " << m_nsamMax << endl;
    return 4;
  }
  for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
    double amp = dft.convAmplitude(ifrq);
    double pha = dft.phase(ifrq);
    double xre = amp*cos(pha);
    double xim = amp*sin(pha);
    m_outData[ifrq][0] = xre;
    m_outData[ifrq][1] = xim;
  }
  Plan& plan = backwardPlan(nsam);
  fftw_execute(plan);
  float nfac = 1.0/nsam;
  sams.resize(nsam);
  for ( Index isam=0; isam<nsam; ++isam ) sams[isam] = nfac*m_inData[isam];
  return 0;
}

//**********************************************************************
