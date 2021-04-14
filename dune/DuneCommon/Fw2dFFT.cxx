// Fw2dFFT.cxx
#include "Fw2dFFT.h"
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

Fw2dFFT::Fw2dFFT(Index ndatMax, Index opt)
: m_ndatMax(ndatMax),
  m_flag(opt==2 ? FFTW_PATIENT : opt==1 ? FFTW_MEASURE : FFTW_ESTIMATE),
  m_inData(reinterpret_cast<DftFloat*>(fftw_malloc(m_ndatMax*sizeof(DftFloat)))),
  m_outData(reinterpret_cast<Complex*>(fftw_malloc(m_ndatMax*sizeof(DftFloat)))) { }

//**********************************************************************

Fw2dFFT::~Fw2dFFT() {
  for ( auto& iplan : m_forwardPlans ) fftw_destroy_plan(iplan.second);
  for ( auto& iplan : m_backwardPlans ) fftw_destroy_plan(iplan.second);
  fftw_free(m_inData);
  fftw_free(m_outData);
}

//**********************************************************************

bool Fw2dFFT::checkDataSize(const IndexArray& nsams) const {
  Index ndatIn = Real2dDftData<DataFloat>::dataSize(nsams);
  Index ndatOut = FftwReal2dDftData<DftFloat>::dftFloatDataSize(nsams);
  if ( ndatIn > m_ndatMax ) return 1;
  if ( ndatOut > m_ndatMax ) return 2;
  if ( ndatIn == 0 ) return 3;
  return 0;
}

//**********************************************************************

bool Fw2dFFT::haveForwardPlan(const IndexArray& nsams) const {
  return m_forwardPlans.count(nsams);
}

//**********************************************************************

bool Fw2dFFT::haveBackwardPlan(const IndexArray& nsams) const {
  return m_backwardPlans.count(nsams);
}

//**********************************************************************

Fw2dFFT::Plan& Fw2dFFT::forwardPlan(const IndexArray& nsams) {
  const string myname = "Fw2dFFT::forwardPlan: ";
  static Plan badplan;
  if ( checkDataSize(nsams) ) {
    cout << myname << "Data cannot be accomodated. Maximum data size is " << m_ndatMax << endl;
    return badplan;
  }
  if ( m_forwardPlans.count(nsams) == 0 ) {
    m_forwardPlans[nsams] = fftw_plan_dft_r2c_2d(nsams[0], nsams[1], m_inData, fftwOutData(), m_flag);
  }
  return m_forwardPlans[nsams];
}

//**********************************************************************

Fw2dFFT::Plan& Fw2dFFT::backwardPlan(const IndexArray& nsams) {
  const string myname = "Fw2dFFT::backwardPlan: ";
  static Plan badplan;
  if ( checkDataSize(nsams) ) {
    cout << myname << "Data cannot be accomodated. Maximum data size is " << m_ndatMax << endl;
    return badplan;
  }
  if ( m_backwardPlans.count(nsams) == 0 ) {
    m_backwardPlans[nsams] = fftw_plan_dft_c2r_2d(nsams[0], nsams[1], fftwOutData(), m_inData, m_flag);
  }
  return m_backwardPlans[nsams];
}

//**********************************************************************

int Fw2dFFT::
fftForward(const Data& dat, DFT& dft, Index logLevel) {
  const string myname = "Fw2dFFT::fftForward: ";
  if ( ! dat.isValid() ) return 1;
  auto nsams = dat.nSamples();
  dft.reset(nsams);
  if ( ! dft.isValid() ) return 2;
  if ( dft.normalization().isPower() ) {
    cout << "ERROR: Power normalization is not (yet) supported." << endl;
    return 3;
  }
  Index ndatIn = dat.size();
  Index ndatOut = DFT::dftFloatDataSize(nsams);
  Plan& plan = forwardPlan(nsams);
  float fndat = ndatIn;
  float nfac = dft.normalization().isStandard()   ? 1.0             :
               dft.normalization().isConsistent() ? 1.0/sqrt(fndat) :
               dft.normalization().isBin()        ? 1.0/fndat       : 0.0;
  if ( dat.copyDataOut(m_inData) ) {
    cout << myname << "ERROR: Copy of input data failed." << endl;
    return 4;
  }
  fftw_execute(plan);
  dft.reset(nsams);
  for ( Index idat=0; idat<ndatOut; ++idat ) dft.floatData()[idat] = nfac*floatOutData()[idat];
  return 0;
}

//**********************************************************************

int Fw2dFFT::
fftBackward(const DFT& dft, Data& dat, Index logLevel) {
  const string myname = "Fw2dFFT::fftBackward: ";
  if ( ! dft.isValid() ) return 1;
  auto nsams = dft.nSamples();
  if ( checkDataSize(nsams) ) {
    cout << myname << "Sample counts are too large. Maximum data size is " << m_ndatMax << endl;
    return 2;
  }
  Index ndatIn = Real2dDftData<DataFloat>::dataSize(nsams);
  Index ndatOut = FftwReal2dDftData<DftFloat>::dftFloatDataSize(nsams);
  if ( dft.normalization().isPower() ) {
    cout << myname << "ERROR: Power normalization is not (yet) supported." << endl;
    return 3;
  }
  dat.reset(nsams);
  if ( ! dat.isValid() ) {
    cout << myname << "ERROR: Unable to initialize output data container." << endl;
    return 4;
  }
  float fndat = ndatIn;
  float nfac = dft.normalization().isStandard()   ? 1.0/fndat       :
               dft.normalization().isConsistent() ? 1.0/sqrt(fndat) :
               dft.normalization().isBin()        ? 1.0             : 0.0;
  Plan& plan = backwardPlan(nsams);
  for ( Index idat=0; idat<ndatOut; ++idat ) floatOutData()[idat] = nfac*dft.floatData()[idat];
  fftw_execute(plan);
  dat.copyDataIn(m_inData);
  return 0;
}

//**********************************************************************
