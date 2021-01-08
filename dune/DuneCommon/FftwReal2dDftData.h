// FftwReal2dDftData.h
//
// David Adams
// January 2021
//
// Concrete class that holds the representation of 2D DFT data used in FFTW.
//
// For real data[i0][i1] with dimensions n0 and n1, the DFT data is stored as 
// complex values dft[i0][i1] with dimensions n0 and (n1/2 + 1).
// The values for i1 > n1/2 may be otained with
//   dft[i0][i1] = CC{dft[n0-i0][n1-i1]}    for i0 > 0
//   dft[0][i1] = CC{dft[0][n1-i1]}
// where CC denotes complex conjugation.

#ifndef FftwReal2dDftData_H
#define FftwReal2dDftData_H

#include "dune/DuneCommon/Real2dDftData.h"
#include "fftw3.h"

//**********************************************************************

template<typename F>
class FftwReal2dDftData : public Real2dDftData<F> {

public:

  using typename Real2dDftData<F>::Index;
  using typename Real2dDftData<F>::IndexArray;
  using typename Real2dDftData<F>::Float;
  using typename Real2dDftData<F>::Complex;     // same memory layout as fftw_complex
  using ComplexVector = std::vector<Complex>;
  using typename Real2dDftData<F>::Norm;

  // FFTW DFT data size for specified sample sizes.
  // Returns the number of floats, i.e. twice the number of complex values.
  template<std::size_t N>
  static Index dftDataSize(const std::array<Index, N>& nsams) {
    Index ndat = 1;
    Index ndim = nsams.size();
    for ( Index idim=0; idim<ndim; ++idim ) {
      Index nsam = nsams[idim];
      Index fac = idim == ndim - 1 ? nsam/2 + 1 : nsam;
      ndat *= fac;
    }
    return 2*ndat;
  }

  // Ctor.
  FftwReal2dDftData(Norm norm, const IndexArray& nsams)
  : m_norm(norm) {
    reset(nsams);
  }

  // Normalization.
  Norm normalization() const override { return m_norm; }

  // Clear data.
  void clear() override {
    for ( Index& nsam : m_nsams ) nsam = 0;
    m_data.clear();
  }

  // Reset the DFT data.
  void reset(const IndexArray& nsams) override {
    clear();
    m_nsams = nsams;
    Index ndat = dftDataSize(nsams);
    m_data.resize(ndat, 0.0);
  }

  // Return the dimensions.
  const IndexArray& nSamples() const override {
    return m_nsams;
  }

  // Access the raw data.
  Complex* data() { return m_data.data(); }

  // Access the raw data by component: real0, imag0, real1, ...
  Float* floatData() { return reinterpret_cast<Float*>(m_data.data()); }
  const Float* floatData() const { return reinterpret_cast<const Float*>(m_data.data()); }

  // Access on entry in the data.
  Complex value(const IndexArray& ifrqs) const override {
    Index ifrq0 = ifrqs[0];
    Index ifrq1 = ifrqs[1];
    Index nfrq0 = m_nsams[0];
    if ( ifrq0 >= nfrq0 ) return this->badValue();
    Index nfrq1 = m_nsams[1];
    if ( ifrq1 >= nfrq1 ) return this->badValue();
    Index n1 = nfrq1/2 + 1;
    if ( ifrq1 >= n1 ) {
      if ( ifrq0 > 0 ) ifrq0 = nfrq0 - ifrq0;
      ifrq1 = nfrq1 - ifrq1;
    }
    Index idat = n1*ifrq0 + ifrq1;
    return m_data[idat];
  }

  // Move data in.
  int moveIn(const IndexArray nsams, ComplexVector& data) {
    m_nsams = nsams;
    m_data = std::move(data);
    return 0;
  }

  // Copy data in.
  int copyIn(const IndexArray nsams, const ComplexVector& data) {
    m_nsams = nsams;
    m_data = data;
    return 0;
  }

  // Move data out.
  int moveOut(ComplexVector& data) {
    data = m_data;
    clear();
    return 0;
  }

  // Copy the data out.
  int copyOut(ComplexVector& data) const {
    data = m_data;
    return 0;
  }

private:

  // Data.
  Norm m_norm;
  IndexArray m_nsams;
  ComplexVector m_data;

};

//**********************************************************************

#endif
