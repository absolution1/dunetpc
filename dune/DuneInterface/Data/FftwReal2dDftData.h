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

#include "dune/DuneInterface/Data/Real2dDftData.h"
#include "fftw3.h"

//**********************************************************************

template<typename F> class FftwReal2dDftData;
using FftwFloat2dDftData = FftwReal2dDftData<float>;
using FftwDouble2dDftData = FftwReal2dDftData<double>;

//**********************************************************************

template<typename F>
class FftwReal2dDftData : public Real2dDftData<F> {

public:

  using typename Real2dDftData<F>::Index;
  using typename Real2dDftData<F>::IndexArray;
  using IndexArrayVector = std::vector<IndexArray>;
  using typename Real2dDftData<F>::Float;
  using typename Real2dDftData<F>::Complex;     // same memory layout as fftw_complex
  using ComplexVector = std::vector<Complex>;
  using typename Real2dDftData<F>::Norm;

  // FFTW DFT data size for specified sample sizes.
  // Returns the number of complex values.
  template<std::size_t N>
  static Index dftComplexDataSize(const std::array<Index, N>& nsams) {
    Index ndat = 1;
    Index ndim = nsams.size();
    for ( Index idim=0; idim<ndim; ++idim ) {
      Index nsam = nsams[idim];
      Index fac = idim == ndim - 1 ? nsam/2 + 1 : nsam;
      ndat *= fac;
    }
    return ndat;
  }

  // FFTW DFT data size for specified sample sizes.
  // Returns the number of floats, i.e. twice the number of complex values.
  template<std::size_t N>
  static Index dftFloatDataSize(const std::array<Index, N>& nsams) {
    return 2*dftComplexDataSize(nsams);
  }

  // Ctor.
  FftwReal2dDftData(Norm norm, const IndexArray& nsams)
  : m_norm(norm) {
    reset(nsams);
  }

  // Normalization.
  Norm normalization() const override { return m_norm; }

  // The number of dimensions.
  Index rank() const { return 2; }

  // Clear data.
  void clear() override {
    for ( Index& nsam : m_nsams ) nsam = 0;
    m_data.clear();
  }

  // Reset the DFT data.
  void reset(const IndexArray& nsams) override {
    clear();
    m_nsams = nsams;
    Index ndat = dftComplexDataSize(nsams);
    m_data.resize(ndat, 0.0);
  }

  // Return the dimensions.
  const IndexArray& nSamples() const override {
    return m_nsams;
  }
  Index size() const { return m_data.size(); }
  Index nSamples(Index idim) const {
    if ( idim > rank() ) return 0;
    return nSamples()[idim];
  }
  Index size(Index idim) const { return nSamples(idim); }

  // Access the raw data.
  Complex* data() { return m_data.data(); }
  const Complex* data() const { return m_data.data(); }

  // Access the raw data by component: real0, imag0, real1, ...
  Float* floatData() { return reinterpret_cast<Float*>(m_data.data()); }
  const Float* floatData() const { return reinterpret_cast<const Float*>(m_data.data()); }

  // Return the global index for an index array.
  // Returns the number of global indices if any index is out of range.
  Index globalIndex(const IndexArray& ifrqs) const {
    Index ifrq0 = ifrqs[0];
    Index ifrq1 = ifrqs[1];
    Index nfrq0 = m_nsams[0];
    Index nfrq1 = m_nsams[1];
    Index n1 = nfrq1/2 + 1;
    if ( ifrq0 >= nfrq0 ) {
      return size();
    }
    if ( ifrq1 >= nfrq1 ) {
      return size();
    }
    if ( ifrq1 >= n1 ) {
      if ( ifrq0 > 0 ) ifrq0 = nfrq0 - ifrq0;
      ifrq1 = nfrq1 - ifrq1;
    }
    return n1*ifrq0 + ifrq1;
  }

  // Return the index arrays for a global index.
  IndexArrayVector indexArrays(Index idat) const  {
    IndexArrayVector arrs;
    Index nfrq0 = m_nsams[0];
    Index nfrq1 = m_nsams[1];
    Index n1 = nfrq1/2 + 1;
    Index ifrq0 = idat/n1;
    if ( ifrq0 >= nfrq0 ) return arrs;
    Index ifrq1 = idat%n1;
    arrs.push_back(IndexArray({ifrq0, ifrq1}));
    Index jfrq0 = ifrq0 > 0 ? nfrq0 - ifrq0 : 0;
    if ( ifrq1 == 0 ) return arrs;
    Index jfrq1 = nfrq1 - ifrq1;
    if ( jfrq1 == ifrq1 ) return arrs;
    arrs.push_back(IndexArray({jfrq0, jfrq1}));
    return arrs;
  }

  // Return the value for a global index.
  Complex value(Index idat) const {
    if ( idat >= size() ) return this->badValue();
    return m_data[idat];
  }
 
  // Return the value for an index array.
  // The index arrays are checked first iff pchk != nullptr.
  // If any of these checks fail or the calculated index is out of range,
  // badValue() is returned.
  Complex value(const IndexArray& isams) const override {
    Index idat = globalIndex(isams);
    return value(idat);
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
