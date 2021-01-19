// Real2dData.h
//
// David Adams
// January 2021
//
// This class provides a 2D array of real numbers.
//
// It is templated so DFT elements can be stored at different levels of
// floating point precision.

#ifndef Real2dData_H
#define Real2dData_H

#include "dune/DuneCommon/RealDftNormalization.h"
#include <complex>
#include <array>

template<typename F>
class Real2dData {

public:

  using Float = F;
  using DataVector = std::vector<F>;
  using Index = unsigned int;
  using IndexArray = std::array<Index,2>;
  using Norm = RealDftNormalization;
  using Complex = std::complex<Float>;     // same memory layout as fftw_complex

  // Data size for specified sample sizes.
  // Returns the number of floats.
  template<std::size_t N>
  static Index dataSize(const std::array<Index, N>& nsams) {
    Index ndat = 1;
    Index ndim = nsams.size();
    for ( Index idim=0; idim<ndim; ++idim ) {
      Index nsam = nsams[idim];
      ndat *= nsam;
    }
    return ndat;
  }

  // Constructor.
  Real2dData() = default;

  // Constructor from dimension sizes.
  // The data is set to zero.
  Real2dData(const IndexArray& nsams) : m_nsams(nsams), m_data(dataSize(nsams), 0.0) { }

  // Constructor from dimension sizes and data.
  Real2dData(const IndexArray& nsams, const DataVector& data)
  : Real2dData(nsams) {
    copyDataIn(data);
  }

  // Clear the data, i.e. zero the # samples.
  void clear() {
    for ( Index& nsam : m_nsams ) nsam = 0;
    m_data.clear();
  }

  // Reset to (nsamx, nsamy) samples and zero the DFT.
  void reset(const IndexArray& nsams)  {
    m_nsams = nsams;
    m_data.clear();
    Index ndat = dataSize(m_nsams);
    m_data.resize(ndat, 0.0);
  }

  // The number of dimensions.
  Index rank() const { return 2; }

  // Number of samples in each dimension.
  // nSample = size = # samples
  const IndexArray& nSamples() const { return m_nsams; }
  Index nSamples(Index idim) const {
    if ( idim > rank() ) return 0;
    return nSamples()[idim];
  }
  Index size(Index idim) const { return nSamples(idim); }

  // Size of the full data.
  Index size() const { return m_data.size(); }

  // Return the status for a data index.
  bool inRange(Index idim, Index idat) const { return idim < rank() && idat < size(idim); }
  bool isZero(Index idat) const { return idat == 0; }

  // Check if this object is valid.
  // Object must have an assigned size so that the frequency status methods work.
  bool isValid() const {
    if ( rank() != 2 ) return false;
    if ( m_nsams.size() != rank() ) return false;
    if ( size() == 0 ) return false;
    return true;
  }

  // Bad value.
  // Test real and imaginary pars with isnan or self != self;
  Float badValue() const {
    double bad = std::numeric_limits<Float>::quiet_NaN();
    return bad;
  }

  // Return the data (row major order).
  const std::vector<F>& data() const { return m_data; }

  // Copy the data in from a vector.
  // Return 0 for success.
  int copyDataIn(const DataVector& data) {
    if ( ! isValid() ) return 1;
    if ( data.size() != m_data.size() ) return 2;
    m_data = data;
    return 0;
  }

  // Copy the data in from an array of any type.
  // Return 0 for success.
  template<typename T>
  int copyDataIn(const T* pdat) {
    if ( ! isValid() ) return 1;
    Index ndat = size();
    for ( Index idat=0; idat<ndat; ++idat ) {
      m_data[idat] = pdat[idat];
    }
    return 0;
  }

  // Copy the data out to an array of any type.
  // Return 0 for success.
  template<typename T>
  int copyDataOut(T* pdat) const {
    if ( ! isValid() ) return 1;
    Index ndat = dataSize(nSamples());
    for ( Index idat=0; idat<ndat; ++idat ) {
      pdat[idat] = m_data[idat];
    }
    return 0;
  }

  // Return the global index for given dimension indices.
  // If pchk is not null, then the range of each dimension is checked.
  Index globalIndex(const IndexArray& isams, Index* pchk =nullptr) const {
    if ( ! isValid() ) {
      if ( pchk != nullptr ) *pchk = 1;
      return size();
    }
    if ( isams.size() != rank() ) {
      if ( pchk != nullptr ) *pchk = 2;
      return size();
    }
    if ( pchk != nullptr ) {
      for ( Index idim=0; idim<rank(); ++idim ) {
        if  ( ! inRange(idim, isams[idim]) ) {
          *pchk = 10 + idim;
          return size();
        }
      }
      *pchk = 0;
    }
    Index idat = isams[0];
    for ( Index idim=1; idim<rank(); ++idim ) {
      idat *= m_nsams[idim];
      idat += isams[idim];
    }
    return idat;
  }
  
  // Return the value for indices isams.
  // The index arrays are checked first iff pchk != nullptr.
  // If any of these checks fail or the calculated index is out of range,
  // badValue(0 is returned.
  Float value(const IndexArray& isams, Index* pchk =nullptr) const {
    Index idat = globalIndex(isams, pchk);
    if ( idat >= size() ) return this->badValue();
    return m_data[idat];
  }

  // Set the value for indices isams.
  // Set any term, i.e. idat < nSample.
  // The index arrays are checked if pchk != nullptr.
  // If any of these checks fail or the calculated index is out of range,
  // the value is not set and the size of the data is returned.
  // If the value is set, the global index is returned.
  // That index is guranteed to be less than the data size.
  Index setValue(const IndexArray& isams, Float val, Index* pchk =nullptr) {
    Index idat = globalIndex(isams, pchk);
    if ( pchk != nullptr && *pchk > 0 ) return size();
    if ( idat >= size() ) return size();
    m_data[idat] = val;
    return idat;
  }

  // Return the total power.
  F power() const {
    F pwr = 0.0;
    IndexArray idats;
    for ( idats[0]=0; idats[0]<size(0); ++idats[0] ) {
      for ( idats[1]=0; idats[1]<size(1); ++idats[1] ) {
        F val = value(idats);
        pwr += val*val;
      }
    }
    return pwr;
  }

private:

  IndexArray m_nsams;
  DataVector m_data;

};

//**********************************************************************

#endif
