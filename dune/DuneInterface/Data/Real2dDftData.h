// Real2dDftData.h
//
// David Adams
// January 2021
//
// This class provides an interface and partial implementation for the data describing
// a 2D DFT (discrete Fourier transform) of real data. Subclasses add the
// normalization and DFT data.
//
// It is templated so DFT elements can be stored at different levels of
// floating point precision.

#ifndef Real2dDftData_H
#define Real2dDftData_H

#include "dune/DuneInterface/Data/RealDftNormalization.h"
#include <complex>
#include <array>

template<typename F>
class Real2dDftData {

public:

  using Float = F;
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

  // Dtor.
  virtual ~Real2dDftData() =default;

  // Clear the DFT, i.e. zero the # samples.
  // The normalization is retained.
  virtual void clear() =0;

  // Reset to (nsamx, nsamy) samples and zero the DFT.
  virtual void reset(const IndexArray& nsams) =0;

  // Return the normalization.
  // Same is used for all dimensions.
  virtual Norm normalization() const =0;

  // The number of dimensions.
  Index rank() const { return 2; }

  // Number of samples in each dimension.
  // nSample = size = # samples
  virtual const IndexArray& nSamples() const =0;
  virtual Index nSample(Index idim) const {
    if ( idim > rank() ) return 0;
    return nSamples()[idim];
  }
  Index size(Index idim) const { return nSample(idim); }

  // Return the status for a frequency index.
  bool inRange(Index idim, Index ifrq) const { return idim < rank() && ifrq < size(idim); }
  bool isZero(Index ifrq) const { return ifrq == 0; }
  bool isNyquist(Index idim, Index ifrq) const { return idim < rank() && 2*ifrq == size(idim); }

  // Check if this object is valid.
  // Object must have an assigned size so that the frequency status methods work.
  virtual bool isValid() const {
    return normalization().isValid() && dataSize(nSamples()) > 0;
  }

  // Bad complex.
  // Test real and imaginary pars with isnan or self != self;
  Complex badValue() const {
    double badDouble = std::numeric_limits<double>::quiet_NaN();
    return Complex(badDouble, badDouble);
  }

  // Return any term, i.e. ifrq < nSample.
  // Values depend on global and term normalizations.
  virtual Complex value(const IndexArray& ifrqs) const =0;
  virtual F power(    const IndexArray& ifrqs) const { return std::norm(value(ifrqs)); };
  virtual F amplitude(const IndexArray& ifrqs) const { return sqrt(power(ifrqs)); };
  virtual F phase(    const IndexArray& ifrqs) const { return std::arg(value(ifrqs)); };
  virtual F real(     const IndexArray& ifrqs) const { return std::real(value(ifrqs)); };
  virtual F imag(     const IndexArray& ifrqs) const { return std::imag(value(ifrqs)); };

  // Return the total power.
  F power() const {
    F pwr = 0.0;
    IndexArray ifrqs;
    for ( ifrqs[0]=0; ifrqs[0]<size(0); ++ifrqs[0] ) {
      for ( ifrqs[1]=0; ifrqs[1]<size(1); ++ifrqs[1] ) {
        F amp = amplitude(ifrqs);
        pwr += amp*amp;
      }
    }
    return pwr;
  }

};

//**********************************************************************

#endif
