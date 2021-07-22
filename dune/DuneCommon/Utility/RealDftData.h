// RealDftData.h
//
// David Adams
// April 2019
//
// This class provides an interface and partial implementation for the data describing
// a 1D DFT (discrete Fourier transform) of real data. Subclasses add the
// normalization and DFT data.
//
// It is templated so DFT elements can be stored at different levels of
// floating point precision.

#ifndef RealDftData_H
#define RealDftData_H

#include "dune/DuneInterface/Data/RealDftNormalization.h"

template<typename F>
class RealDftData {

public:

  using Float = F;

  using Index = unsigned int;

  // Return the normalization.
  virtual const RealDftNormalization& normalization() const =0;

  // Check and return dimension information.
  // nSample = size = # samples
  // nCompact = # amplitudes in the compact representation
  // Subclass must override nSample and may override nCompact.
  virtual Index nSample() const =0;
  Index size() const { return nSample(); }
  virtual Index nCompact() const { return nSample() ? nSample()/2 + 1 : 0; }
  bool isEven() const { return (nSample() + 1) % 2; }
  bool isOdd() const { return nSample() % 2; }

  // Check if this object is valid.
  // Object must have an assigned size so that the frequency status methods work.
  virtual bool isValid() const { return normalization().isValid() && size() > 0; }

  // Clear the DFT, i.e. zero the # samples.
  // The normalization is retained.
  virtual void clear() =0;

  // Reset to nsam samples and zero the DFT.
  virtual void reset(Index nsam) =0;

  // Return the status for a frequency index.
  bool inRange(Index ifrq) const { return ifrq < size(); }
  bool isZero(Index ifrq) const { return ifrq == 0; }
  bool isNyquist(Index ifrq) const { return 2*ifrq == size(); }
  bool isNotAliased(Index ifrq) const { return isZero(ifrq) || isNyquist(ifrq); }
  bool isAliased(Index ifrq) const { return ifrq < size() && !isZero(ifrq) && !isNyquist(ifrq); }

  // Value to return if caller requests a frequency out of range.
  virtual F badValue() const { return 0.0; }

  // Return any term, i.e. ifrq < nSample.
  // Values depend on global and term normalizations.
  virtual F amplitude(Index ifrq) const =0;
  virtual F phase(Index ifrq) const =0;
  virtual F real(Index ifrq) const =0;
  virtual F imag(Index ifrq) const =0;

  // Return the amplitude with convolution normalization, i.e.
  // standard global and simple term normalizations.
  F convAmplitude(Index ifrq) const {
    if ( ifrq >= size() ) return badValue();
    if ( ! isValid() ) return badValue();
    F amp = amplitude(ifrq);
    F len = size();
    const F osq2 = 1.0/sqrt(2.0);
    if ( normalization().isConsistent() ) amp *= sqrt(len);
    if ( normalization().isBin() ) amp *= len;
    if ( normalization().isPower() && isAliased(ifrq) ) amp *= osq2;
    return amp;
  }

  // Return the compact power for any compact term, i.e. adding that from its alias.
  // The sum over the range [0, nCompact) gives the total power.
  F compactPower(Index ifrq) const {
    if ( ifrq >= nCompact() ) return 0.0;
    F amp = amplitude(ifrq);
    F pwr = amp*amp;
    if ( normalization().isUnit() && isAliased(ifrq) ) pwr *= 2.0;
    if ( normalization().isStandard() ) pwr /= size();
    if ( normalization().isBin() ) pwr *= size();
    return pwr;
  }

  // Return the total power.
  F power() const {
    F pwr = 0.0;
    for ( Index ifrq=0; ifrq<nCompact(); ++ifrq ) {
      pwr += compactPower(ifrq);
    }
    return pwr;
  }

};

//**********************************************************************

#endif
