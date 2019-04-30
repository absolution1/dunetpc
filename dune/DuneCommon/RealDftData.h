// RealDftData.h
//
// David Adams
// April 2019
//
// This class provides an interface and partial implementation for the data describing
// a DFT (discrete Fourier transform) of real data. Subclasses add the DFT data.
//
// The DFT for N samples x has N complex entries aka frequency components. With standard
// normalization (see below) the kth component of the DFT is
//   X_k = SUM_m x_m exp(-2*pi*i*m*k/N)
// Here, the 2N values that specify the real and imaginary components of X are called
// the long representation of the DFT because only N values are needed in the case of
// interest here where the original data x is purely real.
//
// The redundancy in the these 2*N values is reflected in the mirror (or alias) relation
//   X[N-i] = X*[i]
// where X* denotes complex conjugation.
//
// We define the compact representation that omits the later aliased entries. There are two
// shorter vectors with one holding amplitudes and the second holding phases.
// For N even, the amplitude vector has length N/2+1 and the phase vector is length N/2.
// For N odd, both vectors have the same length: N/2+1 = (N+1)/2.
// Equivalently for any N, the respective lengths are N/2+1 and (N+1)/2.
// Thus the length of the original data may be deduced from the lengths of the
// amplitude and phase vectors.
// The first phase is always zero or pi (the first term is real) leaving N significant values
// for N even or odd.
//
// There is ambiguity in the sign of the amplitude as flipping the sign is equivalent to
// adding pi to the correponding phase. Thus, for most of the entries, then amplitude may
// be chosen positive. The one exception is the last entry for N even. There is no corresponding
// phase and so that amplitude may be forced to be negative.
//
// The values in the long representation may be obtained from those in the compact.
// Denoting the real and imaginary vectors by real and imag and the amplitude and phase
// vectors by amp and pha with respective lengths Namp and Npha:
//   real[0] =   amp[0]*cos(pha[0]),     imag[0] = 0                                      
//   real[i] = K*amp[i]*cos(pha[i]),     imag[i] = K*amp[i]*sin(pha[i])        i < Npha
//   real[i] =   amp[i]*cos(pha[i]),     imag[i] = 0                           i = Npha = N/2, N even
//   real[i] = K*amp[N-i]*cos(pha[N-i]), imag[i] = -K*amp[N-i]*sin(pha[N-i])   i > Npha
//
// We support two term normalizations in the above expressions:
// termNormalization:
//     1 - Unit: K = 1
//     2 - Power:  K = 1/sqrt(2)
//
// In the former case, the amplitude for any frequency component has the same magnitude as the
// corresponding term in the long representation. Convolution may be accomplished by directly
// multiplying terms in the compact representation. Evaluation of the total power requires a sum
// over all terms in the long representation or a doulbe counting of aliased terms in the
// compact representation.
//
// In the latter case, the amplitude for a component includes the power contribution for its
// alias, if present, i.e. the amplitude is increased by sqrt(2) for terms that have aliases.
// The total power may be obtained directly by summing the (squares of the) amplitudes in the
// compact representation. Aliased terms must be corrected for double counting if both DFTs
// have power normalization.
//
// There are also a few natural conventions for global normalization of the transformed data:
// globalNormalization:
//     1 - Standard ==> integrated raw power of transform is N times larger than input.
//         Normalization of 1/N is applied on inverse to recover original waveform.
//     2 - Consistent ==> Integrated power of transform or inverse is the same as its input.
//         I.e. normalization of 1/sqrt(N) is applied both forward and backward.
//     3 - Bin ==> Integrated power of transform is 1/N times input.
//         Normalization of 1/N on forward transform and one on inverse.
//
// Unlike the term normalization, this normalization is relevant to both both the long
// and compact representations.
//
// Complete specification of the normalization requires both of these values, here denoted
// globalNormalization-termNormalization, e.g. standard-unit for standard global normalization
// and unit term normalization.
//
// Convolution of two DFTs is simple multiplication of their terms and the result has the same
// normalization as the first if the second has standard-unit normalization aka convolution
// normalization.  Corrections to aliased or all terms are required for other normalizations.
//
// We anticipate data transforms will use consistent-power normalization and that convolution
// transforms (response, filter, smearing, etc) will use the standard-uni normalization.
// This way power can be obtained directly from data transforms and the (de)covolution of data
// with a convolution function is obtained by direct multiplication (division) of the compact
// entries.

#ifndef RealDftData_H
#define RealDftData_H

#include <vector>
#include <cmath>

//**********************************************************************

// For convenience, methods relating to the normalization are collected in this
// non-templated base class.

class RealDftNormalization {

public:

  using Index = unsigned int;

  enum GlobalNormalization { InvalidGlobalNormalization=0, Standard=1, Consistent=2, Bin=3 };
  enum TermNormalization { InvalidTermNormalization=0, Unit=1, Power=2 };

  // Subclass that holds both normalizations.
  // Provides conversion from two or one index.
  class FullNormalization {
  public:
    GlobalNormalization global = InvalidGlobalNormalization;
    TermNormalization term = InvalidTermNormalization;
    FullNormalization() =default;
    FullNormalization(GlobalNormalization gnorm, TermNormalization tnorm)
    : global(gnorm), term(tnorm) { }
    // Conversion from indices.
    FullNormalization(Index ignorm, Index itnorm)
    : global(getGlobalNormalization(ignorm)), term(getTermNormalization(itnorm)) { }
    // Conversion from a single index: 10*term + global.
    explicit FullNormalization(Index inorm) : FullNormalization(inorm%10, inorm/10) { };
  };

  // Normalization conversion from an index.
  static GlobalNormalization getGlobalNormalization(Index iarg) {
    if ( iarg == 1 ) return Standard;
    if ( iarg == 2 ) return Consistent;
    if ( iarg == 3 ) return Bin;
    return InvalidGlobalNormalization;
  }
  static TermNormalization getTermNormalization(Index iarg) {
    if ( iarg == 1 ) return Unit;
    if ( iarg == 2 ) return Power;
    return InvalidTermNormalization;
  }
  
  // Return the full normalization for this object.
  virtual FullNormalization fullNormalization() const =0;

  // Return the global normalization.
  GlobalNormalization globalNormalization() const { return fullNormalization().global; }
  bool isStandard() const   { return globalNormalization() == Standard; }
  bool isConsistent() const { return globalNormalization() == Consistent; }
  bool isBin() const        { return globalNormalization() == Bin; }

  // Return the term normalization.
  TermNormalization termNormalization() const { return fullNormalization().term; }
  bool isUnit() const { return termNormalization() == Unit; }
  bool isPower() const  { return termNormalization() == Power; }

  // Check normalization is valid.
  bool hasValidGlobalNormalization() const { return isStandard() || isConsistent() || isBin(); }
  bool hasValidTermNormalization()   const { return isUnit() || isPower(); }
  bool hasValidNormalization() const { return hasValidGlobalNormalization() && hasValidTermNormalization(); }

  // Normalization to use for convolution.
  static FullNormalization convolutionNormalization() { return FullNormalization(Standard, Unit); }

};

//**********************************************************************

// This is the class that (along with the base class) defines the interface
// (and provides some of the implementation) for accessing DFT information.
// It is templated so DFT elements can be stored at different levels of
// floating point precision.

template<typename F>
class RealDftData : public RealDftNormalization {

public:

  using Float = F;

  // Clear the DFT, i.e. zero the # samples.
  // The normalization is retained.
  virtual void clear() =0;

  // Reset to nsam samples and zero the DFT.
  virtual void reset(Index nsam) =0;

  // Check and return dimension information.
  // nSample = size = # samples
  // nCompact = # amplitudes in the compact representation
  // Subclass must override nSample and may override nCompact.
  virtual Index nSample() const =0;
  Index size() const { return nSample(); }
  virtual Index nCompact() const { return nSample() ? nSample()/2 + 1 : 0; }
  bool isEven() const { return (nSample() + 1) % 2; }
  bool isOdd() const { return nSample() % 2; }

  // Return the status for a frequency index.
  bool inRange(Index ifrq) const { return ifrq < size(); }
  bool isZero(Index ifrq) const { return ifrq == 0; }
  bool isNyquist(Index ifrq) const { return 2*ifrq == size(); }
  bool isNotAliased(Index ifrq) const { return isZero(ifrq) || isNyquist(ifrq); }
  bool isAliased(Index ifrq) const { return ifrq < size() && !isZero(ifrq) && !isNyquist(ifrq); }

  // Check if this object is valid.
  // Object must have an assigned size so that the frequency status methods work.
  virtual bool isValid() const { return this->hasValidNormalization() && size() > 0; }

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
    if ( this->isConsistent() ) amp *= sqrt(len);
    if ( this->isBin() ) amp *= len;
    if ( this->isPower() && isAliased(ifrq) ) amp *= osq2;
    return amp;
  }

  // Return the compact power for any compact term, i.e. adding that from its alias.
  // The sum over the range [0, nCompact) gives the total power.
  F compactPower(Index ifrq) const {
    if ( ifrq >= nCompact() ) return 0.0;
    F amp = amplitude(ifrq);
    F pwr = amp*amp;
    if ( isUnit() && isAliased(ifrq) ) pwr *= 2.0;
    if ( isStandard() ) pwr /= size();
    if ( isBin() ) pwr *= size();
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
