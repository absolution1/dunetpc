// RealDftNormalization.h
//
// David Adams
// April 2019, December 2020
//
// This class provides an interface to specify and interpret the normalization
// of one dimension of a DFT (discrete Fourier transform) of real data.
//
// The 1D DFT for N samples x has N complex entries aka frequency components. With standard
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
// We define a compact representation that omits the later aliased entries. There are two
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
// over all terms in the long representation or a double counting of aliased terms in the
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

#ifndef RealDftNormalization_H
#define RealDftNormalization_H

#include <vector>
#include <cmath>

class RealDftNormalization {

public:

  using Index = unsigned int;

  enum GlobalNormalization { InvalidGlobalNormalization=0, Standard=1, Consistent=2, Bin=3 };
  enum TermNormalization { InvalidTermNormalization=0, Unit=1, Power=2 };

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
  
  // Ctors from the data that specifies the two normalizations.
  // The data can be the two enums, the corresponding int or a single packed int:
  // 10*term + global.
  RealDftNormalization(GlobalNormalization gnorm, TermNormalization tnorm)
  : m_global(gnorm), m_term(tnorm) { }
  RealDftNormalization(Index ignorm, Index itnorm)
  : m_global(getGlobalNormalization(ignorm)), m_term(getTermNormalization(itnorm)) { }
    // Conversion from a single index: 10*term + global.
  explicit RealDftNormalization(Index inorm) : RealDftNormalization(inorm%10, inorm/10) { }

  // Return the global normalization.
  GlobalNormalization globalNormalization() const { return m_global; }
  bool isStandard() const   { return globalNormalization() == Standard; }
  bool isConsistent() const { return globalNormalization() == Consistent; }
  bool isBin() const        { return globalNormalization() == Bin; }

  // Return the term normalization.
  TermNormalization termNormalization() const { return m_term; }
  bool isUnit() const { return termNormalization() == Unit; }
  bool isPower() const  { return termNormalization() == Power; }

  // Check normalization is valid.
  bool hasValidGlobalNormalization() const { return isStandard() || isConsistent() || isBin(); }
  bool hasValidTermNormalization()   const { return isUnit() || isPower(); }
  bool isValid() const { return hasValidGlobalNormalization() && hasValidTermNormalization(); }

  // Normalization to use for convolution.
  static RealDftNormalization convolutionNormalization() { return RealDftNormalization(Standard, Unit); }

private:  // data

  GlobalNormalization m_global = InvalidGlobalNormalization;
  TermNormalization m_term = InvalidTermNormalization;

};

//**********************************************************************

#endif
