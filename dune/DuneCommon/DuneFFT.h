// DuneFFT.h

// David Adams
// April 2019
//
// This utility provides wrappers for performing forward and backward DFT (discrete
// Fourier transform) of real (time-domain) data. It uses the the Root TVirtualFFT
// interface which in turn uses FFTW by default.
//
// For a real (time-domain) data of length N, the DFT, i.e the (frequency-domain)
// transform has N complex entries X. This is called the long representation because
// 2*N values are required to hold the real and imaginary parts of this vector.
// For convenience, the methods here return these in two vectors (of lenght N).
//
// The redundancy in the these 2*N values follows from the time-domain data being real.
// This is reflected in the mirror (or alias) relation X[N-i] = X*[i] where
// x* denotes complex conjugation.
//
// Here, the DFT data is (also) expressed in a second compact representation that omits
// the later aliased entries. There are two shorter vectors with one holding magnitudes
// and the second holding phases.
// For N even, the magnitude vector has length N/2+1 and the phase vector is length N/2.
// For N odd, both vectors have the same length: N/2+1 = (N+1)/2.
// Equivalently for any N, the respective lengths are N/2+1 and (N+1)/2.
// Thus the length of the original data may be deduced from the lengths of the
// magnitude and phase vectors.
// The first phase is always zero or pi (the first term is real) leaving N significant values
// for N even or odd.
//
// There is ambiguity in the sign of the magnitude entries as flipping the sign is equivalent to
// adding pi to the correponding phase. Thus, for most of the entries, then magnitude may
// be chosen positive. The one exception is the last entry for N even. There is no corresponding
// phase and so that "magnitude" may be forced to be negative.
//
// The values in the long representation may be obtained from those in the compact.
// Denoting the real and imaginary vectors by real and imag and the magnitude and phase
// vectors by mag and pha with respective lengths Nmag and Npha:
//   real[0] =   mag[0]*cos(pha[0]),     imag[0] = 0                                      
//   real[i] = K*mag[i]*cos(pha[i]),     imag[i] = K*mag[i]*sin(pha[i])        i < Npha
//   real[i] =   mag[i]*cos(pha[i]),     imag[i] = 0                           i = Npha = N/2, N even
//   real[i] = K*mag[N-i]*cos(pha[N-i]), imag[i] = -K*mag[N-i]*sin(pha[N-i])   i > Npha
//
// Note that the first term and, for N even, the last "magnitude" may be negative.
//
// We support two term normalizations in the above expressions:
// normOptTerm:
//     1 - Simple: K = 1
//     2 - Power:  K = 1/sqrt(2)
// In the former case, convolution or deconvolution may be accomplished by directly multiplying
// or dividing the terms in the compact representation but a calculation of the power must
// double the contribution from terms with aliases.
// In the latter case, the power is obtained directly from the (square of) the magnitude for any
// frequency index but the mirrored terms must be corrected for double counting of K if both
// DFTs use the power normalization.
//
// There are also a few natural conventions for global normalization of the transformed data:
// normOptGlobal:
//     0 - Expanded ==> integrated raw power of transform is N times larger than input.
//         Normalization of 1/N is applied on inverse to recover original waveform.
//     1 - Consistent ==> Integrated power of transform or inverse is the same as its input.
//         I.e. normalization of 1/sqrt(N) is applied both forward and backward.
//     2 - Bin ==> Integrated power of transform is 1/N times input.
//         Normalization of 1/N on forward transform and one on inverse.
//
// Note, that unlike the term normalization, this normalization is relevant to both both the long
// and compact representations.
//
// The methods here take a flag that merges these options:
//   normOpt = 10*normOptTerm + normOptGlobal
//
// Allowed values are 10, 11, 12, 20, 12 and 22.
//
// We anticpate data transforms will use consistent-power normalization (21) and that convolution
// transforms (response, filter, smearing, etc) will use consistent-simple normalization (11).
// This way power can be obtained directly from data transforms and the (de)covolution of data
// with a convolution function is obtained by direct multiplication (division) of the compact
// entries.

#ifndef DuneFFT_H
#define DuneFFT_H

#include <vector>

class DuneFFT {

public:

  using Index = unsigned int;
  using FloatVector = std::vector<float>;

  // Return the global part of a normalization flag.
  static Index getNormOptGlobal(Index normOpt) {
    if ( normOpt < 10 ) return 99;
    if ( normOpt > 22 ) return 99;
    Index opt = normOpt % 10;
    if ( opt > 2 ) return 99;
    return opt;
  }

  // Return the term part of a normalization flag.
  static Index getNormOptTerm(Index normOpt) {
    if ( normOpt < 10 ) return 99;
    if ( normOpt > 22 ) return 99;
    Index opt = normOpt / 10;
    if ( opt > 2 ) return 99;
    return opt;
  }
  static bool isSimple(Index normOpt) { return getNormOptTerm(normOpt) == 1; }
  static bool isPower(Index normOpt) { return getNormOptTerm(normOpt) == 2; }

  // Forward transform: real data (ntick starting at psam[0] --> complex freqs).
  //   normOpt - Normalization (see above).
  //   ntick - # ticks to use in transform.
  //   psam - Address of the first element in the data array
  //   xres, xims - Real and imaginary components of the transformation.
  //   xmgs, xphs - Magnitudes and phases for the transformation.
  //   logLevel - 0=silent, 1=init, 2=each event, >2=more
  static int fftForward(Index normOpt, Index ntick, const float* psam,
                        FloatVector& xres, FloatVector& xims, FloatVector& mags, FloatVector& phases,
                        Index logLevel =0);

  // Same for a full sample vector sams.
  static int fftForward(Index normOpt, const FloatVector& sams,
                        FloatVector& xres, FloatVector& xims, FloatVector& mags, FloatVector& phases,
                        Index logLevel =0);

  // Inverse transform: complex freqs mags, phases --> real data sams
  // Thr real and imag freq component are also recorded in xres, xims
  static int fftInverse(Index normOpt, const FloatVector& mags, const FloatVector& phases,
                        FloatVector& xres, FloatVector& xims, FloatVector& sams,
                        Index logLevel =0);

};

#endif
