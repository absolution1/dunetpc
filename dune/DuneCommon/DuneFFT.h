// DuneFFT.h

// David Adams
// April 2019
//
// Utility perform forward and backward FFT.
//
// Normalization options:
//     0 - Unnormalized ==> integrated raw power of transform is N times larger than input.
//         Normalization of 1/N is applied on inverse to recover original waveform.
//     1 - Normalized ==> Integrated power of transform or inverse is the same as its input.
//         I.e. normalization of 1/sqrt(N) is applied both forward and backward.
//     2 - Bin normalized ==> Integrated power of transform is 1/N times input.
//         Normalization of 1/N on forward transform and one on inverse.

#ifndef DuneFFT_H
#define DuneFFT_H

#include <vector>

class DuneFFT {

public:

  using Index = unsigned int;
  using FloatVector = std::vector<float>;

  // Forward transform: real data (ntick starting at psam[0] --> complex freqs).
  //   normOpt - Normalization (see above).
  //   ntick - # ticks to use in transform.
  //   psam - Address of the first element in the data array
  //   xres, xims - Real and imaginary componets of the transformation.
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
