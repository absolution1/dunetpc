// DuneFFT.h

// David Adams
// April 2019
//
// This utility provides wrappers for performing forward and backward DFT (discrete
// Fourier transform) of real (time-domain) data. It uses the the Root TVirtualFFT
// interface which in turn uses FFTW by default.
//
// The real data is held in a vector of floats. The DFT is a vector of complex values
// of the same length and so has a factor of two redundancy. The DFT data is returned
// in an object with the RealDftData interface that provkdes methods to access the real,
// imaginary, amplitude and phase parts of the terms in the complex vector.
//
// The concrete type for the returned data is
//   CompactRealDftData<float>

#ifndef DuneFFT_H
#define DuneFFT_H

#include "dune/DuneCommon/Utility/CompactRealDftData.h"

class DuneFFT {

public:

  using Index = unsigned int;
  using DFT = CompactRealDftData<float>;
  using FloatVector = std::vector<float>;

  // Forward transform: real data (ntick starting at psam[0] --> complex freqs).
  //   ntick - # ticks to use in transform.
  //   psam - Address of the first element in the data array
  //   dft - the DFT
  //   logLevel - 0=silent, 1=init, 2=each event, >2=more
  static int fftForward(Index ntick, const float* psam, DFT& dft, Index logLevel =0);

  // Same for a full sample vector sams.
  static int fftForward(const FloatVector& sams, DFT& dft, Index logLevel =0);

  // Inverse transform: complex freqs mags, phases --> real data sams
  // Thr real and imag freq component are also recorded in xres, xims
  static int fftInverse(const DFT& dft, FloatVector& sams, Index logLevel =0);

};

#endif
