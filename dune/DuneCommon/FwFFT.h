// FwFFT.h

// David Adams
// January 2021
//
// This utility provides wrappers for performing forward and backward DFT (discrete
// Fourier transform) of real (time-domain) data using FFTW.
//
// The real data is held in a vector of floats. The DFT is a vector of complex values
// of the same length and so has a factor of two redundancy. The DFT data is returned
// in an object with the RealDftData interface that provides methods to access the real,
// imaginary, amplitude and phase parts of the terms in the complex vector.
//
// The concrete type for the returned data is
//   CompactRealDftData<float>

#ifndef FwFFT_H
#define FwFFT_H

#include "dune/DuneCommon/CompactRealDftData.h"
#include "fftw3.h"
#include <map>

class FwFFT {

public:

  using Index = unsigned int;
  using Float = double;
  using Complex = fftw_complex;
  using DFT = CompactRealDftData<float>;
  using FloatVector = std::vector<float>;
  using Plan = fftw_plan;
  using PlanMap = std::map<Index, Plan>;

  // Ctor from maximum data size and optimization.
  // opt = 0-2 (FFTW_ESTIMATE, FFTW_PLAN, FFTW_PATIENT)
  //       Larger nubers take longer to build plans but are faster for each transform.
  FwFFT(Index nsamMax, Index opt);

  // Dtor. Frees the data caches and destroys the plans.
  ~FwFFT();

  // Return the plan for a data size.
  // The plan is cretaed if not already existing.
  Plan& forwardPlan(Index nsam);
  Plan& backwardPlan(Index nsam);

  // Return the plan for a given size.
  // If not already existing, the plan is created.
  // Forward transform: real data (ntick starting at psam[0] --> complex freqs).
  //   ntick - # ticks to use in transform.
  //   psam - Address of the first element in the data array
  //   dft - the DFT
  //   logLevel - 0=silent, 1=init, 2=each event, >2=more
  int fftForward(Index nsam, const float* psam, DFT& dft, Index logLevel =0);

  // Same for a full sample vector sams.
  int fftForward(const FloatVector& sams, DFT& dft, Index logLevel =0);

  // Inverse transform: complex freqs mags, phases --> real data sams
  // The real and imag freq component are also recorded in xres, xims
  int fftInverse(const DFT& dft, FloatVector& sams, Index logLevel =0);

private:

  Index m_nsamMax;
  Index m_flag;
  double* m_inData;
  Complex* m_outData;
  PlanMap m_forwardPlans;
  PlanMap m_backwardPlans;

};

#endif
