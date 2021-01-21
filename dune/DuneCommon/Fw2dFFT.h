// Fw2dFFT.h

// David Adams
// January 2021
//
// This utility provides wrappers for performing forward and backward DFT (discrete
// 2D Fourier transform) of real (time-domain) data using FFTW.
//
// Separate forward and backward FFTW plans are created for each array of data dimensions.
// These may be created in advance by calling fowardPlan or backwardPlan or are created
// in the first call to perform an FFT with that set of dimentsions.
//
// Space to hold the transformation data is allocated in the constructor. Use
// checkDataSize to check if the space is sufficient for given data dimensions.
//
// The real data is held in a vector of floats. The DFT is a vector of complex values
// of the same length and so has a factor of two redundancy. The DFT data is returned
// in an object with the Real2dDftData interface.
//
// The concrete type for the returned data is
//   FftwReal2dDftData<double>

#ifndef Fw2dFFT_H
#define Fw2dFFT_H

#include "dune/DuneInterface/Data/Real2dData.h"
#include "dune/DuneCommon/FftwReal2dDftData.h"
#include <map>

class Fw2dFFT {

public:

  using DataFloat = float;
  using DftFloat = double;      // Type used in FFTW
  using Data = Real2dData<DataFloat>;
  using DFT = FftwReal2dDftData<DftFloat>;
  using Index = DFT::Index;
  using IndexArray = DFT::IndexArray;
  using Complex = DFT::Complex;
  using Plan = fftw_plan;
  using PlanMap = std::map<IndexArray, Plan>;
  typedef double FftwComplex[2];

  // Ctor from maximum data size and optimization.
  // For n1xn2 data, ndat must be at least as big as n1*n2 and 2*n1*(n2/2+1)
  // opt = 0-2 (FFTW_ESTIMATE, FFTW_PLAN, FFTW_PATIENT)
  //       Larger numbers take longer to build plans but are faster for each transform.
  Fw2dFFT(Index ndatMax, Index opt);

  // Dtor. Frees the data caches and destroys the plans.
  ~Fw2dFFT();

  // Check an array of sample sizes.
  //   0 - Ok for processing here.
  //   1 - Insufficent size for input data.
  //   2 - Insufficient size for output (DFT) data.
  //   3 - One or more dimensions has size zero.
  bool checkDataSize(const IndexArray& nsams) const;

  // Return if we have  plan for given data sizes.
  bool haveForwardPlan(const IndexArray& nsams) const;
  bool haveBackwardPlan(const IndexArray& nsams) const;

  // Return the plan for a data size.
  // The plan is created if not already existing.
  Plan& forwardPlan(const IndexArray& nsams);
  Plan& backwardPlan(const IndexArray& nsams);

  // Forward transform: real data (ntick starting at psam[0] --> complex freqs).
  //   psam - Address of the first element in input the data array
  //   dft - Output DFT data. Dimension sizes are obtained from this.
  //   logLevel - 0=silent, 1=init, 2=each event, >2=more
  // Returns 0 for success.
  int fftForward(const Data& dat, DFT& dft, Index logLevel =0);

  // Inverse transform: complex freqs mags, phases --> real data sams
  // The real and imag freq component are also recorded in xres, xims
  // Returns 0 for success.
  int fftBackward(const DFT& dft, Data& dat, Index logLevel =0);

  // Return the DFT data as FFTW complex.
  FftwComplex* fftwOutData() { return reinterpret_cast<FftwComplex*>(m_outData); }

  // Return the DFT data as floating point: real0, imag0, real1, ...
  DftFloat* floatOutData() { return reinterpret_cast<DftFloat*>(m_outData); }

private:

  Index m_ndatMax;
  Index m_flag;
  DftFloat* m_inData;
  Complex* m_outData;
  PlanMap m_forwardPlans;
  PlanMap m_backwardPlans;

};

#endif
