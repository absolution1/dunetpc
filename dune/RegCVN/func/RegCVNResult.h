////////////////////////////////////////////////////////////////////////
/// \file    RegCVNResult.h
/// \brief   RegCVNResult for RegCVN modified from Result.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef CVN_RESULT_H
#define CVN_RESULT_H

#include <vector>
#include "dune/RegCVN/func/InteractionType.h"

namespace cvn
{
  /// RegCVNResult, basic output of CVN neural net
  class RegCVNResult
  {
  public:
    RegCVNResult(const float* output, unsigned int& nOutputs);
    // Vector version of the constructor
    RegCVNResult(const std::vector<float> output);
    RegCVNResult();

    /// Index of maximum value in vector
    unsigned int ArgMax();

    /// Maximum value in vector
    float Max();

    /// Return the predicted interaction type
    TFResultType PredictedInteractionType();

    /// Return the numu flavour probability
    float GetNumuProbability();

    /// Return the nue flavour probability
    float GetNueProbability();

    /// Return the nutau flavour probability
    float GetNutauProbability();

    /// Return the NC probability
    float GetNCProbability();

    /// Number of outputs, i.e. size of vector
    unsigned int NOutput();

    std::vector<float> fOutput;  ///< Vector of outputs from neural net

  };
}

#endif  // CVN_RESULT_H

