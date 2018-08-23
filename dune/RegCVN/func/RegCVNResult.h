////////////////////////////////////////////////////////////////////////
/// \file    RegCVNResult.h
/// \brief   RegCVNResult for RegCVN modified from Result.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCVN_RESULT_H
#define REGCVN_RESULT_H

#include <vector>

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

    std::vector<float> fOutput;  ///< Vector of outputs from neural net

  };
}

#endif  // REGCVN_RESULT_H

