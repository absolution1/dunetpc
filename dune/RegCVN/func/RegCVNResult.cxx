////////////////////////////////////////////////////////////////////////
/// \file    RegCVNResult.h
/// \brief   RegCVNResult for RegCVN modified from Result.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <algorithm>

#include "dune/RegCVN/func/RegCVNResult.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cvn
{

  RegCVNResult::RegCVNResult(const float* output, unsigned int& nOutputs):
  fOutput(nOutputs)
  {
    for(size_t i = 0; i < nOutputs; ++i) fOutput[i] = output[i];
  }

  RegCVNResult::RegCVNResult(const std::vector<float> output){
    fOutput = output; 
  }

  RegCVNResult::RegCVNResult():
  fOutput()
  {}

}
