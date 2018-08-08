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

  unsigned int RegCVNResult::ArgMax(){
    // Get the max element iterator and convert to vector index
    return std::distance(fOutput.begin(),std::max_element(fOutput.begin(),fOutput.end()));
  }

  float RegCVNResult::Max(){
    // Get the maximum value by dereferencing the iterator
    return *std::max_element(fOutput.begin(),fOutput.end());
  }

  unsigned int RegCVNResult::NOutput(){
    return fOutput.size();
  }

  /// Return the predicted interaction type
  TFResultType RegCVNResult::PredictedInteractionType(){
    return static_cast<TFResultType>(this->ArgMax());
  }

  /// Return the numu flavour probability
  float RegCVNResult::GetNumuProbability(){
    return fOutput[TFResultType::kTFNumuQE] + fOutput[TFResultType::kTFNumuRes]
         + fOutput[TFResultType::kTFNumuDIS] + fOutput[TFResultType::kTFNumuOther];
  }
  
  /// Return the nue flavour probability
  float RegCVNResult::GetNueProbability(){
    return fOutput[TFResultType::kTFNueQE] + fOutput[TFResultType::kTFNueRes]
         + fOutput[TFResultType::kTFNueDIS] + fOutput[TFResultType::kTFNueOther];
  }

  /// Return the nutau flavour probability
  float RegCVNResult::GetNutauProbability(){
    return fOutput[TFResultType::kTFNutauQE] + fOutput[TFResultType::kTFNutauRes]
         + fOutput[TFResultType::kTFNutauDIS] + fOutput[TFResultType::kTFNutauOther];
  }

  /// Return the NC probability
  float RegCVNResult::GetNCProbability(){

    // The old caffe network didn't give us an NC probability
    // So make sure we have enough values to grab it
    float result = -999;

    if(fOutput.size() > static_cast<unsigned int>(TFResultType::kTFNC)){
      result = fOutput[TFResultType::kTFNC];
    }
    else{
      mf::LogError("cvn::RegCVNResult") << "Output vector too short to include an NC probability" << std::endl;
    }

    return result;
  }

}
