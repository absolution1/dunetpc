////////////////////////////////////////////////////////////////////////
/// \file    TFRegNetHandler.h
/// \brief   TFRegNetHandler for RegCVN modified from TFNetHandler.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCVN_TFNETHANDLER_H
#define REGCVN_TFNETHANDLER_H

#include <vector>
#include <memory>

#include "dune/RegCVN/func/RegPixelMap.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/RegCVN/func/RegCVN_TF_Graph.h"
//#include "larreco/RecoAlg/ImagePatternAlgs/Tensorflow/TF/tf_graph.h"

namespace cvn
{

  /// Wrapper for caffe::Net which handles construction and prediction
  class TFRegNetHandler
  {
  public:

    /// Constructor which takes a pset with DeployProto and ModelFile fields
    TFRegNetHandler(const fhicl::ParameterSet& pset);

    /// Return prediction arrays for RegPixelMap
    std::vector<float> Predict(const RegPixelMap& pm);

    std::vector<float> PredictNuEEnergy(const RegPixelMap& pm);

  private:

    std::string  fLibPath;  ///< Library path (typically dune_pardata...)
    std::string  fTFProtoBuf;  ///< location of the tf .pb file in the above path
    unsigned int fInputs;   ///< Number of tdcs for the network to classify
    std::vector<bool> fReverseViews; ///< Do we need to reverse any views?
    std::unique_ptr<tf::RegCVNGraph> fTFGraph; ///< Tensorflow graph

  };

}

#endif  // CVN_TFNETHANDLER_H
