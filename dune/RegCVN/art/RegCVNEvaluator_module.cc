////////////////////////////////////////////////////////////////////////
// \file    RegCVNEvaluator_module.cc
// \brief   Producer module creating RegCVN neural net results modified from CVNEvaluator_module.cc
// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// ROOT includes
#include "TFile.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVectorD.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"

// NOvASoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"


#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dune/RegCVN/func/RegCVNResult.h"
#include "dune/RegCVN/func/RegPixelMap.h"
#include "dune/RegCVN/art/TFRegNetHandler.h"

namespace cvn {

  class RegCVNEvaluator : public art::EDProducer {
  public:
    explicit RegCVNEvaluator(fhicl::ParameterSet const& pset);
    ~RegCVNEvaluator();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:

    /// Module label for input pixel maps
    std::string fPixelMapInput;
    std::string fResultLabel;

    std::string fCVNType;

    cvn::TFRegNetHandler fTFHandler;

    /// Number of outputs fron neural net
    unsigned int fNOutput;

  };

  //.......................................................................
  RegCVNEvaluator::RegCVNEvaluator(fhicl::ParameterSet const& pset):
    fPixelMapInput (pset.get<std::string>         ("PixelMapInput")),
    fResultLabel (pset.get<std::string>         ("ResultLabel")),
    fCVNType     (pset.get<std::string>         ("CVNType")),
    fTFHandler       (pset.get<fhicl::ParameterSet> ("TFNetHandler"))
  {
    produces< std::vector<cvn::RegCVNResult>   >(fResultLabel);

  }
  //......................................................................
  RegCVNEvaluator::~RegCVNEvaluator()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void RegCVNEvaluator::beginJob()
  {  }

  //......................................................................
  void RegCVNEvaluator::endJob()
  {
  }

  //......................................................................
  void RegCVNEvaluator::produce(art::Event& evt)
  {

    /// Define containers for the things we're going to produce
    std::unique_ptr< std::vector<RegCVNResult> >
                                  resultCol(new std::vector<RegCVNResult>);

    /// Load in the pixel maps
    art::Handle< std::vector< cvn::RegPixelMap > > pixelmapListHandle;
    std::vector< art::Ptr< cvn::RegPixelMap > > pixelmaplist;
    if (evt.getByLabel(fPixelMapInput, "regcvnmap", pixelmapListHandle))
      art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);


    /// Make sure we have a valid name for the CVN type
    if(fCVNType == "TF" || fCVNType == "Tensorflow" || fCVNType == "TensorFlow"){
      // If we have a pixel map then use the TF interface to give us a prediction
      if(pixelmaplist.size() > 0){

        std::vector<float> networkOutput = fTFHandler.Predict(*pixelmaplist[0]);

        // cvn::Result can now take a vector of floats and works out the number of outputs
        resultCol->emplace_back(networkOutput);
      }
    }
    else{
      mf::LogError("RegCVNEvaluator::produce") << "CVN Type not in the allowed list: Tensorflow" << std::endl;
      mf::LogError("RegCVNEvaluator::produce") << "Exiting without processing events" << std::endl;
      return;
    }


    evt.put(std::move(resultCol), fResultLabel);

  }

  DEFINE_ART_MODULE(cvn::RegCVNEvaluator)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







