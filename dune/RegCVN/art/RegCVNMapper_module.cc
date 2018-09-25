////////////////////////////////////////////////////////////////////////
// \file    RegCVNMapper_module.cc
// \brief   Producer module for creating RegCVN PixelMap objects modified from CVNMapper_module.cc
// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// ROOT includes
#include "TTree.h"
#include "TH2F.h"

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
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dune/RegCVN/art/RegPixelMapProducer.h"
#include "dune/RegCVN/func/RegPixelMap.h"


namespace cvn {

  class RegCVNMapper : public art::EDProducer {
  public:
    explicit RegCVNMapper(fhicl::ParameterSet const& pset);
    ~RegCVNMapper();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:
    /// Module lablel for input clusters
    std::string    fHitsModuleLabel;

    /// Instance lablel for cluster pixelmaps
    std::string    fClusterPMLabel;

    /// Minimum number of hits for cluster to be converted to pixel map
    unsigned short fMinClusterHits;

    /// Width of pixel map in tdcs
    unsigned short fTdcWidth;

    /// Length of pixel map in wires
    unsigned short fWireLength;

    /// Length of pixel map in wires
    double fTimeResolution;

    /// Maximum gap in wires at front of cluster to prevent pruning of upstream
    /// hits
    unsigned int fMaxWireGap;

    /// Use unwrapped pixel maps?
    bool fUnwrappedPixelMap;

    /// PixelMapProducer does the work for us
    RegPixelMapProducer fProducer;

  };



  //.......................................................................
  RegCVNMapper::RegCVNMapper(fhicl::ParameterSet const& pset):
  fHitsModuleLabel  (pset.get<std::string>    ("HitsModuleLabel")),
  fClusterPMLabel(pset.get<std::string>    ("ClusterPMLabel")),
  fMinClusterHits(pset.get<unsigned short> ("MinClusterHits")),
  fTdcWidth     (pset.get<unsigned short> ("TdcWidth")),
  fWireLength   (pset.get<unsigned short> ("WireLength")),
  fTimeResolution   (pset.get<unsigned short> ("TimeResolution")),
  fProducer      (fWireLength, fTdcWidth, fTimeResolution)
  {

    produces< std::vector<cvn::RegPixelMap>   >(fClusterPMLabel);

  }

  //......................................................................
  RegCVNMapper::~RegCVNMapper()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void RegCVNMapper::beginJob()
  {  }

  //......................................................................
  void RegCVNMapper::endJob()
  {
  }

  //......................................................................
  void RegCVNMapper::produce(art::Event& evt)
  {

    art::Handle< std::vector< recob::Hit > > hitListHandle;
    std::vector< art::Ptr< recob::Hit > > hitlist;
    if (evt.getByLabel(fHitsModuleLabel, hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
    art::FindManyP<recob::Wire> fmwire(hitListHandle, evt, fHitsModuleLabel);

    unsigned short nhits = hitlist.size();

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cvn::RegPixelMap> >
      pmCol(new std::vector<cvn::RegPixelMap>);

    if(nhits>fMinClusterHits){
      RegPixelMap pm = fProducer.CreateMap(hitlist, fmwire);
      // skip if PixelMap is empty
      if (pm.fInPM) pmCol->push_back(pm);
      //pm.Print();
    }
    //Boundary bound = pm.Bound();
    //}
    evt.put(std::move(pmCol), fClusterPMLabel);
    //std::cout<<"Map Complete!"<<std::endl;
  }

  //----------------------------------------------------------------------



DEFINE_ART_MODULE(cvn::RegCVNMapper)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







