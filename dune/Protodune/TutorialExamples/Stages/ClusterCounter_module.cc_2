////////////////////////////////////////////////////////////////////////
// Class:       ClusterCounter
// Module Type: analyzer
// File:        ClusterCounter_module.cc
//
// Generated at Tue Oct 25 05:26:55 2016 by Robert using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Cluster.h"

#include "TTree.h"

namespace tutorial {

class ClusterCounter : public art::EDAnalyzer {
public:
  explicit ClusterCounter(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  ClusterCounter(ClusterCounter const &) = delete;
  ClusterCounter(ClusterCounter &&) = delete;
  ClusterCounter & operator = (ClusterCounter const &) = delete;
  ClusterCounter & operator = (ClusterCounter &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.  
  void beginJob() override;
  void endJob() override;

private:
  size_t fEvNumber;

  TTree *fEventTree;
  size_t fNClusters;

  TTree *fClusterTree;
  size_t fNHits;

  // ******* fcl parameters *******
  art::InputTag fClusterModuleLabel;
  size_t fMinSize;
};

ClusterCounter::ClusterCounter(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
    fClusterModuleLabel = p.get< std::string >("ClusterModuleLabel");
    fMinSize = p.get< size_t >("MinSize");
}

void ClusterCounter::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)

    fEventTree = tfs->make<TTree>("EventTree", "event by event info");
    fEventTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fEventTree->Branch("nclusters", &fNClusters, "fNClusters/I");

    fClusterTree = tfs->make<TTree>("ClusterTree", "cluster by cluster info");
    fClusterTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fClusterTree->Branch("nhits", &fNHits, "fNHits/I");
}

void ClusterCounter::analyze(art::Event const & evt)
{
    fEvNumber = evt.id().event();
    mf::LogVerbatim("ClusterCounter") << "ClusterCounter module on event #" << fEvNumber;

    // use auto to make the line shorter when you remember all art types,
    // the type here is: art::ValidHandle< std::vector<recob::Cluster> >
    auto clusterHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);

    fNClusters = 0;

    // if you are old c++ granpa (or need cluster index, which may indeed happen!):
    // for (size_t i = 0; i < clusterHandle->size(); ++i)
    // or:
    for (auto const & clu : *clusterHandle) // loop over recob::Cluster's stored in the std::vector 
    {
        fNHits = clu.NHits();
        fClusterTree->Fill();

        if (fNHits >= fMinSize) { ++fNClusters; }
    }
    fEventTree->Fill();
}

void ClusterCounter::endJob()
{
    mf::LogVerbatim("ClusterCounter") << "ClusterCounter finished job";
}

} // tutorial namespace

DEFINE_ART_MODULE(tutorial::ClusterCounter)
