


////////////////////////////////////////////////////////////////////////
// Class:       DetProptest
// Plugin Type: analyzer (art v2_11_02)
// File:        DetProptest_module.cc
//
// Generated at Wed Jul 11 09:18:07 2018 by Tingjun Yang using cetskelgen
// from cetlib version v3_03_01.
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

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"

#include "dune/DuneObj/ProtoDUNEBeamEvent.h"


#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"

namespace proto {
  class DetProptest;
}


class proto::DetProptest : public art::EDAnalyzer {
public:
  explicit DetProptest(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DetProptest(DetProptest const &) = delete;
  DetProptest(DetProptest &&) = delete;
  DetProptest & operator = (DetProptest const &) = delete;
  DetProptest & operator = (DetProptest &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  //void beginJob() override;

private:
  detinfo::DetectorProperties const* fDetProp;
  double tickToDist;


};


proto::DetProptest::DetProptest(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
  
{}

void proto::DetProptest::analyze(art::Event const & evt)
{
    fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  std::cout<<fDetProp->Efield(0)<<std::endl;
  std::cout<<fDetProp->Efield(1)<<std::endl;
  std::cout<<fDetProp->Efield(2)<<std::endl;
  std::cout<<fDetProp->Efield(3)<<std::endl;
  //std::cout<<fDetProp->Efield(3)<<std::endl;


}

// void proto::DetProptest::beginJob()
// {
//   art::ServiceHandle<art::TFileService> tfs;
//   fTree = tfs->make<TTree>("spt","space point tree");
//   fTree->Branch("run",&run,"run/I");
//   fTree->Branch("subrun",&subrun,"subrun/I");
//   fTree->Branch("event",&event,"event/I");
//   fTree->Branch("trigger",&trigger,"trigger/I");
//   fTree->Branch("evttime",&evttime,"evttime/D");
//   fTree->Branch("vx",&vx);
//   fTree->Branch("vy",&vy);
//   fTree->Branch("vz",&vz);
//   fTree->Branch("vcharge",&vcharge);
//   fTree->Branch("vtrackid",&vtrackid);
//   fTree->Branch("beamPosx",&beamPosx);
//   fTree->Branch("beamPosy",&beamPosy);
//   fTree->Branch("beamPosz",&beamPosz);
//   fTree->Branch("beamDirx",&beamDirx);
//   fTree->Branch("beamDiry",&beamDiry);
//   fTree->Branch("beamDirz",&beamDirz);
//   fTree->Branch("beamMomentum",&beamMomentum);
//   fTree->Branch("tof", &tof, "tof/D");
//   fTree->Branch("ckov0status", &ckov0status, "ckov0status/S");
//   fTree->Branch("ckov1status", &ckov1status, "ckov1status/S");
// }

DEFINE_ART_MODULE(proto::DetProptest)
