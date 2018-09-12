////////////////////////////////////////////////////////////////////////
// Class:       SaveSpacePoints
// Plugin Type: analyzer (art v2_11_02)
// File:        SaveSpacePoints_module.cc
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

#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"

namespace proto {
  class SaveSpacePoints;
}


class proto::SaveSpacePoints : public art::EDAnalyzer {
public:
  explicit SaveSpacePoints(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SaveSpacePoints(SaveSpacePoints const &) = delete;
  SaveSpacePoints(SaveSpacePoints &&) = delete;
  SaveSpacePoints & operator = (SaveSpacePoints const &) = delete;
  SaveSpacePoints & operator = (SaveSpacePoints &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  const art::InputTag fSpacePointModuleLabel;
  const art::InputTag fBeamModuleLabel;
  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;
  double evttime;

  // space point information
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;
  std::vector<double> vcharge;

  // beam information
  std::vector<double> trkStartx;
  std::vector<double> trkStarty;
  std::vector<double> trkStartz;

  std::vector<double> trkEndx;
  std::vector<double> trkEndy;
  std::vector<double> trkEndz;

};


proto::SaveSpacePoints::SaveSpacePoints(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fSpacePointModuleLabel(p.get< art::InputTag >("SpacePointModuleLabel")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel"))
{}

void proto::SaveSpacePoints::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();
  vx.clear();
  vy.clear();
  vz.clear();
  vcharge.clear();
  trkStartx.clear();
  trkStarty.clear();
  trkStartz.clear();
  trkEndx.clear();
  trkEndy.clear();
  trkEndz.clear();

  art::Handle< std::vector<recob::SpacePoint> > spsHandle;
  std::vector< art::Ptr<recob::SpacePoint> > sps;
  if (evt.getByLabel(fSpacePointModuleLabel, spsHandle))
    art::fill_ptr_vector(sps, spsHandle);

  art::Handle< std::vector<recob::PointCharge> > pcsHandle;
  std::vector< art::Ptr<recob::PointCharge> > pcs;
  if (evt.getByLabel(fSpacePointModuleLabel, pcsHandle))
    art::fill_ptr_vector(pcs, pcsHandle);

  for (size_t i = 0; i<sps.size(); ++i){
    vx.push_back(sps[i]->XYZ()[0]);
    vy.push_back(sps[i]->XYZ()[1]);
    vz.push_back(sps[i]->XYZ()[2]);
    vcharge.push_back(pcs[i]->charge());
  }

  art::Handle< std::vector<beam::ProtoDUNEBeamEvent> > pdbeamHandle;
  std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
  if (evt.getByLabel(fBeamModuleLabel, pdbeamHandle))
    art::fill_ptr_vector(beaminfo, pdbeamHandle);
  
  if (beaminfo.size()){
    auto & tracks = beaminfo[0]->GetBeamTracks();
    for (size_t i = 0; i<tracks.size(); ++i){
      //std::cout<<i<<" "<<tracks[i].NPoints()<<std::endl;
      trkStartx.push_back(tracks[i].Vertex().X());
      trkStarty.push_back(tracks[i].Vertex().Y());
      trkStartz.push_back(tracks[i].Vertex().Z());
      trkEndx.push_back(tracks[i].End().X());
      trkEndy.push_back(tracks[i].End().Y());
      trkEndz.push_back(tracks[i].End().Z());
    }
  }

  fTree->Fill();
}

void proto::SaveSpacePoints::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("spt","space point tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("vx",&vx);
  fTree->Branch("vy",&vy);
  fTree->Branch("vz",&vz);
  fTree->Branch("vcharge",&vcharge);
  fTree->Branch("trkStartx",&trkStartx);
  fTree->Branch("trkStarty",&trkStarty);
  fTree->Branch("trkStartz",&trkStartz);
  fTree->Branch("trkEndx",&trkEndx);
  fTree->Branch("trkEndy",&trkEndy);
  fTree->Branch("trkEndz",&trkEndz);
}

DEFINE_ART_MODULE(proto::SaveSpacePoints)
