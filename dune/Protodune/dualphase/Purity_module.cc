////////////////////////////////////////////////////////////////////////
// Class:       Purity
// Module Type: analyzer
// File:        Purity_module.cc
//
// Generated at Wed Aug  2 13:55:30 2017 by Andrea Scarpelli,,, using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"

#include "TTree.h"
#include "TH1F.h"
#include "TInterpreter.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>

namespace pdunedp
{
    class Purity;
}

class pdunedp::Purity : public art::EDAnalyzer {
public:

  struct Config{
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> CalWireModuleLabel{
      Name ("CalWireModuleLabel"), Comment("Calwire data product name")
    };

    fhicl::Atom<art::InputTag> HitModuleLabel{
      Name ("HitModuleLabel"), Comment("Hit data product name")
    };

    fhicl::Atom<art::InputTag> ClusterModuleLabel{
      Name ("ClusterModuleLabel"), Comment("Cluster data product name")
    };

    fhicl::Atom<art::InputTag> TrackModuleLabel{
      Name ("TrackModuleLabel"), Comment("Track data product name")
    };
  }; // end struct

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit Purity(Parameters const & config);

  Purity(Purity const &) = delete;
  Purity(Purity &&) = delete;
  Purity & operator = (Purity const &) = delete;
  Purity & operator = (Purity &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  void beginJob();
  void endJob();

private:
  art::ServiceHandle<sim::LArG4Parameters> larParameters;
  detinfo::DetectorProperties const *detProperties = nullptr;


  art::InputTag fCalWireModuleLabel;
  art::InputTag fHitModuleLabel;
  art::InputTag fClusterModuleLabel;
  art::InputTag fTrackModuleLabel;

  int Events=0; int SkipEvents=0; int BadEvents=0;
  int fRun; int fEvent;
  int fNtotHits, fNtotTracks, fNHitsTrack;

  double fTrackLength; double fChi2Ndof;

  double fElectronsToGeV;
  double fADCToElectrons;

  TTree *fTree;
};//end class

pdunedp::Purity::Purity(Parameters const & config) : EDAnalyzer(config),
  fCalWireModuleLabel(config().CalWireModuleLabel()),
  fHitModuleLabel(config().HitModuleLabel()),
  fClusterModuleLabel(config().ClusterModuleLabel()),
  fTrackModuleLabel(config().TrackModuleLabel())
{}

void pdunedp::Purity::beginJob(){
  fElectronsToGeV = 1./larParameters->GeVToElectrons();
  fADCToElectrons = 1./detProperties->ElectronsToADC();

  //Tfile Services
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("purity", "LAr purity analysis information");
  fTree->Branch("fRun", &fRun,"fRun/I");
  fTree->Branch("fEvent", &fEvent, "fEvent/I");
  fTree->Branch("fNtotHits", &fNtotHits, "fNtotHits/I");
  fTree->Branch("fNtotTracks", &fNtotTracks, "fNtotTracks/I");
  fTree->Branch("fNHitsTrack", &fNHitsTrack, "fNHitsTrack/I");
  fTree->Branch("fTrackLength", &fTrackLength, "fTrackLength/D");
  fTree->Branch("fChi2Ndof", &fChi2Ndof, "fChi2Ndof/D");
}

void pdunedp::Purity::analyze(art::Event const & e){
  fRun = e.run();
  fEvent = e.id().event();   Events++;

  //Require data product handles
  auto CalwireHandle = e.getValidHandle<std::vector<recob::Wire> >(fCalWireModuleLabel);
  auto HitHandle = e.getValidHandle<std::vector<recob::Hit> >(fHitModuleLabel);
  auto ClusterHandle = e.getValidHandle<std::vector<recob::Cluster> >(fClusterModuleLabel);
  auto TrackHandle = e.getValidHandle<std::vector<recob::Track> >(fTrackModuleLabel);

  //check data products exists
  if(!CalwireHandle || !HitHandle || !ClusterHandle || !TrackHandle){
    BadEvents++;
    return;
  }

  //first selection of data based on ADC counts (separate energetic showers)
  float SumWireADC=0;
  for(auto const& Calwire : *CalwireHandle){
    for(float ADC : Calwire.Signal()){
      SumWireADC += ADC;
    }
  }
  float Edep = SumWireADC*fADCToElectrons*fElectronsToGeV;
  mf::LogVerbatim("pdunedp::Purity") << "Non cal energy deposit: " << Edep;
  if(Edep > 5.0){
    SkipEvents++;
    return;
  }

  fNtotHits = HitHandle->size();
  fNtotTracks = TrackHandle->size();

  art::FindManyP< recob::Hit > HitsFromTrack(TrackHandle, e, fTrackModuleLabel);
  for(size_t t=0; t<(size_t)fNtotTracks; t++){
    //caracteristics of tracks
    auto track = TrackHandle->at(t);
    fTrackLength = track.Length();
    fChi2Ndof = track.Chi2PerNdof();

    //hits track
    fNHitsTrack = HitsFromTrack.at(t).size();
  }

  //Select only interesting events()
  fTree->Fill();
}

void pdunedp::Purity::endJob(){
  mf::LogVerbatim("pdunedp::Purity") << "Total Events: " << Events;
  mf::LogVerbatim("pdunedp::Purity") << "Bad Events: " << BadEvents;
  mf::LogVerbatim("pdunedp::Purity") << "Skipped Events: " << SkipEvents;
}

DEFINE_ART_MODULE(pdunedp::Purity)
