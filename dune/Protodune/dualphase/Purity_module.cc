////////////////////////////////////////////////////////////////////////
// Class:       Purity
// Module Type: analyzer
// File:        Purity_module.cc
//
// Generated at Wed Aug  2 13:55:30 2017 by Andrea Scarpelli andrea.scarpelli@cern.ch using artmod
// from cetpkgsupport v1_11_00.
// Analyzer for purity measurements in protodunedp (and 3x1x1 prototype)
//TODO: GainView to be configured from service  
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
#include "dune/DetSim/Service/DPhaseSimChannelExtractService.h"

#include "TTree.h"
#include "TH1F.h"
#include "TInterpreter.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
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

    fhicl::Atom<double> TotalGain{
      Name ("TotalGain"), Comment("Total Gain of the detector") //<---TODO insert it form service
    };

    fhicl::Atom<double> EnergyCut{
      Name ("EnergyCut"), Comment("Cut over the event energy")
    };

    fhicl::Atom<double> DriftGap{
      Name ("DriftGap"), Comment("Gap the mip track start can be found in")
    };

    fhicl::Atom<double> Length{
      Name ("Length"), Comment("minimal length to define a mip")
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
  void Clear();

private:
  //art::ServiceHandle<sim::LArG4Parameters> larParameters;
  detinfo::DetectorProperties const *detProperties = nullptr;

  art::InputTag fCalWireModuleLabel;
  art::InputTag fHitModuleLabel;
  art::InputTag fClusterModuleLabel;
  art::InputTag fTrackModuleLabel;

  int Events=0; int SkipEvents=0; int BadEvents=0;
  int fRun; int fEvent;
  int fNtotHits, fNtotTracks, fNHitsTrack;

  double fTrackLength; double fChi2Ndof;

  //double fElectronsToGeV;
  double fADCToElectrons;
  double fElectronCharge = 1.60217662e-19/1e-15; //in fC
  double fMipChargeCm = 10; //in fC/cm

  double fTotalGain;
  double fECut;
  double fDriftGap;
  double fLength;

  std::map<size_t, recob::Track > fTrackList;

  TTree *fTree;
};//end class

pdunedp::Purity::Purity(Parameters const & config) : EDAnalyzer(config),
  fCalWireModuleLabel(config().CalWireModuleLabel()),
  fHitModuleLabel(config().HitModuleLabel()),
  fClusterModuleLabel(config().ClusterModuleLabel()),
  fTrackModuleLabel(config().TrackModuleLabel()),
  fTotalGain(config().TotalGain()),
  fECut(config().EnergyCut()),
  fDriftGap(config().DriftGap()),
  fLength(config().Length())
{}

void pdunedp::Purity::beginJob(){
  //fElectronsToGeV = 1./larParameters->GeVToElectrons();
  fADCToElectrons = 1./detProperties->ElectronsToADC();
  //auto simChannelExtract = &*art::ServiceHandle<detsim::DPhaseSimChannelExtractService>();
  //fTotalGain = simChannelExtract->GainPerView()*2;

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

  Clear();

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
  float Edep = SumWireADC*fADCToElectrons*fElectronCharge;
  float Efrac = Edep/(fMipChargeCm*350*fTotalGain);
  if(Efrac > fECut){
    SkipEvents++;
    return;
  }

  fNtotHits = HitHandle->size();
  fNtotTracks = TrackHandle->size();

  art::FindManyP< recob::Hit > HitsFromTrack(TrackHandle, e, fTrackModuleLabel);
  int skippedTracks=0;
  for(size_t t=0; t<(size_t)fNtotTracks; t++){
    //caracteristics of tracks
    auto track = TrackHandle->at(t);
    fTrackLength = track.Length();
    fChi2Ndof = track.Chi2PerNdof();

//selecting mips
    if( track.Start().X() > fDriftGap ){ fTrackList[t]= track; }
    else if( track.End().X() > fDriftGap ){ fTrackList[t]= track; }
    else if( (TrackHandle->size() ==1) && (track.Length() < fLength) ){ fTrackList[t]= track; }
    else{ skippedTracks++;}

    //hits track
    fNHitsTrack = HitsFromTrack.at(t).size();
  }//end loop tracks

  if(skippedTracks - fNtotTracks){
    SkipEvents++;
    return;
  }

  //do more analysis here

  //Select only interesting events()
  fTree->Fill();
}

void pdunedp::Purity::Clear(){
  fTrackList.clear();
}

void pdunedp::Purity::endJob(){
  mf::LogVerbatim("pdunedp::Purity") << "Total Events: " << Events;
  mf::LogVerbatim("pdunedp::Purity") << "Bad Events: " << BadEvents;
  mf::LogVerbatim("pdunedp::Purity") << "Skipped Events: " << SkipEvents;
}

DEFINE_ART_MODULE(pdunedp::Purity)
