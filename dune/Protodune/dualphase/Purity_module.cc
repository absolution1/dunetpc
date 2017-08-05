////////////////////////////////////////////////////////////////////////
// Class:       Purity
// Module Type: analyzer
// File:        Purity_module.cc
//
// Generated at Wed Aug  2 13:55:30 2017 by Andrea Scarpelli
//(andrea.scarpelli@cern.ch) using artmod
// from cetpkgsupport v1_11_00.
// Analyzer for purity measurements in protodunedp (and 3x1x1 prototype)
//TODO: GainView to be configured from service
//TODO: Add energy of tracks
//TODO: Add Code for the purity
//TODO: Add track stitching
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
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
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

    fhicl::Table<calo::CalorimetryAlg::Config> CalorimetryAlg {
        Name("CalorimetryAlg"), Comment("Used to calculate electrons from ADC area.")
    };

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

    fhicl::Sequence<double>VolCut{
      Name("VolCut"), Comment("Volume Cut to select a going trougth muon")
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
  void StitchTracks(std::map<size_t, recob::Track > TrackList,
        std::vector<std::map<size_t, recob::Track > > & MipCandidate);
  bool IsMip(recob::Track track, std::map<size_t, recob::Track > & TrackList,
          size_t t, const std::vector<art::Ptr<recob::Hit> >  HitsTrk  , double ChargeDep);
  double GetCharge(std::vector<recob::Hit> hits);
  double GetCharge(const std::vector<art::Ptr<recob::Hit> >  hits);


private:
  calo::CalorimetryAlg fCalorimetryAlg;

  art::InputTag fCalWireModuleLabel;
  art::InputTag fHitModuleLabel;
  art::InputTag fClusterModuleLabel;
  art::InputTag fTrackModuleLabel;

  int Events=0; int SkipEvents=0; int BadEvents=0;
  int fRun; int fEvent;
  int fNtotHits, fNtotTracks, fNHitsTrack;
  int fTrkNum; int fNumOfMips;

  double fTrackLength; double fChi2Ndof;
  double fHitsCharge; //in fC
  double fTrackCharge; //in fC
  double fHitsTrackCharge; //in fC

  //double fElectronsToGeV;
  double fADCToElectrons;
  double fElectronCharge = 1.60217662e-19/1e-15; //in fC
  double fMipChargeCm = 10; //in fC/cm

  double fTotalGain;
  double fECut;
  double fDriftGap;
  double fLength;
  std::vector<double> fVolCut;

  std::map<size_t, recob::Track > fTrackList;
  std::vector<double> goodevents;

  TTree *fTree; TTree *fTreeTrk;
};//end class

pdunedp::Purity::Purity(Parameters const & config) : EDAnalyzer(config),
  fCalorimetryAlg(config().CalorimetryAlg()),
  fCalWireModuleLabel(config().CalWireModuleLabel()),
  fHitModuleLabel(config().HitModuleLabel()),
  fClusterModuleLabel(config().ClusterModuleLabel()),
  fTrackModuleLabel(config().TrackModuleLabel()),
  fTotalGain(config().TotalGain()),
  fECut(config().EnergyCut()),
  fDriftGap(config().DriftGap()),
  fLength(config().Length()),
  fVolCut(config().VolCut())
{}

void pdunedp::Purity::beginJob(){
  auto const* dp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fADCToElectrons = 1./dp->ElectronsToADC();
  //auto simChannelExtract = &*art::ServiceHandle<detsim::DPhaseSimChannelExtractService>();
  //fTotalGain = simChannelExtract->GainPerView()*2;

  //Tfile Services
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("purity", "LAr purity analysis information");
  fTree->Branch("fRun", &fRun,"fRun/I");
  fTree->Branch("fEvent", &fEvent, "fEvent/I");
  fTree->Branch("fNtotHits", &fNtotHits, "fNtotHits/I");
  fTree->Branch("fNtotTracks", &fNtotTracks, "fNtotTracks/I");
  fTree->Branch("fNumOfMips", &fNumOfMips, "fNumOfMips/I");
  fTree->Branch("fNHitsTrack", &fNHitsTrack, "fNHitsTrack/I");
  fTree->Branch("fTrackLength", &fTrackLength, "fTrackLength/D");
  fTree->Branch("fChi2Ndof", &fChi2Ndof, "fChi2Ndof/D");
  fTree->Branch("fHitsCharge", &fHitsCharge, "fHitsCharge/D");
  fTree->Branch("fHitsTrackCharge", &fHitsTrackCharge, "fHitsTrackCharge/D");

  fTreeTrk = tfs->make<TTree>("TrkInfo", "Information on tracks");
  fTreeTrk->Branch("fRun", &fRun,"fRun/I");
  fTreeTrk->Branch("fEvent", &fEvent, "fEvent/I");
  fTreeTrk->Branch("fTrkNum", &fTrkNum, "fTrkNum/I");
  fTreeTrk->Branch("fTrackCharge", &fTrackCharge, "fTrackCharge/D");
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
  mf::LogVerbatim("pdunedp::Purity")<< "Energy: " << Edep << " " << Efrac;
  if(Efrac > fECut){
    SkipEvents++;
    return;
  }

  //Get Hits charge
  fNtotHits = HitHandle->size();
  fHitsCharge = GetCharge(*HitHandle);

  fNtotTracks = (int)TrackHandle->size();
  if(fNtotTracks == 0){
    SkipEvents++;
    return;
  }

  art::FindManyP< recob::Hit > HitsFromTrack(TrackHandle, e, fTrackModuleLabel);
  double ChargeTrk=0;
  for(size_t t=0; t<(size_t)fNtotTracks; t++){
    //retrive information for every track (not mip selected yet)
    auto track = TrackHandle->at(t);
    fTrackList[t] = track;
    fTrackLength = track.Length();
    fChi2Ndof = track.Chi2PerNdof();
    fNHitsTrack = HitsFromTrack.at(t).size();
    fTrkNum = (int)t;
    fTrackCharge = GetCharge(HitsFromTrack.at(t));
    ChargeTrk+= GetCharge(HitsFromTrack.at(t));
    fTreeTrk->Fill();
  }
    fHitsTrackCharge = ChargeTrk; //<-- Sum charge from all event

  std::vector<std::map<size_t, recob::Track > >  fMipCandidate;
  StitchTracks(fTrackList, fMipCandidate);

  //is mip



  //purity stuff goes here



  goodevents.push_back(fEvent);
  fTree->Fill();
}

double pdunedp::Purity::GetCharge(const std::vector<art::Ptr<recob::Hit> >  hits){
  //It returns the uncalibrated charge in the detector given a list of hits (summed on both views)
  if(!hits.size()){ return 0.0;}

  double charge=0;
  for(auto const  hit : hits){
    unsigned short plane = hit->WireID().Plane;
    double dqadc = hit->Integral();
    if (!std::isnormal(dqadc) || (dqadc < 0)) continue;
    charge += dqadc*fCalorimetryAlg.ElectronsFromADCArea(dqadc, plane)*fElectronCharge;
  }
    return charge;
}

double pdunedp::Purity::GetCharge(std::vector<recob::Hit> hits){
  //It returns the uncalibrated charge in the detector given a list of hits (summed on both views)
  if(!hits.size()){ return 0.0;}

  double charge=0;
  for(auto hit : hits){
    unsigned short plane = hit.WireID().Plane;
    double dqadc = hit.Integral();
    if (!std::isnormal(dqadc) || (dqadc < 0)) continue;
    charge += dqadc*fCalorimetryAlg.ElectronsFromADCArea(dqadc, plane)*fElectronCharge;
  }
    return charge;
}

void pdunedp::Purity::StitchTracks(std::map<size_t, recob::Track > TrackList,
    std::vector<std::map<size_t, recob::Track > > & MipCandidate){
  //Attempt to stitch tracks together
  if(TrackList.size()==0){ return; }
  //std::vector<size_t> StitchList;



  return;
}

bool pdunedp::Purity::IsMip(recob::Track track, std::map<size_t, recob::Track > & TrackList,
              size_t t, const std::vector<art::Ptr<recob::Hit> >  HitsTrk  , double ChargeDep){
  bool isMip = false;
  art::ServiceHandle<geo::Geometry> geom;

  //check start and end position of the track
  double StartPos[3] = {track.Start().X(), track.Start().Y(), track.Start().Z()};
  double EndPos[3] = {track.End().X(), track.End().Y(), track.End().Z()};


//All the tracks starting close to the anode and with a length above a certain safe length can be used for the purity analysis
//Single tracks going trougth the detector can be used as well (check fraction of energy on that track)
  if( (StartPos[0] > fDriftGap) && (track.Length() > fLength) ){
    TrackList[t] = track;
    return isMip = true;
  }
  else if( (EndPos[0] > fDriftGap) && (track.Length() > fLength)){
    TrackList[t] = track;
    return isMip = true;
  }else if( (GetCharge(HitsTrk)/ChargeDep > 0.30) && (track.Length() > fLength)){
    TrackList[t] = track;
    return isMip = true;
  /*
    geo::TPCID idtpc = geom->FindTPCAtPosition(StartPos);

    if (geom->HasTPC(idtpc)) // <----Assuming there is only one TPC
	  {
		  const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
		  double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
		  double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
		  double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

		  if( (fabs(fabs(EndPos[0]  - StartPos[0]) - fabs(maxx -minx)) > 2*fVolCut[0])
        || (fabs(fabs(EndPos[1] - StartPos[1]) - fabs(maxy -miny)) > 2*fVolCut[1])
        || (fabs(fabs(EndPos[2] - StartPos[2]) - fabs(maxz -minz)) > 2*fVolCut[2])){
          TrackList[t] = track;
          return isMip = true;
      }else{
        return isMip = false;
      }
	  }
  */
  }else{
    return isMip = false;
  }
  return isMip;
}

void pdunedp::Purity::Clear(){
  fTrackList.clear();
}

void pdunedp::Purity::endJob(){
  mf::LogVerbatim("pdunedp::Purity") << "Total Events: " << Events;
  mf::LogVerbatim("pdunedp::Purity") << "Bad Events: " << BadEvents;
  mf::LogVerbatim("pdunedp::Purity") << "Skipped Events: " << SkipEvents;
  mf::LogVerbatim("pdunedp::Purity") << "selected Event List: ";
  for(int goodevent : goodevents){
      mf::LogVerbatim("pdunedp::Purity") << ".. " << goodevent ; //<----write this list on file (?)
  }
}

DEFINE_ART_MODULE(pdunedp::Purity)
