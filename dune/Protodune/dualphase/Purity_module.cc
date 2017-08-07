
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
//TODO: Stitch to multiple tracks (more than 2)
//TODO: generalization to multiple TPCs
//TODO: Add Code for the purity
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
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "dune/DetSim/Service/DPhaseSimChannelExtractService.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <stdio.h>
#include <math.h>       /* acos */


namespace pdunedp
{
    struct bHitInfo;
    class  Purity;
}

struct pdunedp::bHitInfo
{
	bHitInfo(size_t i, double x, double e, int w) :
		Index(i), dE(e), dx(x), wire(w)
	{ }
	size_t Index;
	double dE, dx;
	int wire;
};

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

    fhicl::Atom<double> Length{
      Name ("Length"), Comment("minimal length to define a mip")
    };

    fhicl::Atom<double>StitchAngle{
      Name("StitchAngle"), Comment("Max. stitching angle")
    };

    fhicl::Atom<double>StitchDistance{
      Name("StitchDistance"), Comment("Max. stitching distance")
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
  double GetCharge(std::vector<recob::Hit> hits);
  double GetCharge(const std::vector<art::Ptr<recob::Hit> >  hits);
  bool IsCrossing(TVector3 Start, TVector3 End);
  bool StitchTracks(recob::Track Track1, recob::Track Track2, TVector3 & Edge1, TVector3 & Edge2);
  void Make_dEdx(std::vector< double > & dEdx, std::vector< double > & range,
                            const std::vector< pdunedp::bHitInfo > & hits, double rmax, int Flip);
  void FillTajectoryGraph(std::map<size_t, recob::Track > MipCandidate,
                                                             art::FindManyP<recob::Hit>  HitsTrk);
  void FindMipInfo(recob::Track mip, std::vector<art::Ptr<recob::Hit>> vhits,
    const std::vector<const recob::TrackHitMeta*, std::allocator<const recob::TrackHitMeta*> > vmeta);

private:
  calo::CalorimetryAlg fCalorimetryAlg;
  art::ServiceHandle<geo::Geometry> geom;

  art::InputTag fCalWireModuleLabel;
  art::InputTag fHitModuleLabel;
  art::InputTag fClusterModuleLabel;
  art::InputTag fTrackModuleLabel;

  int Events=0; int SkipEvents=0; int BadEvents=0;
  int fRun; int fEvent;
  int fNtotHits, fNtotTracks, fNHitsTrack;
  int fTrkNum; int fNumOfMips;
  int fMipNum; int fMipIndex;
  int fPlane; int fIndex;

  double fTrackLength; double fChi2Ndof;
  double fHitsCharge; //in fC
  double fTrackCharge; //in fC
  double fHitsTrackCharge; //in fC
  double fMipCharge, fMipLength, fMipPhi, fMipTheta;
  double fCosX; double fCosY; double fCosZ;
  double fADC; double fCorrectedChargeX; double fCharge;
  double fDrift;
  double fdEdx; double fRange;

  //double fElectronsToGeV;
  double fADCToElectrons;
  double fElectronCharge = 1.60217662e-19/1e-15; //in fC
  double fMipChargeCm = 10; //in fC/cm

  double fTotalGain;
  double fECut;
  double fLength;
  double fStitchAngle; double fStitchDistance;
  std::vector<double> fVolCut;

  std::map< size_t, std::vector< pdunedp::bHitInfo >[3] > fTrk2InfoMap; // hits info sorted by views
  std::map<size_t, recob::Track > TrackList;
  std::vector<double> goodevents;

  TTree *fTree; TTree *fTreeTrk; TTree *fTreeMip; TTree *fTreePurity; TTree *fTreeCalib;
  TH2D *hTrkTrajectory_0; TH2D *hTrkTrajectory_1;
};//end class

pdunedp::Purity::Purity(Parameters const & config) : EDAnalyzer(config),
  fCalorimetryAlg(config().CalorimetryAlg()),
  fCalWireModuleLabel(config().CalWireModuleLabel()),
  fHitModuleLabel(config().HitModuleLabel()),
  fClusterModuleLabel(config().ClusterModuleLabel()),
  fTrackModuleLabel(config().TrackModuleLabel()),
  fTotalGain(config().TotalGain()),
  fECut(config().EnergyCut()),
  fLength(config().Length()),
  fStitchAngle(config().StitchAngle()),
  fStitchDistance(config().StitchDistance()),
  fVolCut(config().VolCut())
{}

void pdunedp::Purity::beginJob(){
  auto const* dp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fADCToElectrons = 1./dp->ElectronsToADC();
  //auto simChannelExtract = &*art::ServiceHandle<detsim::DPhaseSimChannelExtractService>();
  //fTotalGain = simChannelExtract->GainPerView()*2;

  //Tfile Services
  art::ServiceHandle<art::TFileService> tfs;

  //one entry for every event
  fTree = tfs->make<TTree>("purity", "LAr purity analysis information");
  fTree->Branch("fRun", &fRun,"fRun/I");
  fTree->Branch("fEvent", &fEvent, "fEvent/I");
  fTree->Branch("fNtotHits", &fNtotHits, "fNtotHits/I");
  fTree->Branch("fNtotTracks", &fNtotTracks, "fNtotTracks/I");
  fTree->Branch("fNumOfMips", &fNumOfMips, "fNumOfMips/I");
  fTree->Branch("fHitsCharge", &fHitsCharge, "fHitsCharge/D");
  fTree->Branch("fHitsTrackCharge", &fHitsTrackCharge, "fHitsTrackCharge/D");

  //one entry for every track
  fTreeTrk = tfs->make<TTree>("TrkInfo", "Information on tracks");
  fTreeTrk->Branch("fRun", &fRun,"fRun/I");
  fTreeTrk->Branch("fEvent", &fEvent, "fEvent/I");
  fTreeTrk->Branch("fTrkNum", &fTrkNum, "fTrkNum/I");
  fTreeTrk->Branch("fChi2Ndof", &fChi2Ndof, "fChi2Ndof/D");
  fTreeTrk->Branch("fTrackLength", &fTrackLength, "fTrackLength/D");
  fTreeTrk->Branch("fTrackCharge", &fTrackCharge, "fTrackCharge/D");
  fTreeTrk->Branch("fNHitsTrack", &fNHitsTrack, "fNHitsTrack/I");

  //one entry for every selected mip (stitched mips are considered separately)
  fTreeMip = tfs->make<TTree>("MipInfo", "Information on mips");
  fTreeMip->Branch("fRun", &fRun,"fRun/I");
  fTreeMip->Branch("fEvent", &fEvent, "fEvent/I");
  fTreeMip->Branch("fMipIndex", &fMipIndex, "fMipIndex/I");
  fTreeMip->Branch("fMipNum", &fMipNum, "fMipNum/I");
  fTreeMip->Branch("fMipCharge", &fMipCharge, "fMipcharge/D");
  fTreeMip->Branch("fMipLength", &fMipLength, "fMipLength/D");
  fTreeMip->Branch("fMipLength", &fMipLength, "fMipLength/D");
  fTreeMip->Branch("fMipPhi", &fMipPhi, "fMipPhi/D");
  fTreeMip->Branch("fMipTheta", &fMipTheta, "fMipTheta/D");

  //branch for purity measuremtns usign mips: one entry per mip
  fTreePurity =tfs->make<TTree>("Purity","Corrected charge for purity analysis");
  fTreePurity->Branch("fRun", &fRun,"fRun/I");
  fTreePurity->Branch("fEvent", &fEvent, "fEvent/I");
  fTreePurity->Branch("fMipIndex", &fMipIndex, "fMipIndex/I");
  fTreePurity->Branch("fMipNum", &fMipNum, "fMipNum/I");
  fTreePurity->Branch("fPlane", &fPlane, "fPlane/I");
  fTreePurity->Branch("fCosX", &fCosX, "fCosX/D");
  fTreePurity->Branch("fCosY", &fCosY, "fCosY/D");
  fTreePurity->Branch("fCosZ", &fCosZ, "fCosZ/D");
  fTreePurity->Branch("fADC", &fADC, "fADC/D");
  fTreePurity->Branch("fCharge", &fCharge, "fCharge/D");
  fTreePurity->Branch("fCorrectedChargeX", &fCorrectedChargeX, "fCorrectedChargeX/D");
  fTreePurity->Branch("fDrift", &fDrift, "fDrift/D");

  fTreeCalib = tfs->make<TTree>("calibration", "dE/dx info");
  fTreeCalib->Branch("fRun", &fRun, "fRun/I");
  fTreeCalib->Branch("fEvent", &fEvent, "fEvent/I");
  fTreeCalib->Branch("fIndex", &fIndex, "fIndex/I");
  fTreeCalib->Branch("fdEdx", &fdEdx, "fdEdx/D");
  fTreeCalib->Branch("fRange", &fRange, "fRange/D");
  fTreeCalib->Branch("fPlane", &fPlane, "fPlane/I");


  //TODO<<--Generalize these histograms to different geometries
  hTrkTrajectory_0 = tfs->make<TH2D>("hTrkTrajectory_0", "Selected mips hit position view 0;Channel;Ticks", 320, 0, 319, 1667, 0, 1666);
  hTrkTrajectory_1 = tfs->make<TH2D>("hTrkTrajectory_1", "Selected mips hit position view 1;Channel;Ticks", 960, 0, 959, 1667, 0, 1666);
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
    mf::LogVerbatim("pdunedp::Purity") << "Too much energy energy deposited! Not a mip";
    SkipEvents++;
    return;
  }

  //Get Hits charge
  fNtotHits = HitHandle->size();
  fHitsCharge = GetCharge(*HitHandle);

  fNtotTracks = (int)TrackHandle->size();
  if(fNtotTracks == 0){
    mf::LogVerbatim("pdunedp::Purity") << "Event has no 3D tracks";
    SkipEvents++;
    return;
  }

  art::FindManyP< recob::Hit > HitsFromTrack(TrackHandle, e, fTrackModuleLabel);
  art::FindManyP< recob::Hit, recob::TrackHitMeta > TrackHitMeta(TrackHandle, e, fTrackModuleLabel);

  double ChargeTrk=0;
  //int skippedTrk=0;
  for(size_t t=0; t<(size_t)fNtotTracks; t++){
    //retrive information for every track (not mip selected yet)
    auto track = TrackHandle->at(t);
    TrackList[t] = track;
    fTrackLength = track.Length();
    fChi2Ndof = track.Chi2PerNdof();
    fNHitsTrack = HitsFromTrack.at(t).size();
    fTrkNum = (int)t;
    fTrackCharge = GetCharge(HitsFromTrack.at(t));
    ChargeTrk+= GetCharge(HitsFromTrack.at(t));
    fTreeTrk->Fill();
  }
  fHitsTrackCharge = ChargeTrk; //<-- Sum charge from all event

  std::map<size_t, recob::Track > fMipCandidate;
  /*Try to select mips. Particle passing trougth the detector (Selecting a fiducial
  volume to keep into account the inefficiencies of the lems ).
  Tracks that doesn't meet the criteria of crossing the detector after a first pass
  are stitched togheter and the crossing criteria are evaluated again.
  Stitch is done only on two consecuitve tracks. TODO: Stitch more than two tracks*/

  std::map<size_t, recob::Track >::iterator It;
  for(It = TrackList.begin(); It !=TrackList.end(); It++){
    //first check if the track is already inside the candidate vector
    if (fMipCandidate.find(It->first) != fMipCandidate.end()){
      continue;
    };

    TVector3 Start = (It->second).Vertex(); TVector3 End = (It->second).End();
    if(!IsCrossing(Start, End) && (TrackList.size() > 1)){
      //It is not crossing, let's see if it can be stitched to something else
      for(std::map<size_t, recob::Track >::iterator It2=TrackList.begin(); It2!=TrackList.end(); It2++){
        if(It->first == It2->first){ continue; }
        TVector3 newStart; TVector3 newEnd;
        if(!StitchTracks(It->second, It2->second, newStart, newEnd)){
          continue; //Not stitching: continue searcing
        }else{
          if(IsCrossing(newStart, newEnd)){
            fMipCandidate[It->first] = It->second;
            fMipCandidate[It2->first] = It2->second;
          }else{
            break; //Stop the loop: we stitched and still not crossing...
          }
        }
      }//end while
    }else if(!IsCrossing(Start, End) && (TrackList.size() == 1)){
      continue;
    }else{
      fMipCandidate[It->first] = It->second;
    }
  }

  if(!fMipCandidate.size()){
    mf::LogVerbatim("pdunedp::Purity") << "No mips selected! ";
    SkipEvents++;
    return;
  }
  //print the list of mips:
  mf::LogVerbatim("pdunedp::Purity") << "List of mips candidates in event " << fEvent;
  size_t num=0;
  for(std::map<size_t, recob::Track >::iterator mip = fMipCandidate.begin(); mip !=fMipCandidate.end(); mip++){
    mf::LogVerbatim("pdunedp::Purity") << "Track: " << mip->first;
    fMipNum = (int)num++;
    fMipIndex = mip->first;
    fMipCharge = GetCharge(HitsFromTrack.at(mip->first));
    fMipLength = (mip->second).Length();
    fMipPhi = (mip->second).Phi();      //polar angle
    fMipTheta = (mip->second).Theta();  //Azimutal angle
    auto vhit = TrackHitMeta.at(mip->first);
    auto vmeta = TrackHitMeta.data(mip->first);
    FindMipInfo(mip->second, vhit, vmeta);

    fTreeMip->Fill();
  }

  for (auto const & trkEntry : fTrk2InfoMap)
  {
    //check flip
    int Flip = 0; // tag track if it is reconstructed in the opposite direction
    if (fMipCandidate[trkEntry.first].End().Z() < fMipCandidate[trkEntry.first].Vertex().Z()){ Flip = 1; }

    for(unsigned int plane=0; plane < geom->Nplanes(0); plane++){ //TODO<<--Extend to different geometries
      std::vector< double > dEdx, range;
      auto const & info = trkEntry.second;
      Make_dEdx(dEdx, range, info[plane], fMipCandidate[trkEntry.first].Length(), Flip);
      for (size_t i = 0; i < dEdx.size(); ++i){
        fdEdx = dEdx[i];
        fPlane = plane;
        fIndex = trkEntry.first;
        fRange = range[i];
        fTreeCalib->Fill();
      }
    }
  }

  FillTajectoryGraph(fMipCandidate, HitsFromTrack);



  goodevents.push_back(fEvent);
  fTree->Fill();
}//end analyzer

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

bool pdunedp::Purity::StitchTracks(recob::Track Track1, recob::Track Track2, TVector3 & Edge1, TVector3 & Edge2){
  //Attempt to sticth two tracks togheter. return the non stitched verteces of the two tracks.
  //  mf::LogVerbatim("pdunedp::Purity") << "Doing Stitch";
  bool Stitch = false;

  double Dx = sqrt(pow(Track1.End().X()-Track2.Vertex().X(),2) + pow(Track1.End().Y()-Track2.Vertex().Y(),2)+pow(Track1.End().Z()-Track2.Vertex().Z(),2));
  double Angle= (180.0/3.14159)*Track1.EndDirection().Angle(Track2.VertexDirection());
  int Criteria = 1;

  if(Dx > sqrt(pow(Track1.End().X()-Track2.End().X(),2) + pow(Track1.End().Y()-Track2.End().Y(),2) + pow(Track1.End().Z()-Track2.End().Z(),2)))
	{
	    Dx= sqrt(pow(Track1.End().X()-Track2.End().X(),2) + pow(Track1.End().Y()-Track2.End().Y(),2) + pow(Track1.End().Z()-Track2.End().Z(),2));
	    Angle = 180. - (180.0/3.14159)*Track1.EndDirection().Angle(Track2.EndDirection());
	    Criteria= 2;
	 }

  if(Dx > sqrt(pow(Track1.Vertex().X()-Track2.End().X(),2) + pow(Track1.Vertex().Y()-Track2.End().Y(),2) + pow(Track1.Vertex().Z()-Track2.End().Z(),2)))
  {
	    Dx = sqrt(pow(Track1.Vertex().X()-Track2.End().X(),2) + pow(Track1.Vertex().Y()-Track2.End().Y(),2) + pow(Track1.Vertex().Z()-Track2.End().Z(),2));
	    Angle = (180.0/3.14159)*Track1.VertexDirection().Angle(Track2.EndDirection());
	    Criteria = 3;
	 }

  if(Dx > sqrt(pow(Track1.Vertex().X()-Track2.Vertex().X(),2) + pow(Track1.Vertex().Y()-Track2.Vertex().Y(),2) + pow(Track1.Vertex().Z()-Track2.Vertex().Z(),2)))
	{
	    Dx = sqrt(pow(Track1.Vertex().X()-Track2.Vertex().X(),2) + pow(Track1.Vertex().Y()-Track2.Vertex().Y(),2) + pow(Track1.Vertex().Z()-Track2.Vertex().Z(),2));
	    Angle = 180. - (180.0/3.14159)*Track1.VertexDirection().Angle(Track2.VertexDirection());
	    Criteria = 4;
	}

  //check the stitching criteria
  if( (fStitchAngle > Angle) && (fStitchDistance > Dx)){
    Stitch = true;
    switch (Criteria) {
      case 1:
        Edge1 = Track1.Vertex(); Edge2 = Track2.End();
        break;
      case 2:
        Edge1 = Track1.Vertex(); Edge2 = Track2.Vertex();
        break;
      case 3:
        Edge1 = Track1.End(); Edge2 = Track2.Vertex();
        break;
      case 4:
        Edge1 = Track1.End(); Edge2 = Track2.End();
        break;
      default:
        mf::LogError("pdunedp::Purity") << "Unknown error!";
        break;
    }
  }
  return Stitch;
}

bool pdunedp::Purity::IsCrossing(TVector3 Start, TVector3 End){
  bool isCrossing = false;
  //define the boundaries (for the moment consider only 1 TPC) //<--TODO Generalisaztion to Multiple TPCs
  double vtx[3] = {0, 0, 150}; //<<--TODO Generalisaztion
  art::ServiceHandle<geo::Geometry> geom;
  geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);
  if (geom->HasTPC(idtpc)) // <----Assuming there is only one TPC
  {
    const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
    double minx = tpcgeo.MinX() + fVolCut[0]; double maxx = tpcgeo.MaxX() - fVolCut[1];
    double miny = tpcgeo.MinY() + fVolCut[2]; double maxy = tpcgeo.MaxY() - fVolCut[3];
    double minz = tpcgeo.MinZ() + fVolCut[4]; double maxz = tpcgeo.MaxZ() - fVolCut[5];

  //check if both start and end are both outside the fiducial volume (in at least one direction )
    if(   (   ((minx > Start.X()) || (maxx < Start.X()))
           || ((miny > Start.Y()) || (maxy < Start.Y()))
           || ((minz > Start.Z()) || (maxz < Start.Z())) ) //condiotion on start
        &&(   ((minx > End.X()  ) || (maxx < End.X()  ))
           || ((miny > End.Y()  ) || (maxy < End.Y()  ))
           || ((minz > End.Z()  ) || (maxz < End.Z()  )) ) //condiotion on end
      ){
        isCrossing=true;
      }else{
        isCrossing=false;
      }
  //from the point of view of X only particles that goes trougth the whole detector can be
  //considered particles on trigger (check all the possible cases)
    TVector3 Diff = Start-End;
    double ProjX = fabs(Diff.Unit().X());
    if( acos(1./ProjX) < 30*(TMath::Pi()/180) ){ //<<-- only on very vertical tracks
      if((maxx < Start.X())||(maxx < End.X())){
        if((minx < Start.X())||(minx < End.X()))
          isCrossing=false;
      }
      if((minx > Start.X())||(minx > End.X())){
        if((maxx > Start.X())||(maxx > End.X()))
          isCrossing=false;
      }
    }
  //Mips must be above certain length (for not be confused with random hadronic processes)
    if( (Diff).Mag() < fLength){
      isCrossing=false;
    }
  }
  return isCrossing;
}

void pdunedp::Purity::FillTajectoryGraph(std::map<size_t, recob::Track > MipCandidate,
                                           art::FindManyP<recob::Hit> HitsTrk){
  //Merge the 2D hit position view in a single graph to evaluate defects in mip selection
  //<<--TODO Generic geometry
  std::map<size_t, recob::Track >::iterator It;
  for(It = MipCandidate.begin(); It !=MipCandidate.end(); It++){
    auto hits = HitsTrk.at(It->first);
    if(!hits.size()){ continue; }
    for(auto const  hit : hits){
      unsigned short plane = hit->WireID().Plane;
      switch (plane) {
        case 0:
          hTrkTrajectory_0->Fill(hit->WireID().Wire, hit->PeakTime());
          break;
        case 1:
          hTrkTrajectory_1->Fill(hit->WireID().Wire, hit->PeakTime());
          break;
        default:
          mf::LogError("pdunedp::Purity") << "Invalid Plane;";
          break;
      }
    }
  }
  return;
}

void pdunedp::Purity::FindMipInfo(recob::Track mip, std::vector<art::Ptr<recob::Hit>> vhits,
  const std::vector<const recob::TrackHitMeta*, std::allocator<const recob::TrackHitMeta*> > vmeta){
  /*This function is intended to caluclate the most important quantities from a mip
  and fill a tree for further analysis. There will be an entry for every hit in the mip.
  Stitched mips are considered independently. T0 is assumed to be 0, as consequence of
  the selection done*/

  if(!vhits.size()){return;}
  for(size_t h =0; h< vhits.size(); h++){
    unsigned int plane= vhits[h]->WireID().Plane;
    int wire = vhits[h]->WireID().Wire;
    //double pitch = geom->WirePitch(vhits[h]->View());
    double dQadc = vhits[h]->Integral();
    if (!std::isnormal(dQadc) || (dQadc < 0)) continue;
    size_t idx = vmeta[h]->Index();
    double dx =  vmeta[h]->Dx();
    mf::LogVerbatim("pdunedp::Purity") << "dx: " << dx << " idx " << idx;

    double dQ = dQadc*fCalorimetryAlg.ElectronsFromADCArea(dQadc, plane)*fElectronCharge;
    fTrk2InfoMap[fMipIndex][plane].emplace_back(idx, dx, dQ, wire);

    TVector3 Dir = mip.Vertex() - mip.End();
    //Direction cosines of the track
    double CosX = 1./Dir.Unit().X();
    double CosY = 1./Dir.Unit().Y();
    double CosZ = 1./Dir.Unit().Z();

    double dQ_corrX = dQ*fabs(1/CosX) ;//charge corrected over the effective pitch
    double PeakTime = vhits[h]->PeakTime();

    //Fill Purity tree
    fPlane  = plane;
    fCosX   = CosX;
    fCosY   = CosY;
    fCosZ   = CosZ;
    fADC    = dQadc;
    fCharge = dQ;
    fCorrectedChargeX = dQ_corrX;
    fDrift = PeakTime;
    //fHitX = //drift coordinate of my hit expressed in
    fTreePurity->Fill();
  }
  return;
}

void pdunedp::Purity::Make_dEdx(std::vector< double > & dEdx, std::vector< double > & range,
                                        const std::vector< pdunedp::bHitInfo > & hits, double rmax, int Flip){
  if (!hits.size()) return;

	dEdx.clear(); range.clear();

	int i0 = hits.size() - 1; int i1 = -1; int di = -1;
	if (Flip) {i0 = 0; i1 = hits.size(); di = 1;}

	double de = 0.0;
	double dx = 0.0;
	double r0 = 0.0; double r1 = 0.0; double r = 0.0;

	double minDx = 0.1; // can be a parameter

	while ((i0 != i1) && (r < rmax))
	{
		dx = 0.0; de = 0.0;
		while ((i0 != i1) && (dx <= minDx))
		{
			de += hits[i0].dE;
			dx += hits[i0].dx;
			i0 += di;
		}

		r0 = r1;
		r1 += dx;
		r = 0.5 * (r0 + r1);

		if ((de > 0.0) && (dx > 0.0) && (r < rmax))
		{
			dEdx.push_back(de/dx);
			range.push_back(r);
		}
	}
  return;
}

void pdunedp::Purity::Clear(){
  TrackList.clear();
  fTrk2InfoMap.clear();
}

void pdunedp::Purity::endJob(){
  mf::LogVerbatim("pdunedp::Purity") << "====== Run " << fRun << " summary ======";
  mf::LogVerbatim("pdunedp::Purity") << "Total Events: " << Events;
  mf::LogVerbatim("pdunedp::Purity") << "Bad Events: " << BadEvents;
  mf::LogVerbatim("pdunedp::Purity") << "Skipped Events: " << SkipEvents;
  mf::LogVerbatim("pdunedp::Purity") << "selected Events: " << Events-BadEvents-SkipEvents <<" \nList: ";
  for(int goodevent : goodevents){
      mf::LogVerbatim("pdunedp::Purity") << ".. " << goodevent ; //<----write this list on file (?)
  }
  mf::LogVerbatim("pdunedp::Purity") << "===========================";

}

DEFINE_ART_MODULE(pdunedp::Purity)
