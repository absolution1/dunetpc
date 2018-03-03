///////////////////////////////////////////////////////////////////////
// Class:       DEdxdp
// Module Type: analyzer
// File:        DEdxdp_module.cc
//
// This is a simple module that makes plot of dE/dx using 3D reconstructed tracks.
// It cab be run before each MC production on the sample of generated muon tracks.
// DP calibrations, quality control, prompt measure of effective gain
//
// Author: andrea.scarpelli@cern.ch
//
////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/Utilities/DatabaseUtil.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/ArtDataHelper/MVAReader.h"

#include "TH1.h"
#include "TTree.h"

namespace pdunedp
{
  class DEdxdp;
}

class pdunedp::DEdxdp : public art::EDAnalyzer {
public:
  explicit DEdxdp(fhicl::ParameterSet const & p);

  //plugins should not be copied or assigned.
  DEdxdp(DEdxdp const &) = delete;
  DEdxdp(DEdxdp &&) = delete;
  DEdxdp & operator = (DEdxdp const &) = delete;
  DEdxdp & operator = (DEdxdp &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  void beginJob() override;
  void endJob() override;

  void reconfigure(fhicl::ParameterSet const& p) ;

private:

  // Declare member data here.

  // void CountdEdx(const std::vector < art::Ptr< recob::Hit > > & hits, const std::vector < recob::TrackHitMeta const* > & data); // MeV/cm
  void CountdEdx(const std::vector < art::Ptr< recob::Hit > > & hits, const std::vector< recob::TrackHitMeta const* > & data);
  void ResetVars();

  TTree *fTree;
  TTree *fTreere;
  int fRun;
  int fSubRun;
  int fEvent;
  int fCryo;
  int fTPC;
  int fView;
  int fTrackID;
  size_t fNumberOfTracks;

  double fdQdx;
  double fdEdx;
  double fdQ;
  double fdx;

  // Module labels to get data products
  std::string fHitModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  calo::CalorimetryAlg fCalorimetryAlg;
  int fGainPerView;

};

pdunedp::DEdxdp::DEdxdp(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
  // More initializers here.
{
  reconfigure(p);
}

void pdunedp::DEdxdp::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTreere = tfs->make<TTree>("entries","entries tree");
  fTree = tfs->make<TTree>("DEdx","dedx tree");

  fTreere->Branch("fRun", &fRun, "fRun/I");
  fTreere->Branch("fEvent", &fEvent, "fEvent/I");
  fTreere->Branch("fNumberOfTracks", &fNumberOfTracks, "fNumberOfTracks/I");

  fTree->Branch("fRun", &fRun, "fRun/I");
  fTree->Branch("fSubRun", &fSubRun, "fSubRun/I");
  fTree->Branch("fEvent", &fEvent, "fEvent/I");
  fTree->Branch("fCryo", &fCryo, "fCryo/I");
  fTree->Branch("fTPC", &fTPC, "fTPC/I");
  fTree->Branch("fTrackID", &fTrackID, "fTrackID/I");
  fTree->Branch("fView", &fView, "fView/I");
  fTree->Branch("fdQdx", &fdQdx, "fdQdx/D");
  fTree->Branch("fdEdx", &fdEdx, "fdEdx/D");
  fTree->Branch("fdQ", &fdQ, "fdQ/D");
  fTree->Branch("fdx", &fdx, "fdx/D");
}

void pdunedp::DEdxdp::endJob()
{
}

void pdunedp::DEdxdp::reconfigure(fhicl::ParameterSet const & p)
{
  fHitModuleLabel = p.get< std::string >("HitModuleLabel");
  fTrackModuleLabel = p.get< std::string >("TrackModuleLabel");
  fGainPerView      = p.get<int>("GainPerView");
}

void pdunedp::DEdxdp::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  ResetVars();
  fRun = e.run();
  fSubRun = e.subRun();
  fEvent = e.id().event();

  // tracks
  auto trkHandle = e.getValidHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  fNumberOfTracks = trkHandle->size();

  const art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trkHandle, e, fTrackModuleLabel);

  if (fmthm.isValid())
  {
    // loop over tracks
    for (size_t t = 0; t < trkHandle->size(); ++t)
    {
	     auto vhit = fmthm.at(t);
	     auto vmeta = fmthm.data(t);

       fTrackID = t;
	     CountdEdx(vhit, vmeta);
    }
  }
  fTreere->Fill();
}

void pdunedp::DEdxdp::CountdEdx(const std::vector < art::Ptr< recob::Hit > > & hits, const std::vector< recob::TrackHitMeta const* > & data) // MeV/cm
{
  for (size_t h = 0; h < hits.size(); ++h)
  {
        fCryo = hits[h]->WireID().Cryostat;
        fTPC =  hits[h]->WireID().TPC;
        fView = hits[h]->WireID().Plane;

	unsigned short plane = hits[h]->WireID().Plane;
      	unsigned short time = hits[h]->PeakTime();

      	double dqadc = hits[h]->Integral();
      	if (!std::isnormal(dqadc) || (dqadc < 0)) dqadc = 0.0;

      	double t0 = 0;

        fdQdx = 0.0;

        fdQ = dqadc;
      	fdx = data[h]->Dx();
      	if ((fdx > 0) && (fdQ > 0))
        {
        	 fdQdx = fdQ/fdx;
        	 fdQdx /= fGainPerView;
         	 fdEdx = fCalorimetryAlg.dEdx_AREA(fdQdx, time, plane, t0);
         	 fdQdx *= fCalorimetryAlg.LifetimeCorrection(time, t0);
         	 //if (fdEdx > 35) fdEdx = 35;
        }
        fTree->Fill();
    }
}

void pdunedp::DEdxdp::ResetVars()
{
  fdQdx = 0.0;
  fdEdx = 0.0;
  fdQ = 0.0;
  fdx = 0.0;
  fNumberOfTracks = 0;
}


DEFINE_ART_MODULE(pdunedp::DEdxdp)
