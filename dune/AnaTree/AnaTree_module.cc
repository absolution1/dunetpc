////////////////////////////////////////////////////////////////////////
// Class:       AnaTree
// Module Type: analyzer
// File:        AnaTree_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
//
//  ** modified by Muhammad Elnimr to access track information and clusters as well
//  mmelnimr@as.ua.edu
//  August 2014
//  ** Numerous changes since Feb 2015, Karl Warburton k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////
// Framework includes
//#include "art/Framework/Core/EDProducer.h"
//#include "art/Framework/Core/FindManyP.h"
//#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"
//#include "art/Utilities/InputTag.h"
//#include <iterator>
//#include "SimpleTypesAndConstants/PhysicalConstants.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/OpFlash.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/RawData/ExternalTrigger.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/AnalysisBase/Calorimetry.h"
#include "lardata/AnalysisBase/T0.h"
#include "lardata/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackUtils.h" // lar::TrackPitchInView()

#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TLorentzVector.h"
#include "TH2F.h"
#include "TFile.h"

//standard library includes
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits> // std::numeric_limits<>

const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxHits       = 10000; //maximum number of hits
const int kMaxClust      = 10000; //maximum number of clusters
const int kMaxTrackHits  = 1000;  //maximum number of space points
const int kMaxFlash      = 1000;  //maximum number of flashes
const int kMaxTrig       = 1000;  //maximum number of triggers




namespace AnaTree {
  class AnaTree;
}
namespace {

  // Local functions.
  
}
class AnaTree::AnaTree : public art::EDAnalyzer {
public:
  explicit AnaTree(fhicl::ParameterSet const & p);
  virtual ~AnaTree();

  void analyze(art::Event const & e) override;

  void beginJob() override;
  void endJob();
  //void reconfigure(fhicl::ParameterSet const & p) override;

private:

  void ResetVars();
  
  TTree* fTree;

  // Run information
  int run;
  int subrun;
  int event;
  double evttime;
  float efield[3];
  int t0;

  int ntracks_reco;         //number of reconstructed tracks
  float trkstartx[kMaxTrack];
  float trkstarty[kMaxTrack];
  float trkstartz[kMaxTrack];
  float trkendx[kMaxTrack];
  float trkendy[kMaxTrack];
  float trkendz[kMaxTrack];

  int nMCParticles;
  int nTPCHits_MC[kMaxTrack];
  float StartTime_MC[kMaxTrack];
  float trkstartx_MC[kMaxTrack];
  float trkstarty_MC[kMaxTrack];
  float trkstartz_MC[kMaxTrack];
  float trkendx_MC[kMaxTrack];
  float trkendy_MC[kMaxTrack];
  float trkendz_MC[kMaxTrack];
  float TPCLen_MC[kMaxTrack];
  float EnergyDeposited_MC[kMaxTrack];
  float Hits_posx_MC[kMaxTrack][10000];
  float Hits_posy_MC[kMaxTrack][10000];
  float Hits_posz_MC[kMaxTrack][10000];
  float Hits_inTPC_MC[kMaxTrack][10000];
  float Hits_mom_MC[kMaxTrack][10000];
  float Hits_momx_MC[kMaxTrack][10000];
  float Hits_momy_MC[kMaxTrack][10000];
  float Hits_momz_MC[kMaxTrack][10000];
  float Hits_E_MC[kMaxTrack][10000];
  float Hits_T_MC[kMaxTrack][10000];
  float Hits_dEdx_MC[kMaxTrack][10000];
  float Hits_dE_MC[kMaxTrack][10000];
  float Hits_dx_MC[kMaxTrack][10000];
 
  int    trkid_MC[kMaxTrack];
  int    trkpdg_MC[kMaxTrack];
  int    trkndaughters_MC[kMaxTrack];
  int    StartInTPC_MC[kMaxTrack];
  int    EndInTPC_MC[kMaxTrack];
  float trkmom_MC[kMaxTrack];
  float trkmom_XMC[kMaxTrack];
  float trkmom_YMC[kMaxTrack];
  float trkmom_ZMC[kMaxTrack];
  float trkenergy_MC[kMaxTrack];
  float trkstartdoc_XMC[kMaxTrack];
  float trkstartdoc_YMC[kMaxTrack];
  float trkstartdoc_ZMC[kMaxTrack];
  int    trkMother_MC[kMaxTrack];
  int    trkNumDaughters_MC[kMaxTrack];
  int    trkFirstDaughter_MC[kMaxTrack];
  int    trkPrimary_MC[kMaxTrack];
  float trktheta_xz_MC[kMaxTrack];
  float trktheta_yz_MC[kMaxTrack];
  float trketa_xy_MC[kMaxTrack];
  float trketa_zy_MC[kMaxTrack];
  float trktheta_MC[kMaxTrack];
  float trkphi_MC[kMaxTrack];

  float trkd2[kMaxTrack];
  float trklen[kMaxTrack];
  float trklen_L[kMaxTrack];
  int trkid[kMaxTrack];
  
  int    trkMCTruthTrackID[kMaxTrack];
  float trkMCTruthT0[kMaxTrack];
  int    trkPhotonCounterID[kMaxTrack];
  float trkPhotonCounterT0[kMaxTrack];
  float trkPhotonCounterConf[kMaxTrack];

  float trktheta_xz[kMaxTrack];
  float trketa_xy[kMaxTrack];
  float trktheta_yz[kMaxTrack];
  float trketa_zy[kMaxTrack];
  float trktheta[kMaxTrack];
  float trkphi[kMaxTrack];
  float trkdQdxSum[kMaxTrack];
  float trkdQdxAverage[kMaxTrack];
  float trkdEdxSum[kMaxTrack];
  float trkdEdxAverage[kMaxTrack];

  float trkdedx2[kMaxTrack][3][1000];
  float trkdqdx[kMaxTrack][3][1000];
  float trkpitchHit[kMaxTrack][3][1000];
  float trkkinE[kMaxTrack][3];
  float trkrange[kMaxTrack][3];
  float trkTPC[kMaxTrack][3][1000];
  float trkplaneid[kMaxTrack][3][1000];
  float trkresrg[kMaxTrack][3][1000];
  float trkPosx[kMaxTrack][3][1000];
  float trkPosy[kMaxTrack][3][1000];
  float trkPosz[kMaxTrack][3][1000]; 
  float mcang_x[kMaxTrack];
  float mcang_y[kMaxTrack];
  float mcang_z[kMaxTrack];
  float mcpos_x[kMaxTrack];
  float mcpos_y[kMaxTrack];
  float mcpos_z[kMaxTrack];

  float trkstartdcosx[kMaxTrack];
  float trkstartdcosy[kMaxTrack];
  float trkstartdcosz[kMaxTrack];
  float trkenddcosx[kMaxTrack];
  float trkenddcosy[kMaxTrack];
  float trkenddcosz[kMaxTrack];
  int    ntrkhits[kMaxTrack];
  float trkx[kMaxTrack][kMaxTrackHits];
  float trky[kMaxTrack][kMaxTrackHits];
  float trkz[kMaxTrack][kMaxTrackHits];
  
  float trkpitch[kMaxTrack][3];
  int    nhits;
  int    nhits2;
  int    nclust;

  int    hit_tpc[kMaxHits];
  int    hit_plane[kMaxHits];
  int    hit_wire[kMaxHits];
  int    hit_channel[kMaxHits];
  float hit_peakT[kMaxHits];
  float hit_charge[kMaxHits];
  float hit_ph[kMaxHits];
  int    hit_trkid[kMaxHits];

  int    flash_total;
  float flash_time[kMaxFlash];
  float flash_width[kMaxFlash];
  float flash_abstime[kMaxFlash];
  float flash_YCentre[kMaxFlash];
  float flash_YWidth[kMaxFlash];
  float flash_ZCentre[kMaxFlash];
  float flash_ZWidth[kMaxFlash];
  float flash_TotalPE[kMaxFlash];

  int ntrigs;
  int trig_time[kMaxTrig];
  int trig_id[kMaxTrig];

  std::string fTrigModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fClusterModuleLabel;
  std::string fTrkSpptAssocModuleLabel;
  std::string fHitSpptAssocModuleLabel;
  std::string fSimulationProducerLabel; 
  std::string fCalorimetryModuleLabel; 
  std::string fMCTruthT0ModuleLabel;
  std::string fPhotonT0ModuleLabel;
  std::string fFlashModuleLabel;

  float fElectronsToGeV; // conversion factor
  art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

  //  std::map<int, MCHists> fMCHistMap;       // Indexed by pdg id.
  //  std::map<int, RecoHists> fRecoHistMap;   // Indexed by pdg id.

};

AnaTree::AnaTree::AnaTree(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fTrigModuleLabel         ( pset.get< std::string >("TrigModuleLabel"))
  , fHitsModuleLabel         ( pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel        ( pset.get< std::string >("TrackModuleLabel"))
  , fClusterModuleLabel      ( pset.get< std::string >("ClusterModuleLabel"))
  , fTrkSpptAssocModuleLabel ( pset.get< std::string >("TrkSpptAssocModuleLabel"))
  , fHitSpptAssocModuleLabel ( pset.get< std::string >("HitSpptAssocModuleLabel"))
  , fSimulationProducerLabel ( pset.get< std::string >("SimulationLabel"))
  , fCalorimetryModuleLabel  ( pset.get< std::string >("CalorimetryModuleLabel"))
  , fMCTruthT0ModuleLabel    ( pset.get< std::string >("MCTruthT0ModuleLabel"))
  , fPhotonT0ModuleLabel     ( pset.get< std::string >("PhotonT0ModuleLabel"))
  , fFlashModuleLabel        ( pset.get< std::string >("FlashModuleLabel"))
{

}

AnaTree::AnaTree::~AnaTree()
{
  // Clean up dynamic memory and other resources here.

}

void AnaTree::AnaTree::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  ResetVars();

  art::ServiceHandle<geo::Geometry> geom;
  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const *timeservice = lar::providerFrom<detinfo::DetectorClocksService>();
  //fClock = timeservice->TPCClock();
  art::ServiceHandle<cheat::BackTracker> bktrk;

  
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();

  efield[0] = detprop->Efield(0);
  efield[1] = detprop->Efield(1);
  efield[2] = detprop->Efield(2);
  
  t0 = detprop->TriggerOffset();
  
  art::Handle< std::vector<raw::ExternalTrigger> > trigListHandle;
  std::vector<art::Ptr<raw::ExternalTrigger> > triglist;
  if (evt.getByLabel(fTrigModuleLabel,trigListHandle))
    art::fill_ptr_vector(triglist, trigListHandle);

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::Handle< std::vector<recob::OpFlash> > flashListHandle;
  std::vector<art::Ptr<recob::OpFlash> > flashlist;
  if (evt.getByLabel(fFlashModuleLabel, flashListHandle))
    art::fill_ptr_vector(flashlist, flashListHandle);

  art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
  evt.getByLabel(fSimulationProducerLabel, simChannelHandle);

  const sim::ParticleList& plist = bktrk->ParticleList();
  simb::MCParticle *particle=0;

  //  Event.getByLabel(fSimulationProducerLabel, particleHandle);
  //  art::Handle< std::vector<simb::MCParticle> > particleHandle;

  art::Handle< std::vector<recob::Track> > trackh;
  evt.getByLabel(fTrackModuleLabel, trackh);

  art::Handle< std::vector< art::PtrVector < recob::Track > > > trackvh;
  evt.getByLabel(fTrackModuleLabel, trackvh);
  
  // Protect against invalid art::Handle (both here and when trackvh is used below)
  // TODO: How to do this for art::PtrVector without the uninitialised iterator?
  std::vector< art::PtrVector<recob::Track> >::const_iterator cti; 
  if (trackvh.isValid()) cti = trackvh->begin(); 

  //track information
  ntracks_reco=tracklist.size();

  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  
  // Get Cryostat information.....
  int c=0;//only one cryo   
  double boundaries[6];
  double TempBounds[6];
  int NCryo = geom->Ncryostats();
  for ( int b=0; b<6; b++) boundaries[b] = 0;
  for ( int c1 = 0; c1 < NCryo; c1++ ) {
    geom->CryostatBoundaries(TempBounds, c); // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )
    if ( boundaries[0] > TempBounds [0] ) boundaries[0] = TempBounds [0];
    if ( boundaries[1] < TempBounds [1] ) boundaries[1] = TempBounds [1];
    if ( boundaries[2] > TempBounds [2] ) boundaries[2] = TempBounds [2];
    if ( boundaries[3] < TempBounds [3] ) boundaries[3] = TempBounds [3];
    if ( boundaries[4] > TempBounds [4] ) boundaries[4] = TempBounds [4];
    if ( boundaries[5] < TempBounds [5] ) boundaries[5] = TempBounds [5];
  }
  //std::cout << boundaries[0] << " " << boundaries[1] << " " << boundaries[2] << " " << boundaries[3] << " " <<boundaries[4] << " " << boundaries[5] << std::endl;


  // ------------------------------------
  //  EXTERNAL TRIGGER STUFF 
  // ------------------------------------

  ntrigs = std::min(int(triglist.size()),kMaxTrig);

  for (int i = 0; i < ntrigs; ++i){
    trig_time[i] = triglist[i]->GetTrigTime();
    trig_id[i] = triglist[i]->GetTrigID();
  }
 
  // ------------------------------------
  //  END EXTERNAL TRIGGER STUFF 
  // ------------------------------------

  // ------------------------------------
  //  NOW DO THE TRACK LIST STUFF 
  // ------------------------------------
  
  if ( trackListHandle.isValid() ) {
    art::FindManyP<recob::SpacePoint> fmsp  (trackListHandle, evt, fTrackModuleLabel);
    art::FindManyP<recob::Hit>        fmht  (trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal (trackListHandle, evt, fCalorimetryModuleLabel);
    art::FindMany<anab::T0>           fmt0  (trackListHandle, evt, fMCTruthT0ModuleLabel);
    art::FindMany<anab::T0>           fmphot(trackListHandle, evt, fPhotonT0ModuleLabel);
    for(int i=0; i<std::min(int(tracklist.size()),kMaxTrack);++i){

      //***************************************************
      //   T0 stuff - So can correct X Start/End positions
      //***************************************************
      double TickT0 = 0;
      if ( fmt0.isValid() ) {
	std::vector<const anab::T0*> T0s = fmt0.at(i);
	//std::cout << int(T0s.size()) << std::endl;
	for (size_t t0size =0; t0size < T0s.size(); t0size++) {
	  //std::cout << "Using Trigger Type Monte Carlo Truth!!! T0 " << T0s[t0size]->Time() << ", GEANT4 track id " << T0s[t0size]->TriggerBits() << ", Size " << T0s[t0size]->ID() << std::endl;
	  trkMCTruthT0[i]  = T0s[t0size]->Time();
	  TickT0 = trkMCTruthT0[i] / detprop->SamplingRate();
	  trkMCTruthTrackID[i] = T0s[t0size]->TriggerBits();
	} // T0 size
      } // T0 valid
      if ( fmphot.isValid() ) {
        std::vector<const anab::T0*> PhotT0 = fmphot.at(i);
        for (size_t T0it=0; T0it<PhotT0.size(); ++T0it) {
          trkPhotonCounterT0[i]   = PhotT0[T0it]->Time();
          trkPhotonCounterID[i]   = PhotT0[T0it]->TriggerBits();
	  trkPhotonCounterConf[i] = PhotT0[T0it]->TriggerConfidence();
        }
      }       
      // Add other T0 handles...with same structure...
      //**************
      // END T0 stuff
      //***************

      //******************************
      // Hit Level Stuff - Correct X
      //******************************
      std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(i);
      double Hit_Size = allHits.size();
      //*********************************
      // End Hit Level Stuff - Correct X
      //*********************************

      trkid[i]=i;
      trackStart.clear();
      trackEnd.clear();
      memset(larStart, 0, 3);
      memset(larEnd, 0, 3);
      tracklist[i]->Extent(trackStart,trackEnd);
      tracklist[i]->Direction(larStart,larEnd);
      trkstartx[i]      = trackStart[0] - detprop->ConvertTicksToX( TickT0, allHits[Hit_Size-1]->WireID().Plane, allHits[Hit_Size-1]->WireID().TPC, allHits[Hit_Size-1]->WireID().Cryostat ); // Correct X, last entry is first 'hit'
      trkstarty[i]      = trackStart[1];
      trkstartz[i]      = trackStart[2];
      trkendx[i]        = trackEnd[0] - detprop->ConvertTicksToX( TickT0, allHits[0]->WireID().Plane, allHits[0]->WireID().TPC, allHits[0]->WireID().Cryostat ); // Correct X, first entry is last 'hit'
      trkendy[i]        = trackEnd[1];
      trkendz[i]        = trackEnd[2];
      trkstartdcosx[i]  = larStart[0];
      trkstartdcosy[i]  = larStart[1];
      trkstartdcosz[i]  = larStart[2];
      trkenddcosx[i]    = larEnd[0];
      trkenddcosy[i]    = larEnd[1];
      trkenddcosz[i]    = larEnd[2];
      TLorentzVector v1(trackStart[0],trackStart[1],trackStart[2],0);
      TLorentzVector v2(trackEnd[0],trackEnd[1],trackEnd[2],0);
      trklen[i]=(v2-v1).Rho();

      art::Ptr<recob::Track> ptrack(trackh, i);
      const recob::Track& track = *ptrack;
      trklen_L[i]=track.Length();

      //double recotime = 0.;
      //double trackdx = recotime * 1.e-3 * detprop->DriftVelocity();  // cm
      // Fill histograms involving reco tracks only.
      int ntraj = track.NumberTrajectoryPoints();
      if(ntraj > 0) {
	TVector3 dir = track.VertexDirection();
	trktheta_xz[i] = std::atan2(dir.X(), dir.Z());
	trktheta_yz[i] = std::atan2(dir.Y(), dir.Z());
	trketa_xy[i] = std::atan2(dir.X(), dir.Y());
	trketa_zy[i] = std::atan2(dir.Z(), dir.Y());
	trktheta[i]=dir.Theta();
	trkphi[i]=dir.Phi();
      }

      double distance_squared=0;
      TVector3 V1(trackStart[0],trackStart[1],trackStart[2]);
      TVector3 V2(trackEnd[0],trackEnd[1],trackEnd[2]);
      TVector3 vOrth=(V2-V1).Orthogonal();
      TVector3 pointVector=V1;

      // *********************
      //  Space Point stuff:
      // *********************
      if(fmsp.isValid() ){
	ntrkhits[i] = fmsp.at(i).size();
	//std::cout << "Space Points " <<  ntrkhits[i] << std::endl;
	//	double distance_squared=0;
	double distance=0;

	std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
	for (size_t j = 0; j<spts.size(); ++j){
	  TVector3 sptVector(spts[j]->XYZ()[0],spts[j]->XYZ()[1],spts[j]->XYZ()[2]);
	  TVector3 vToPoint=sptVector-pointVector;
	  distance=(vOrth.Dot(vToPoint))/vOrth.Mag();
	  distance_squared+=distance *distance;
	  trkx[i][j] = spts[j]->XYZ()[0];
	  trky[i][j] = spts[j]->XYZ()[1];
	  trkz[i][j] = spts[j]->XYZ()[2];

	}
	distance_squared=distance_squared/spts.size();
	trkd2[i]=distance_squared;
      }
      // ***********************
      // END Space Point stuff:
      // ***********************

      // *********************
      //  Calorimetric stuff:
      // *********************
      if (fmcal.isValid()){
	std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
	//std::cout<<"calo size "<<calos.size()<<std::endl;
	for (size_t jj = 0; jj<calos.size(); ++jj){
	  trkrange[i][jj] = calos[jj]->Range();
	  trkkinE[i][jj]  = (calos[jj] -> KineticEnergy());
	  //std::cout << trkkinE[i][jj] << std::endl;
	  int tt= calos[jj] -> dEdx().size();
	  for(int k = 0; k < tt; ++k) {
	    trkdedx2[i][jj][k]  = (calos[jj] -> dEdx())[k];
	    trkdqdx[i][jj][k]   = (calos[jj] -> dQdx())[k];
	    trkpitchHit[i][jj][k]  = (calos[jj] -> TrkPitchVec())[k];
	    trkresrg[i][jj][k]  = (calos[jj] -> ResidualRange())[k];
	    trkTPC[i][jj][k]    =(calos[jj]->PlaneID()).TPC;
	    trkplaneid[i][jj][k]=(calos[jj]->PlaneID()).Plane;

	    trkPosx[i][jj][k]    = (calos[jj]->XYZ())[k].x();
	    trkPosy[i][jj][k]    = (calos[jj]->XYZ())[k].y();
	    trkPosz[i][jj][k]    = (calos[jj]->XYZ())[k].z();

	    trkdQdxSum[i] += trkdqdx[i][jj][k];
	    trkdEdxSum[i] += trkdedx2[i][jj][k];
	  }
	}
	trkdQdxAverage[i] = trkdQdxSum[i] / calos.size();
	trkdEdxAverage[i] = trkdEdxSum[i] / calos.size();
      }
      // ************************
      //  End Calorimetric stuff
      // ************************
    
      // ---------FIND THE TRACK PITCH IN EACH PLANE------------
      for (int j = 0; j<3; ++j){
	try {
	  if (j==0)
	    trkpitch[i][j] = lar::TrackPitchInView(*(tracklist[i]), geo::kU);
	  else if (j==1)
	    trkpitch[i][j] = lar::TrackPitchInView(*(tracklist[i]), geo::kV);
	  else if (j==2)
	    trkpitch[i][j] = lar::TrackPitchInView(*(tracklist[i]), geo::kZ);
	}
	catch( cet::exception &e) {
	  mf::LogWarning("AnaTree")<<"caught exception "<<e<<"\n setting pitch to 0";
	  trkpitch[i][j] = 0;
	}
      }
    
    } // TRACK LOOP
  }
  // -----------------------------------------------
  // NOW DO THE CLUSTER/HIT STUFF.......
  // -----------------------------------------------
  if ( hitListHandle.isValid() ) {
    art::FindMany<recob::Track>       fmtk (hitListHandle  , evt, fTrackModuleLabel);
    nhits  = hitlist.size();
    nhits2 = std::min(int(hitlist.size()),kMaxHits);
    nclust = clusterlist.size();
    for (int i = 0; i < nhits2; ++i){
      unsigned int channel = hitlist[i]->Channel();
      geo::WireID wireid = hitlist[i]->WireID();
      hit_tpc[i]     = wireid.TPC;
      hit_plane[i]   = wireid.Plane;
      hit_wire[i]    = wireid.Wire;
      hit_channel[i] = channel;
      hit_peakT[i]   = hitlist[i]->PeakTime();
      hit_charge[i]  = hitlist[i]->Integral();
      hit_ph[i]      = hitlist[i]->PeakAmplitude();
      if (fmtk.at(i).size()!=0){
      	hit_trkid[i] = fmtk.at(i)[0]->ID();
      }
      if ( i == kMaxHits ) break;
    }
  }
  // -----------------------------------------------
  // NOW DO THE FLASH STUFF.......
  // -----------------------------------------------
  if ( flashListHandle.isValid() ) {
    flash_total = flashlist.size();
    std::cout << "Total Number of flashes for this event...." << flash_total << std::endl;
    for ( int f = 0; f < std::min(flash_total,kMaxHits); ++f ) {
      flash_time[f]      = flashlist[f]->Time();
      flash_width[f]     = flashlist[f]->TimeWidth();
      flash_abstime[f]   = flashlist[f]->AbsTime();
      flash_YCentre[f]   = flashlist[f]->YCenter();
      flash_YWidth[f]    = flashlist[f]->YWidth();
      flash_ZCentre[f]   = flashlist[f]->ZCenter();
      flash_ZWidth[f]    = flashlist[f]->ZWidth();
      flash_TotalPE[f]   = flashlist[f]->TotalPE();
    }
  }

  // -----------------------------------------------
  // NOW DO ALL THE MONTE CARLO TRUTH STUFF.......
  // -----------------------------------------------
 
  int i=0; // particle index
  for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    particle = ipar->second;

    // Line below selects only primaries, have removed so get all particles.
    // So want to add if Process=="Primary" in macro, this way can see particles other than primaries too.
    //if(!(particle->Process()=="primary" && abs(particle->PdgCode())== abs(fPdg))) continue;
    
    if (particle->Process() == "primary" ) trkPrimary_MC[i] = 1;
    else trkPrimary_MC[i] = 0;
    trkid_MC[i]=particle->TrackId();
    trkpdg_MC[i]=particle->PdgCode();
    trkndaughters_MC[i]=particle->NumberDaughters();
    StartTime_MC[i] = particle->T();                                 // nsec
    trkMother_MC[i]=particle->Mother();
    trkNumDaughters_MC[i]=particle->NumberDaughters();
    trkFirstDaughter_MC[i]=particle->FirstDaughter();
          
    double xyztArray[4];
    int zz =0;
    int zz2 =0;
    bool insideActiveVolume    = false;
    int nTPCHitsCounts         = 0;
    double XPlanePosition      = 0;
    double DriftTimeCorrection = 0;
    double TimeAtPlane         = 0;
    double XDriftVelocity      = detprop->DriftVelocity()*1e-3; //cm/ns
    double WindowSize          = detprop->NumberTimeSamples() * timeservice->TPCClock().TickPeriod() * 1e3;
  
    int numberTrajectoryPoints = particle->NumberTrajectoryPoints(); // Looking at each MC hit
    for(int ii=0;ii<numberTrajectoryPoints;++ii) {
      const TLorentzVector& tmpPosition=particle->Position(ii);
      double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
      Hits_posx_MC[i][ii] = tmpPosition[0];
      Hits_posy_MC[i][ii] = tmpPosition[1];
      Hits_posz_MC[i][ii] = tmpPosition[2];
      Hits_mom_MC[i][ii]  = particle->P(ii);
      Hits_momx_MC[i][ii] = particle->Px(ii);
      Hits_momy_MC[i][ii] = particle->Py(ii);
      Hits_momz_MC[i][ii] = particle->Pz(ii);
      Hits_E_MC[i][ii]    = particle->E(ii);
      Hits_T_MC[i][ii]    = particle->T(ii);

      if (ii != 0) { // Work out MCTruth dEdx
	Hits_dx_MC[i][ii] = pow ( (Hits_posx_MC[i][ii-1]-Hits_posx_MC[i][ii])*(Hits_posx_MC[i][ii-1]-Hits_posx_MC[i][ii])
				  + (Hits_posy_MC[i][ii-1]-Hits_posy_MC[i][ii])*(Hits_posy_MC[i][ii-1]-Hits_posy_MC[i][ii])
				  + (Hits_posz_MC[i][ii-1]-Hits_posz_MC[i][ii])*(Hits_posz_MC[i][ii-1]-Hits_posz_MC[i][ii])
				  , 0.5 );
	Hits_dE_MC[i][ii] = Hits_E_MC[i][ii-1] - Hits_E_MC[i][ii];
	Hits_dEdx_MC[i][ii] = Hits_dE_MC[i][ii]/ Hits_dx_MC[i][ii];
      }
      
      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
      if (tpcid.isValid) { // Check if hit is in TPC
	if ( ii == 0 ) StartInTPC_MC[i] = 1; // Particle starts in TPC
	if ( ii == numberTrajectoryPoints-1 ) EndInTPC_MC[i] = 1; // Particle stops in TPC
	geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat);
	geo::TPCGeo      const& tpc  = cryo.TPC(tpcid.TPC);
	XPlanePosition      = tpc.PlaneLocation(0)[0];
	DriftTimeCorrection = fabs( Hits_posx_MC[i][ii] - XPlanePosition ) / XDriftVelocity;
	TimeAtPlane         = Hits_T_MC[i][ii] + DriftTimeCorrection;
	//std::cout << "Time at the APA " << TimeAtPlane << ", Compared to Gen Time " << Hits_T_MC[i][ii] << std::endl;
	if ( TimeAtPlane < detprop->TriggerOffset() || 
	     TimeAtPlane > detprop->TriggerOffset() + WindowSize ) {
	  //std::cout <<"!!!!!THIS TIME IS OUTSIDE READOUT WINDOW!!!"<< WindowSize << std::endl;
	  continue;
	} else { // If outside drift window...
	  //std::cout << "Got an entry with correct timing." << std::endl;
	  if (!insideActiveVolume) { // Check if first hit in TPC
	    zz = ii;
	    insideActiveVolume=true;
	  }		  
	  tmpPosition.GetXYZT(xyztArray);
	  zz2 = ii;
	  Hits_inTPC_MC[i][ii]=1;
	  ++nTPCHitsCounts; //Count MCHits within the TPC
	} // If inside Drift Window
      }
      //std::cout << ii << " " << Hits_inTPC_MC[i][ii] << " " << Hits_posx_MC[i][ii] << " " << Hits_posy_MC[i][ii] << " " << Hits_posz_MC[i][ii] << " " << zz << " " << zz2 << std::endl;
    }
    nTPCHits_MC[i] = nTPCHitsCounts;
    const TLorentzVector& positionStart = particle->Position(zz);
    TLorentzVector& positionEnd  =( TLorentzVector&)particle->Position(zz2);     
    TLorentzVector& momentumStart  =( TLorentzVector&)particle->Momentum(zz);

    mcang_x[i]=particle->Px(zz);
    mcang_y[i]=particle->Py(zz);
    mcang_z[i]=particle->Pz(zz);
    
    TLorentzVector tmpVec;
    TLorentzVector vectx(1,0,0,0);
    TLorentzVector vecty(0,1,0,0);
    TLorentzVector vectz(0,0,1,0);
    tmpVec=particle->Position(zz);
    mcpos_x[i]=(1./TMath::DegToRad())*tmpVec.Angle(vectx.Vect());
    mcpos_y[i]=(1./TMath::DegToRad())*tmpVec.Angle(vecty.Vect());
    mcpos_z[i]=(1./TMath::DegToRad())*tmpVec.Angle(vectz.Vect());
    
    float fMC_startXYZT[1000][4];
    float fMC_endXYZT[1000][4];  

    trkenergy_MC[i]=particle->E(zz);
    EnergyDeposited_MC[i] = particle->E(zz) - particle->E(zz2);
    trkmom_MC[i]=momentumStart.P();
    trkmom_XMC[i]=momentumStart.Px();
    trkmom_YMC[i]=momentumStart.Py();
    trkmom_ZMC[i]=momentumStart.Pz();
    trktheta_xz_MC[i] =  std::atan2(trkmom_XMC[i], trkmom_ZMC[i]);
    trktheta_yz_MC[i] = std::atan2(trkmom_YMC[i], trkmom_ZMC[i]);
    trketa_xy_MC[i] = std::atan2(trkmom_XMC[i], trkmom_YMC[i]);
    trketa_zy_MC[i] = std::atan2(trkmom_ZMC[i], trkmom_YMC[i]);
    TVector3 mcstartmom = particle->Momentum(zz).Vect();
    trktheta_MC[i]=mcstartmom.Theta();
    trkphi_MC[i]=mcstartmom.Phi();

    trkstartdoc_XMC[i]= pow ( (momentumStart.Px()*momentumStart.Px()) / ( trkmom_MC[i]* trkmom_MC[i]) , 0.5);
    trkstartdoc_YMC[i]= pow ( (momentumStart.Py()*momentumStart.Py()) / ( trkmom_MC[i]* trkmom_MC[i]) , 0.5);
    trkstartdoc_ZMC[i]= pow ( (momentumStart.Pz()*momentumStart.Pz()) / ( trkmom_MC[i]* trkmom_MC[i]) , 0.5);
    if ( trkmom_XMC[i] < 0 ) trkstartdoc_XMC[i] = -trkstartdoc_XMC[i];
    if ( trkmom_YMC[i] < 0 ) trkstartdoc_YMC[i] = -trkstartdoc_YMC[i];
    if ( trkmom_ZMC[i] < 0 ) trkstartdoc_ZMC[i] = -trkstartdoc_ZMC[i];
    positionStart.GetXYZT(fMC_startXYZT[i]); // In TPC
    positionEnd.GetXYZT(fMC_endXYZT[i]);     // In TPC
    trkstartx_MC[i]=fMC_startXYZT[i][0];
    trkstarty_MC[i]=fMC_startXYZT[i][1];
    trkstartz_MC[i]=fMC_startXYZT[i][2];
    trkendx_MC[i]=fMC_endXYZT[i][0];
    trkendy_MC[i]=fMC_endXYZT[i][1];
    trkendz_MC[i]=fMC_endXYZT[i][2];
    //std::cout << trkstartx_MC[i] << " " << trkstarty_MC[i] << " " << trkstartz_MC[i] << " " << trkendx_MC[i] << " " << trkendy_MC[i] << " " << trkendz_MC[i] << std::endl;

    for ( int LenLoop = zz+1; LenLoop<=zz2; ++LenLoop) {
      TPCLen_MC[i] +=  Hits_dx_MC[i][LenLoop];
    }
        
    ++i; // Increment number of MCParticles
    if (i == kMaxTrack) break;
 
  }// ipar
  nMCParticles = i;
  // -------------------------------------------
  // FINALLLY Fill Tree:
  // -------------------------------------------
  fTree->Fill();
  
}


void AnaTree::AnaTree::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("efield",efield,"efield[3]/F");
  fTree->Branch("t0",&t0,"t0/I");

  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("ntrkhits",ntrkhits,"ntrkhits[ntracks_reco]/I");
  fTree->Branch("trkid",trkid,"trkid[ntracks_reco]/I");  
  fTree->Branch("trkstartx",trkstartx,"trkstartx[ntracks_reco]/F");
  fTree->Branch("trkstarty",trkstarty,"trkstarty[ntracks_reco]/F");
  fTree->Branch("trkstartz",trkstartz,"trkstartz[ntracks_reco]/F");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/F");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/F");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/F");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/F");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/F");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/F");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/F");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/F");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/F");
  fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/F");
  fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/F");
  fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/F");
  fTree->Branch("trktheta_xz",trktheta_xz,"trktheta_xz[ntracks_reco]/F");
  fTree->Branch("trktheta_yz",trktheta_yz,"trktheta_yz[ntracks_reco]/F");
  fTree->Branch("trketa_xy",trketa_xy,"trketa_xy[ntracks_reco]/F");
  fTree->Branch("trketa_zy",trketa_zy,"trketa_zy[ntracks_reco]/F");
  fTree->Branch("trktheta",trktheta,"trktheta[ntracks_reco]/F");
  fTree->Branch("trkphi",trkphi,"trkphi[ntracks_reco]/F");
  fTree->Branch("trkd2",trkd2,"trkd2[ntracks_reco]/F");
  fTree->Branch("trkdedx2",trkdedx2,"trkdedx2[ntracks_reco][3][1000]/F");
  fTree->Branch("trkdqdx",trkdqdx,"trkdqdx[ntracks_reco][3][1000]/F");
  fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][3]/F");
  fTree->Branch("trkpitchHit",trkpitchHit,"trkpitchHit[ntracks_reco][3][1000]/F"); 
  fTree->Branch("trkkinE",trkkinE,"trkkinE[ntracks_reco][3]/F"); 
  fTree->Branch("trkrange",trkrange,"trkrange[ntracks_reco][3]/F"); 
  fTree->Branch("trkTPC",trkTPC,"trkTPC[ntracks_reco][3][1000]/F");
  fTree->Branch("trkplaneid",trkplaneid,"trkplaneid[ntracks_reco][3][1000]/F");
  fTree->Branch("trkresrg",trkresrg,"trkresrg[ntracks_reco][3][1000]/F");
  fTree->Branch("trkPosx",trkPosx,"trkPosx[ntracks_reco][3][1000]/F");
  fTree->Branch("trkPosy",trkPosy,"trkPosy[ntracks_reco][3][1000]/F");
  fTree->Branch("trkPosz",trkPosz,"trkPosz[ntracks_reco][3][1000]/F");
  fTree->Branch("trklen",trklen,"trklen[ntracks_reco]/F");
  fTree->Branch("trklen_L",trklen_L,"trklen_L[ntracks_reco]/F");
  fTree->Branch("trkdQdxSum",trkdQdxSum,"trkdQdxSum[ntracks_reco]/F");
  fTree->Branch("trkdQdxAverage",trkdQdxAverage,"trkdQdxAverage[ntracks_reco]/F");
  fTree->Branch("trkdEdxSum",trkdEdxSum,"trkdEdxSum[ntracks_reco]/F");
  fTree->Branch("trkdEdxAverage",trkdEdxAverage,"trkdEdxAverage[ntracks_reco]/F");
  
  fTree->Branch("trkMCTruthT0",trkMCTruthT0,"trkMCTruthT0[ntracks_reco]/F");
  fTree->Branch("trkMCTruthTrackID",trkMCTruthTrackID,"trkMCTruthTrackID[ntracks_reco]/I");
  fTree->Branch("trkPhotonCounterT0",trkPhotonCounterT0,"trkPhotonCounterT0[ntracks_reco]/F");
  fTree->Branch("trkPhotonCounterID",trkPhotonCounterID,"trkPhotonCounterID[ntracks_reco]/I");
  fTree->Branch("trkPhotonCounterConf",trkPhotonCounterConf,"trkPhotonCounterConf[ntracks_reco]/I");

  fTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
  fTree->Branch("trkid_MC",trkid_MC,"trkid_MC[nMCParticles]/I");
  fTree->Branch("trkpdg_MC",trkpdg_MC,"trkpdg_MC[nMCParticles]/I");
  fTree->Branch("trkndaughters_MC",trkndaughters_MC,"trkndaughters_MC[nMCParticles]/I");
  fTree->Branch("nTPCHits_MC",nTPCHits_MC,"nTPCHits_MC[nMCParticles]/I");
  fTree->Branch("StartInTPC_MC",StartInTPC_MC,"StartInTPC_MC[nMCParticles]/I");
  fTree->Branch("EndInTPC_MC",EndInTPC_MC,"EndInTPC_MC[nMCParticles]/I");
  fTree->Branch("trkMother_MC",trkMother_MC,"trkMother_MC[nMCParticles]/I");
  fTree->Branch("trkNumDaughters_MC",trkNumDaughters_MC,"trkNumDaughters_MC[nMCParticles]/I");
  fTree->Branch("trkFirstDaughter_MC",trkFirstDaughter_MC,"trkFirstDaughter_MC[nMCParticles]/I");
  fTree->Branch("trkPrimary_MC",trkPrimary_MC,"trkPrimarys_MC[nMCParticles]/I");
  fTree->Branch("StartTime_MC",StartTime_MC,"StartTime_MC[nMCParticles]/F");
  fTree->Branch("trkstartx_MC",trkstartx_MC,"trkstartx_MC[nMCParticles]/F");
  fTree->Branch("trkstarty_MC",trkstarty_MC,"trkstarty_MC[nMCParticles]/F");
  fTree->Branch("trkstartz_MC",trkstartz_MC,"trkstartz_MC[nMCParticles]/F");
  fTree->Branch("trkendx_MC",trkendx_MC,"trkendx_MC[nMCParticles]/F");
  fTree->Branch("trkendy_MC",trkendy_MC,"trkendy_MC[nMCParticles]/F");
  fTree->Branch("trkendz_MC",trkendz_MC,"trkendz_MC[nMCParticles]/F");
  fTree->Branch("trkenergy_MC",trkenergy_MC,"trkenergy_MC[nMCParticles]/F");
  fTree->Branch("EnergyDeposited_MC",EnergyDeposited_MC,"EnergyDeposited_MC[nMCParticles]/F");
  /* // MC information per hit....uses too much space though?
  fTree->Branch("Hits_posx_MC",Hits_posx_MC,"Hits_posx_MC[nMCParticles][1000]/F");
  fTree->Branch("Hits_posy_MC",Hits_posy_MC,"Hits_posy_MC[nMCParticles][1000]/F");
  fTree->Branch("Hits_posz_MC",Hits_posz_MC,"Hits_posz_MC[nMCParticles][1000]/F");
  fTree->Branch("Hits_inTPC_MC",Hits_inTPC_MC,"Hits_inTPC_MC[nMCParticles][1000]/I");
  fTree->Branch("Hits_mom_MC",Hits_mom_MC,"Hits_mom_MC[nMCParticles][1000]/F");
  fTree->Branch("Hits_momx_MC",Hits_momx_MC,"Hits_momx_MC[nMCParticles][1000]/F");
  fTree->Branch("Hits_momy_MC",Hits_momy_MC,"Hits_momy_MC[nMCParticles][1000]/F");
  fTree->Branch("Hits_momz_MC",Hits_momz_MC,"Hits_momz_MC[nMCParticles][1000]/F");
  fTree->Branch("Hits_E_MC",Hits_E_MC,"Hits_E_MC[nMCParticles][1000]/F");
  fTree->Branch("Hits_dEdx_MC",Hits_dEdx_MC,"HitsdEdx_MC[nMCParticles][1000]/F");
  */
  fTree->Branch("trkmom_MC",trkmom_MC,"trkmom_MC[nMCParticles]/F");
  fTree->Branch("trkmom_XMC",trkmom_XMC,"trkmom_XMC[nMCParticles]/F");
  fTree->Branch("trkmom_YMC",trkmom_YMC,"trkmom_YMC[nMCParticles]/F");
  fTree->Branch("trkmom_ZMC",trkmom_ZMC,"trkmom_ZMC[nMCParticles]/F");
  fTree->Branch("trkstartdoc_XMC",trkstartdoc_XMC,"trkstartdoc_XMC[nMCParticles]/F");
  fTree->Branch("trkstartdoc_YMC",trkstartdoc_YMC,"trkstartdoc_YMC[nMCParticles]/F");
  fTree->Branch("trkstartdoc_ZMC",trkstartdoc_ZMC,"trkstartdoc_ZMC[nMCParticles]/F");
  fTree->Branch("mcpos_x",&mcpos_x,"mcpos_x[nMCParticles]/F");
  fTree->Branch("mcpos_y",&mcpos_y,"mcpos_y[nMCParticles]/F");
  fTree->Branch("mcpos_z",&mcpos_z,"mcpos_z[nMCParticles]/F"); 
  fTree->Branch("mcang_x",&mcang_x,"mcang_x[nMCParticles]/F");
  fTree->Branch("mcang_y",&mcang_y,"mcang_y[nMCParticles]/F");
  fTree->Branch("mcang_z",&mcang_z,"mcang_z[nMCParticles]/F");
  fTree->Branch("trktheta_xz_MC",trktheta_xz_MC,"trktheta_xz_MC[nMCParticles]/F");
  fTree->Branch("trktheta_yz_MC",trktheta_yz_MC,"trktheta_yz_MC[nMCParticles]/F");
  fTree->Branch("trktheta_MC",trktheta_MC,"trktheta_MC[nMCParticles]/F");
  fTree->Branch("trkphi_MC",trkphi_MC,"trkphi_MC[nMCParticles]/F");
  fTree->Branch("trketa_xy_MC",trketa_xy_MC,"trketa_xy_MC[nMCParticles]/F");
  fTree->Branch("trketa_zy_MC",trketa_zy_MC,"trketa_zy_MC[nMCParticles]/F");
  fTree->Branch("trkTPCLen_MC",TPCLen_MC,"trkTPCLen_MC[nMCParticles]/F");
 
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("nhits2",&nhits2,"nhits2/I");
  fTree->Branch("nclust",&nclust,"nclust/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits2]/I");
  fTree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits2]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits2]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits2]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits2]/F");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits2]/F");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits2]/F");
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits2]/I");

  fTree->Branch("flash_total"  ,&flash_total ,"flash_total/I");
  fTree->Branch("flash_time"   ,flash_time   ,"flash_time[flash_total]/F");
  fTree->Branch("flash_width"  ,flash_width  ,"flash_width[flash_total]/F");
  fTree->Branch("flash_abstime",flash_abstime,"flash_abstime[flash_total]/F");
  fTree->Branch("flash_YCentre",flash_YCentre,"flash_YCentre[flash_total]/F");
  fTree->Branch("flash_YWidth" ,flash_YWidth ,"flash_YWidth[flash_total]/F");
  fTree->Branch("flash_ZCentre",flash_ZCentre,"flash_ZCentre[flash_total]/F");
  fTree->Branch("flash_ZWidth" ,flash_ZWidth ,"flash_ZWidth[flash_total]/F");
  fTree->Branch("flash_TotalPE",flash_TotalPE,"flash_TotalPE[flash_total]/F");

  fTree->Branch("ntrigs",&ntrigs,"ntrigs/I");
  fTree->Branch("trig_time",trig_time,"trig_time[ntrigs]/I");
  fTree->Branch("trig_id",trig_id,"trig_id[ntrigs]/I");
}

//void AnaTree::AnaTree::reconfigure(fhicl::ParameterSet const & p)
//{
//  // Implementation of optional member function here.
//}
void AnaTree::AnaTree::ResetVars(){

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  for (int i = 0; i<3; ++i){
    efield[i] = -99999;
  }
  t0 = -99999;
  ntracks_reco = -99999;
  nMCParticles = -99999;

  for (int i = 0; i < kMaxTrack; ++i){
    trkstartx[i] = -99999;
    trkstarty[i] = -99999;
    trkstartz[i] = -99999;
    trkendx[i] = -99999;
    trkendy[i] = -99999;
    trkendz[i] = -99999;
    trkstartx_MC[i] = -99999;
    trkstarty_MC[i] = -99999;
    trkstartz_MC[i] = -99999;
    trkendx_MC[i] = -99999;
    trkendy_MC[i] = -99999;
    trkendz_MC[i] = -99999;
    TPCLen_MC[i] = 0;
    trkenergy_MC[i] = -99999;
    EnergyDeposited_MC[i] = -99999;
    trkmom_MC[i] = -99999;
    trkmom_XMC[i] = -99999;
    trkmom_YMC[i] = -99999;
    trkmom_ZMC[i] = -99999;
    trkstartdoc_XMC[i] = -99999;
    trkstartdoc_YMC[i] = -99999;
    trkstartdoc_ZMC[i] = -99999;
    trkid_MC[i] = -99999;
    trkpdg_MC[i] = -99999;
    trkndaughters_MC[i] = -99999;
    trkMother_MC[i] = -99999;
    trkNumDaughters_MC[i] = -99999;
    trkFirstDaughter_MC[i] = -99999;
    trkPrimary_MC[i] = 0;
    trkd2[i] = -99999;
    trktheta_xz_MC[i] = -99999;
    trktheta_yz_MC[i] = -99999;
    trktheta_MC[i] = -99999;
    trkphi_MC[i] = -99999;
    trketa_xy_MC[i] = -99999;
    trketa_zy_MC[i] = -99999;
    trktheta_xz[i] = -99999;
    trktheta_yz[i] = -99999;
    trketa_xy[i] = -99999;
    trketa_zy[i] = -99999;
    trktheta[i] = -99999;
    trkphi[i] = -99999;
    mcpos_x[i] = -99999;
    mcpos_y[i] = -99999;
    mcpos_z[i] = -99999;
    mcang_x[i] = -99999;
    mcang_y[i] = -99999;
    mcang_z[i] = -99999;
    trkdQdxSum[i] = 0;
    trkdEdxSum[i] = 0;
    nTPCHits_MC[i]=0;
    StartTime_MC[i]=0;
    StartInTPC_MC[i] = 0;
    EndInTPC_MC[i] = 0;
    for(int ii=0;ii<3;ii++)
      {
	trkkinE[i][ii] = -99999;
	trkrange[i][ii] = -99999;
	for(int k=0;k<1000;k++)
	  {
	    trkdedx2[i][ii][k] = -99999;
	    trkdqdx[i][ii][k] = -99999;
	    trkpitchHit[i][ii][k] = -99999;
	    trkTPC[i][ii][k] = -99999;
	    trkplaneid[i][ii][k] = -99999;
	    trkresrg[i][ii][k] = -99999;
	    trkPosx[i][ii][k]  = -99999;
	    trkPosy[i][ii][k]  = -99999;
	    trkPosz[i][ii][k]  = -99999;
	  }
      }
    trkstartdcosx[i] = -99999;
    trkstartdcosy[i] = -99999;
    trkstartdcosz[i] = -99999;
    trklen[i] = -99999;
    trklen_L[i] = -99999;
    trkid[i] = -99999;
    trkMCTruthT0[i] = -99999;
    trkMCTruthTrackID[i] = -99999;
    trkPhotonCounterT0[i]= -99999;
    trkPhotonCounterID[i]= -99999;
    trkPhotonCounterConf[i]= -99999;

    trkenddcosx[i] = -99999;
    trkenddcosy[i] = -99999;
    trkenddcosz[i] = -99999;
    ntrkhits[i] = -99999;
    for (int j = 0; j<kMaxTrackHits; ++j){
      trkx[i][j] = -99999;
      trky[i][j] = -99999;
      trkz[i][j] = -99999;
      Hits_posx_MC[i][j] = -99999;
      Hits_posy_MC[i][j] = -99999;
      Hits_posz_MC[i][j] = -99999;
      Hits_inTPC_MC[i][j] = 0;
      Hits_mom_MC[i][j] = -99999;
      Hits_momx_MC[i][j] = -99999;
      Hits_momy_MC[i][j] = -99999;
      Hits_momz_MC[i][j] = -99999;
      Hits_E_MC[i][j] = -99999;
      Hits_T_MC[i][j] = -99999;
      Hits_dEdx_MC[i][j] = -9999;
    }
    for (int j = 0; j<3; ++j){
      trkpitch[i][j] = -99999;
    }
  }
  nhits = -99999;
  nhits2= 0;
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -99999;
    hit_tpc[i] = -99999;
    hit_wire[i] = -99999;
    hit_channel[i] = -99999;
    hit_peakT[i] = -99999;
    hit_charge[i] = -99999;
    hit_ph[i] = -99999;
    hit_trkid[i] = -99999;
  }
  flash_total = 0;
  for (int f = 0; f < kMaxFlash; ++f) {
    flash_time[f]    = -9999;
    flash_width[f]   = -9999;
    flash_abstime[f] = -9999;
    flash_YCentre[f] = -9999;
    flash_YWidth[f]  = -9999;
    flash_ZCentre[f] = -9999;
    flash_ZWidth[f]  = -9999;
    flash_TotalPE[f] = -9999;
  }
}

void AnaTree::AnaTree::endJob()
{
}


DEFINE_ART_MODULE(AnaTree::AnaTree)
