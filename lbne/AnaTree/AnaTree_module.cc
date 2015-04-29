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
//
////////////////////////////////////////////////////////////////////////
// Framework includes
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
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RawData/ExternalTrigger.h"
#include "MCCheater/BackTracker.h"
#include "AnalysisBase/Calorimetry.h"
#include "AnalysisBase/T0.h"
#include "AnalysisBase/ParticleID.h"



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

namespace AnaTree {
  class AnaTree;
}
namespace {

  // Local functions.

  // Calculate distance to boundary.

  //----------------------------------------------------------------------------
  // removed

  // Length of reconstructed track.
  //----------------------------------------------------------------------------
  /*
  double length(const recob::Track& track)
  {
    double result = 0.;
    TVector3 disp = track.LocationAtPoint(0);
    int n = track.NumberTrajectoryPoints();
    for(int i = 1; i < n; ++i) {
      const TVector3& pos = track.LocationAtPoint(i);

      disp -= pos;

      result += disp.Mag();

      disp = pos;
    }
    //    mf::LogVerbatim("output") << " length (track) " << result;
    return result;
  }
  */
  // Length of MC particle.
  //----------------------------------------------------------------------------
  double length(const simb::MCParticle& part, double dx,
                TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
                unsigned int /*tpc*/ = 0, unsigned int /*cstat*/ = 0)
  {

    // Get services.
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;

    double result = 0.;

    TVector3 disp;

    int n = part.NumberTrajectoryPoints();

    bool first = true;
    //    std::cout<< " n is " << n << std::endl;
    for(int i = 0; i < n; ++i) {

      TVector3 pos = part.Position(i).Vect();

      // Make fiducial cuts.  Require the particle to be within the physical volume of
      // the tpc, and also require the apparent x position to be within the expanded
      // readout frame.

      double const tmpArray[]={pos.X(),pos.Y(),pos.Z()};
      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpArray);
      if (!tpcid.isValid) continue;
      
      pos[0] += dx;
      
      double ticks;
      
      ticks = detprop->ConvertXToTicks(pos[0], 0, tpcid.TPC, tpcid.Cryostat);

	//if(ticks >5e+12)
	// continue;
        if(ticks >= 0. && ticks < detprop->ReadOutWindowSize()) {
          if(first) {
            start = pos;
            startmom = part.Momentum(i).Vect();
	   
          }
          else {
            disp -= pos;
            result += disp.Mag();
          }
          first = false;
          disp = pos;
          end = pos;
          endmom = part.Momentum(i).Vect();
        }
	
    }

    //    mf::LogVerbatim("output") << " length (MCParticle) " << result;

    return result;

  }

  // Fill efficiency histogram assuming binomial errors.

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
  double efield[3];
  int t0;
  int trigtime[16];
  int ntracks_reco;         //number of reconstructed tracks
  double trkstartx[kMaxTrack];
  double trkstarty[kMaxTrack];
  double trkstartz[kMaxTrack];
  double trkendx[kMaxTrack];
  double trkendy[kMaxTrack];
  double trkendz[kMaxTrack];

  int nMCParticles;
  int nCryoHits_MC[kMaxTrack];
  double StartTime_MC[kMaxTrack];
  double trkstartx_MC[kMaxTrack];
  double trkstarty_MC[kMaxTrack];
  double trkstartz_MC[kMaxTrack];
  double trkendx_MC[kMaxTrack];
  double trkendy_MC[kMaxTrack];
  double trkendz_MC[kMaxTrack];
  double trklen_MC[kMaxTrack];
  double trklen_cut_MC[kMaxTrack];
  double Hits_posx_MC[kMaxTrack][1000];
  double Hits_posy_MC[kMaxTrack][1000];
  double Hits_posz_MC[kMaxTrack][1000];
  double Hits_mom_MC[kMaxTrack][1000];
  double Hits_E_MC[kMaxTrack][1000];

  int    trkid_MC[kMaxTrack];
  int    trkpdg_MC[kMaxTrack];
  double trkmom_MC[kMaxTrack];
  double trkmom_XMC[kMaxTrack];
  double trkmom_YMC[kMaxTrack];
  double trkmom_ZMC[kMaxTrack];
  double trkenergy_MC[kMaxTrack];
  double trkstartdoc_XMC[kMaxTrack];
  double trkstartdoc_YMC[kMaxTrack];
  double trkstartdoc_ZMC[kMaxTrack];
  int    trkMother_MC[kMaxTrack];
  int    trkNumDaughters_MC[kMaxTrack];
  int    trkFirstDaughter_MC[kMaxTrack];
  int    trkLastDaughter_MC[kMaxTrack];
  int    trkPrimary_MC[kMaxTrack];
  double trkd2[kMaxTrack];
  double trkcolin[kMaxTrack];
  double trklen[kMaxTrack];
  double trklen_L[kMaxTrack];
  int trkid[kMaxTrack];
  double trktheta_xz_MC[kMaxTrack];
  double trktheta_yz_MC[kMaxTrack];
  double trketa_xy_MC[kMaxTrack];
  double trketa_zy_MC[kMaxTrack];
  double trktheta_MC[kMaxTrack];
  double trkphi_MC[kMaxTrack];

  double trkMCTruthT0[kMaxTrack];
  double trkMCTruthPdG[kMaxTrack];

  double trktheta_xz[kMaxTrack];
  double trketa_xy[kMaxTrack];
  double trktheta_yz[kMaxTrack];
  double trketa_zy[kMaxTrack];
  double trktheta[kMaxTrack];
  double trkphi[kMaxTrack];
  double trkdQdxSum[kMaxTrack];
  double trkdQdxAverage[kMaxTrack];
  double trkdEdxSum[kMaxTrack];
  double trkdEdxAverage[kMaxTrack];

  double trkdedx[kMaxTrack];
  double trkdedx2[kMaxTrack][3][1000];
  double trkdqdx[kMaxTrack][3][1000];
  double trkpitchHit[kMaxTrack][3][1000];
  double trkkinE[kMaxTrack][3];
  double trkrange[kMaxTrack][3];
  double trkTPC[kMaxTrack][3][1000];
  double trkplaneid[kMaxTrack][3][1000];
  double trkresrg[kMaxTrack][3][1000];
  double trkPosx[kMaxTrack][3][1000];
  double trkPosy[kMaxTrack][3][1000];
  double trkPosz[kMaxTrack][3][1000]; 
  double trkdedx_MC[kMaxTrack];
  double trkdq_MC[kMaxTrack];
  double mcang_x[kMaxTrack];
  double mcang_y[kMaxTrack];
  double mcang_z[kMaxTrack];
  double mcpos_x[kMaxTrack];
  double mcpos_y[kMaxTrack];
  double mcpos_z[kMaxTrack];

  double trkang[kMaxTrack];
  double trkcolinearity[kMaxTrack];
  double trkmatchdisp[kMaxTrack];
  double trkwmatchdisp[kMaxTrack];
  double trklenratio[kMaxTrack];
  double trkstartdcosx[kMaxTrack];
  double trkstartdcosy[kMaxTrack];
  double trkstartdcosz[kMaxTrack];
  double trkenddcosx[kMaxTrack];
  double trkenddcosy[kMaxTrack];
  double trkenddcosz[kMaxTrack];
  int    ntrkhits[kMaxTrack];
  double trkx[kMaxTrack][kMaxTrackHits];
  double trky[kMaxTrack][kMaxTrackHits];
  double trkz[kMaxTrack][kMaxTrackHits];
  
  double trkpitch[kMaxTrack][3];
  int    nhits;
  int nclust;

  int  hit_tpc[kMaxHits];
  int    hit_plane[kMaxHits];
  int    hit_wire[kMaxHits];
  int    hit_channel[kMaxHits];
  double hit_peakT[kMaxHits];
  double hit_charge[kMaxHits];
  double hit_ph[kMaxHits];
  int    hit_trkid[kMaxHits];


  std::string fTrigModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fClusterModuleLabel;
  std::string fTrkSpptAssocModuleLabel;
  std::string fHitSpptAssocModuleLabel;
  std::string fSimulationProducerLabel; 
  std::string fCalorimetryModuleLabel; 
  std::string fMCTruthT0Modulelabel;





  int fDump;                 // Number of events to dump to debug message facility.
  int fPdg;
  double fMinMCKE;           // Minimum MC particle kinetic energy (GeV).
  double fMinMCLen;          // Minimum MC particle length in tpc (cm).
  double fMatchColinearity;  // Minimum matching colinearity.
  double fMatchDisp;         // Maximum matching displacement.
  double fWMatchDisp;        // Maximum matching displacement in the w direction.
  bool fIgnoreSign;          // Ignore sign of mc particle if true.
  bool fStitchedAnalysis;    // if true, do the whole drill-down from stitched track to assd hits

  double fElectronsToGeV; // conversion factor
  art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

  //  std::map<int, MCHists> fMCHistMap;       // Indexed by pdg id.
  //  std::map<int, RecoHists> fRecoHistMap;   // Indexed by pdg id.

};


AnaTree::AnaTree::AnaTree(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fTrigModuleLabel       (pset.get< std::string >("TrigModuleLabel"))
  , fHitsModuleLabel       (pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel       (pset.get< std::string >("TrackModuleLabel"))
  , fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel"))
  , fTrkSpptAssocModuleLabel    (pset.get< std::string >("TrkSpptAssocModuleLabel"))
  , fHitSpptAssocModuleLabel    (pset.get< std::string >("HitSpptAssocModuleLabel"))
  , fSimulationProducerLabel ( pset.get< std::string >("SimulationLabel"))
  , fCalorimetryModuleLabel ( pset.get< std::string >("CalorimetryModuleLabel"))
  , fMCTruthT0Modulelabel     ( pset.get< std::string >("MCTruthT0Modulelabel"))
  , fDump              (pset.get<int>("Dump"))
  , fPdg              (pset.get<int>("pdg"))
  , fMinMCKE            (pset.get<double>("MinMCKE"))
  , fMinMCLen           (pset.get<double>("MinMCLen"))
  , fMatchColinearity       (pset.get<double>("MatchColinearity"))
  , fMatchDisp             (pset.get<double>("MatchDisp"))
  , fWMatchDisp             (pset.get<double>("WMatchDisp"))
  , fIgnoreSign             (pset.get<bool>("IgnoreSign"))
  , fStitchedAnalysis       (pset.get<bool>("StitchedAnalysis",false))
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
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();


  efield[0] = larprop->Efield(0);
  efield[1] = larprop->Efield(1);
  efield[2] = larprop->Efield(2);
  
  t0 = detprop->TriggerOffset();
  
  art::Handle< std::vector<raw::ExternalTrigger> > trigListHandle;
  std::vector<art::Ptr<raw::ExternalTrigger> > triglist;
  if (evt.getByLabel(fTrigModuleLabel,trigListHandle))
    art::fill_ptr_vector(triglist, trigListHandle);
  
  for (size_t i = 0; i<triglist.size(); ++i){
    trigtime[i] = triglist[i]->GetTrigTime();
  }
  
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);    
  
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);
  
  
  art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
  evt.getByLabel(fSimulationProducerLabel, simChannelHandle);
 
  //  Event.getByLabel(fSimulationProducerLabel, particleHandle);
  //  art::Handle< std::vector<simb::MCParticle> > particleHandle;

  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  art::FindMany<recob::Track>       fmtk(hitListHandle, evt, fTrackModuleLabel);

  art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  art::FindMany<anab::T0>           fmt0( trackListHandle, evt, fMCTruthT0Modulelabel);

  //track information
  ntracks_reco=tracklist.size();

  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  
  // Get Cryostat information.........Only taking one cryo atm....
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
  std::cout << boundaries[0] << " " << boundaries[1] << " " << boundaries[2] << " " << boundaries[3] << " " <<boundaries[4] << " " << boundaries[5] << std::endl;

  // **********************
  // **********************
  //
  //  Trackh:
  //  Trackvh:
  //
  //
  // *********************
  // *********************
  
  art::Handle< std::vector<recob::Track> > trackh;
  evt.getByLabel(fTrackModuleLabel, trackh);
  
  
  art::Handle< std::vector< art::PtrVector < recob::Track > > > trackvh;
  evt.getByLabel(fTrackModuleLabel, trackvh);
  
  // Protect against invalid art::Handle (both here and when trackvh is used below)
  // TODO: How to do this for art::PtrVector without the uninitialised iterator?
  std::vector< art::PtrVector<recob::Track> >::const_iterator cti; 
  if (trackvh.isValid()) cti = trackvh->begin();                   
  
  
  
  // **********************
  // **********************
  //
  //  TrackList:
  //
  //
  //
  // *********************
  // *********************
  
  for(int i=0; i<std::min(int(tracklist.size()),kMaxTrack);++i){
    
    trkid[i]=i;
    trackStart.clear();
    trackEnd.clear();
    memset(larStart, 0, 3);
    memset(larEnd, 0, 3);
    tracklist[i]->Extent(trackStart,trackEnd); 
    tracklist[i]->Direction(larStart,larEnd);
    trkstartx[i]        = trackStart[0];
    trkstarty[i]        = trackStart[1];
    trkstartz[i]        = trackStart[2];
    trkendx[i]        = trackEnd[0];
    trkendy[i]        = trackEnd[1];
    trkendz[i]        = trackEnd[2];
    trkstartdcosx[i]  = larStart[0];
    trkstartdcosy[i]  = larStart[1];
    trkstartdcosz[i]  = larStart[2];
    TLorentzVector v1(trackStart[0],trackStart[1],trackStart[2],0);
    TLorentzVector v2(trackEnd[0],trackEnd[1],trackEnd[2],0);
    trklen[i]=(v2-v1).Rho();
    //trkdedx[i]=trkde/trklen[i];
    //    trkang[i]=TMath::Cos((v2-v1).Angle(tmpVec.Vect()));
    //trkang[i]=TMath::Cos((v2-v1).Angle(tmpVec.Vect()));
    //    trklen[i]=v.Mag();
    trkenddcosx[i]    = larEnd[0];
    trkenddcosy[i]    = larEnd[1];
    trkenddcosz[i]    = larEnd[2];
    //    ntrkhits[i] = fmsp.at(i).size();
    
    //    std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
    //   art::Ptr<recob::Track> pptrack(trackh, i);
    //    auto pp { pptrack };
    //    art::FindManyP<recob::SpacePoint> spptAssns(pp, evt, fTrackModuleLabel); 
    //    ntrkhits[i] = spptAssns.at(0).size();   
    //    int nnn = spptAssns.at(0).size();
    //    std::cout << " Number of clumps ----> " << nnn<< std::endl;
    
    double distance_squared=0;
    TVector3 V1(trackStart[0],trackStart[1],trackStart[2]);
    TVector3 V2(trackEnd[0],trackEnd[1],trackEnd[2]);
    TVector3 vOrth=(V2-V1).Orthogonal();
    TVector3 pointVector=V1;
    
    /* if(trackvh.isValid())
      {
      int k=i;
      int ntrackhits=0;
      const art::PtrVector<recob::Track> pvtrack(*(cti++));
      //        auto it = pvtrack.begin();
      
      int ntrack = pvtrack.size();
      art::FindManyP<recob::SpacePoint> fs( pvtrack, evt, fTrkSpptAssocModuleLabel);
      double distance=0;
      for(int ii=0;ii<ntrack;ii++){
      
      for (size_t j = 0; j<fs.at(ii).size(); ++j){
      
      ntrackhits++;
      TVector3 sptVector(fs.at(ii).at(j)->XYZ()[0],fs.at(ii).at(j)->XYZ()[1],fs.at(ii).at(j)->XYZ()[2]);
      TVector3 vToPoint=sptVector-pointVector;
      distance=(vOrth.Dot(vToPoint))/vOrth.Mag();
      if(isnan(distance)){
      mf::LogVerbatim("output") <<"is nan" <<(vOrth.Dot(vToPoint));
      mf::LogVerbatim("output") <<"is nan" <<vOrth.Mag();
      }
      distance_squared+=distance *distance;
      trkx[k][ntrackhits] = fs.at(ii).at(j)->XYZ()[0];
      trky[k][ntrackhits] = fs.at(ii).at(j)->XYZ()[1];
      trkz[k][ntrackhits] = fs.at(ii).at(j)->XYZ()[2];
      }
      }
      ntrkhits[k]=ntrackhits;
      distance_squared=distance_squared/ntrkhits[k];
      if(!isnan(distance_squared))
      trkd2[k]=distance_squared;
      }*/
    //    else
    if(fmsp.isValid() ){
      ntrkhits[i] = fmsp.at(i).size();
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
    
    // *********************
    //  Calorimetric stuff:
    //  
    // *********************
    //
    //
    //
    //
    //
    //
    
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
        
    // *********************
    //  End Calorimetric stuff
    //  
    // *********************
    
    //***************
    //   T0 stuff 
    //***************

    if ( fmt0.isValid() ) {
      std::vector<const anab::T0*> T0s = fmt0.at(i);
      for (size_t t0size =0; t0size < T0s.size(); t0size++) { // Possibly unnecessary as for MCTruth at least size=0
	if (T0s[t0size] -> TriggerType() == 0 ) { // If Trigger 0
	} // Trig 0
	if (T0s[t0size] -> TriggerType() == 1 ) { // If Counters
	  //std::cout << "Using Counter Information!!!! T0 " << T0s[t0size]->Time() << ", Bits " << T0s[t0size]->TriggerBits() << ", Size " << T0s[t0size]->ID() << std::endl;	  
	} // Trig 1
	if (T0s[t0size] -> TriggerType() == 2 ) { // If MCTruth
	  //std::cout << "Using Trigger Type Monte Carlo Truth!!! T0 " << T0s[t0size]->Time() << ", PdG " << T0s[t0size]->TriggerBits() << ", Size " << T0s[t0size]->ID() << std::endl;
	  trkMCTruthT0[i]  = T0s[t0size]->Time();
	  trkMCTruthPdG[i] = T0s[t0size]->TriggerBits();
	} // Trig 2
      } // T0 size
    } // T0 valid

    // *********************
    //   Cuts specific quantities
    //  
    // *********************
    //
    //
    //
    //
    TMatrixD rot(3,3);
    int start_point =0;
    tracklist[i]->GlobalToLocalRotationAtPoint(start_point, rot);
    //  int ntraj = tracklist[i]->NumberTrajectoryPoints();
    // if(ntraj > 0) {
    TVector3 pos = tracklist[i]->Vertex();
    art::ServiceHandle<cheat::BackTracker> bktrk;
    const   sim::ParticleList& ppplist=bktrk->ParticleList();
    std::vector<const simb::MCParticle*> plist2;
    plist2.reserve(ppplist.size());
    
    art::Ptr<recob::Track> ptrack(trackh, i);
    const recob::Track& track = *ptrack;
    //trklenratio[i] = length(track)/plen;
    //trklen_L[i]=length(track);
    trklen_L[i]=track.Length(); // exactly the same as old function delcared at top, but this accesses track information.
    //
    //**********************
    //
    // End Cut specific quantities
    //  
    //***********************

    
    
    
    for (int j = 0; j<3; ++j){
      try {
	if (j==0)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kU);
	else if (j==1)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kV);
	else if (j==2)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kZ);
      }
      catch( cet::exception &e) {
	mf::LogWarning("AnaTree")<<"caught exeption "<<e<<"\n setting pitch to 0";
	trkpitch[i][j] = 0;
      }
    }
    
    /*   simb::MCParticle* particle=0;
    for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar) {
      particle = ipar->second;
    }
    int pdg = particle->PdgCode();
    if (abs(pdg)!=abs(fPdg)) continue;
    TVector3 startmom;
    startmom=particle->Momentum(0).Vect();
    TVector3 mcmomltemp=rot * startmom;
    trkcolin[i]=mcmomltemp.Z()/mcmomltemp.Mag();
    */
    double recotime = 0.;
    double trackdx = recotime * 1.e-3 * larprop->DriftVelocity();  // cm
    // Fill histograms involving reco tracks only.
    int ntraj = track.NumberTrajectoryPoints();
    if(ntraj > 0) {
      TVector3 pos = track.Vertex();
      TVector3 dir = track.VertexDirection();
      TVector3 end = track.End();
      pos[0] += trackdx;
      end[0] += trackdx;
      trktheta_xz[i] = std::atan2(dir.X(), dir.Z());
      trktheta_yz[i] = std::atan2(dir.Y(), dir.Z());
      trketa_xy[i] = std::atan2(dir.X(), dir.Y());
      trketa_zy[i] = std::atan2(dir.Z(), dir.Y());
      trktheta[i]=dir.Theta();
      trkphi[i]=dir.Phi();
    }
  }
  
  // -----------------------------------------------
  // NOW DO THE CLUSTER/HIT STUFF.......
  // -----------------------------------------------
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  nhits = hitlist.size();
  nclust=clusterlist.size();
  for (int i = 0; i<std::min(int(hitlist.size()),kMaxHits); ++i){
    unsigned int channel = hitlist[i]->Channel();
    geo::WireID wireid = hitlist[i]->WireID();
    hit_tpc[i]     =wireid.TPC;
    hit_plane[i]   = wireid.Plane;
    hit_wire[i]    = wireid.Wire;
    hit_channel[i] = channel;
    hit_peakT[i]   = hitlist[i]->PeakTime();
    hit_charge[i]  = hitlist[i]->Integral();
    hit_ph[i]      = hitlist[i]->PeakAmplitude();
    if (fmtk.at(i).size()!=0){
      hit_trkid[i] = fmtk.at(i)[0]->ID();
    }
  }
  
  // -----------------------------------------------
  // NOW DO ALL THE MONTE CARLO TRUTH STUFF.......
  // -----------------------------------------------
  art::ServiceHandle<cheat::BackTracker> bt;
  //  int trkid=1;
  const sim::ParticleList& plist = bt->ParticleList();
  simb::MCParticle *particle=0;
  //  const simb::MCParticle *particle = bt->TrackIDToParticle(trkid);
  int i=0; // particle index
  for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    particle = ipar->second;

    // Line below selects only primaries, have removed so get all particles.
    // So want to add if Process=="Primary" in macro, this way can see particles other than primaries too.
    //if(!(particle->Process()=="primary" && abs(particle->PdgCode())== abs(fPdg))) continue;
        
    mcang_x[i]=particle->Px();
    mcang_y[i]=particle->Py();
    mcang_z[i]=particle->Pz();
    
    TLorentzVector tmpVec;
    TLorentzVector vectx(1,0,0,0);
    TLorentzVector vecty(0,1,0,0);
    TLorentzVector vectz(0,0,1,0);
    tmpVec=particle->Position();
    mcpos_x[i]=(1./TMath::DegToRad())*tmpVec.Angle(vectx.Vect());
    mcpos_y[i]=(1./TMath::DegToRad())*tmpVec.Angle(vecty.Vect());
    mcpos_z[i]=(1./TMath::DegToRad())*tmpVec.Angle(vectz.Vect());
    
   
    //    int fMC_Ntrack = particles2.size();
    float fMC_startXYZT[1000][4];
    float fMC_endXYZT[1000][4];
    
    //double trkde=0;
    
    
    //for (auto const& particle: particles2 ) {
    //        int Ndaughters = particle->NumberDaughters();
    //        vector<int> daughters;
    //        for (int i=0; i<Ndaughters; i++) {
    //            daughters.push_back(particle->Daughter(i));
    //        }
    //        fMC_daughters.push_back(daughters);
    size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
    trkid_MC[i]=particle->TrackId();
    trkpdg_MC[i]=particle->PdgCode();
    trkMother_MC[i]=particle->Mother();
    trkNumDaughters_MC[i]=particle->NumberDaughters();
    trkFirstDaughter_MC[i]=particle->FirstDaughter();
    //trkLastDaughter_MC[i]=particle->LastDaughter();
    if (particle->Process() == "primary" ) trkPrimary_MC[i] = 1;
    else trkPrimary_MC[i] = 0;
    
    
    //int trackID=particle->TrackId();
    
    /*
      for ( auto const& channel : (*simChannelHandle) ) {
      auto channelNumber = channel.Channel();
      if ( fGeometry->SignalType( channelNumber ) == geo::kCollection ) {
      auto const& timeSlices = channel.TDCIDEMap();
      for ( auto const& timeSlice : timeSlices ) {
      auto const& energyDeposits = timeSlice.second;
      for ( auto const& energyDeposit : energyDeposits ) {
      if ( energyDeposit.trackID == trackID ) {
      
      trkde+=energyDeposit.numElectrons*1./( 1e6/23.6);//MeV
      }//TrackID
      }//energyDeposit
      }//timeSlice
      }//CollectionPlane
      }//simChannel
      trkdq_MC[i]=trkde*1e6/23.6;//back to q
    */
    double xyztArray[4];
    int zz =0;
    int zz2 =0;
    bool insideActiveVolume=false;
    int nCryoHitsCounts =0;
    /*    
    double origin[3] = {0.};
    double world[3] = {0.};
    double cryoBound_pos[3];
    double cryoBound_neg[3];
    int c=0;//only one cryo
    
    geom->Cryostat(c).LocalToWorld(origin, world);
    cryoBound_neg[0]=world[0] - geom->Cryostat(c).HalfWidth();
    cryoBound_neg[1]=world[1] - geom->Cryostat(c).HalfWidth();
    cryoBound_neg[2]=world[2] - geom->Cryostat(c).HalfWidth();
    
    cryoBound_pos[0]=world[0] + geom->Cryostat(c).HalfWidth();
    cryoBound_pos[1]=world[1] + geom->Cryostat(c).HalfWidth();
    cryoBound_pos[2]=world[2] + geom->Cryostat(c).HalfWidth();
    
    for(size_t ii=0;ii<numberTrajectoryPoints;ii++) {
      const TLorentzVector& tmpPosition=particle->Position(ii);
      //tmpPosition.GetXYZT(xyztArray);   
      if((tmpPosition[0]<cryoBound_pos[0]) && (tmpPosition[0]>cryoBound_neg[0])) {
	if((tmpPosition[1]<cryoBound_pos[1]) && (tmpPosition[1]>cryoBound_neg[1])) {
	  if((tmpPosition[2]<cryoBound_pos[2]) && (tmpPosition[2]>cryoBound_neg[2])) {
    */
    for(size_t ii=0;ii<numberTrajectoryPoints;ii++) {
      const TLorentzVector& tmpPosition=particle->Position(ii);
      //tmpPosition.GetXYZT(xyztArray);
      if((tmpPosition[0]>boundaries[0]) && (tmpPosition[0]<boundaries[1])) {
	if((tmpPosition[1]>boundaries[2]) && (tmpPosition[1]<boundaries[3])) {
	  if((tmpPosition[2]>boundaries[4]) && (tmpPosition[2]<boundaries[5])) {
	    if (!insideActiveVolume) {
	      zz = ii;
	      //std::cout << "Now particle is in cryostat " << std::endl;
	      insideActiveVolume=true;
	    }		  
	    //std::cout << "Temp Pos " << tmpPosition[0] << ", " <<  tmpPosition[1] << ", " <<  tmpPosition[2] << std::endl;
	    tmpPosition.GetXYZT(xyztArray);
	    zz2 = ii;
	    ++nCryoHitsCounts; //Count MCHits within the cryostat - note this does not mean they are all in TPC's! 
	    Hits_posx_MC[i][ii] = tmpPosition[0];
	    Hits_posy_MC[i][ii] = tmpPosition[1];
	    Hits_posz_MC[i][ii] = tmpPosition[2];
	    Hits_mom_MC[i][ii]  = particle->P(ii);
	    Hits_E_MC[i][ii]    = particle->E(ii);
	  }
	}
      }
      
      if ( (insideActiveVolume) && (zz2 != (int)ii) ) break;
    }
    nCryoHits_MC[i] = nCryoHitsCounts;
    const TLorentzVector& positionStart = particle->Position(zz);
    TLorentzVector& positionEnd  =( TLorentzVector&)particle->Position(zz2);     
    //        const TLorentzVector& momentumStart = particle->Momentum(0);
    //        const TLorentzVector& momentumEnd   = particle->Momentum(last);
    TLorentzVector& momentumStart  =( TLorentzVector&)particle->Momentum(zz);
    trkenergy_MC[i]=particle->E();
    trkmom_MC[i]=momentumStart.P();
    trkmom_XMC[i]=momentumStart.Px();
    trkmom_YMC[i]=momentumStart.Py();
    trkmom_ZMC[i]=momentumStart.Pz();
    trkstartdoc_XMC[i]= pow ( (momentumStart.Px()*momentumStart.Px()) / ( trkmom_MC[i]* trkmom_MC[i]) , 0.5);
    trkstartdoc_YMC[i]= pow ( (momentumStart.Py()*momentumStart.Py()) / ( trkmom_MC[i]* trkmom_MC[i]) , 0.5);
    trkstartdoc_ZMC[i]= pow ( (momentumStart.Pz()*momentumStart.Pz()) / ( trkmom_MC[i]* trkmom_MC[i]) , 0.5);
    if ( trkmom_XMC[i] < 0 ) trkstartdoc_XMC[i] = -trkstartdoc_XMC[i];
    if ( trkmom_YMC[i] < 0 ) trkstartdoc_YMC[i] = -trkstartdoc_YMC[i];
    if ( trkmom_ZMC[i] < 0 ) trkstartdoc_ZMC[i] = -trkstartdoc_ZMC[i];
    positionStart.GetXYZT(fMC_startXYZT[i]);
    positionEnd.GetXYZT(fMC_endXYZT[i]);
    trkstartx_MC[i]=fMC_startXYZT[i][0];
    trkstarty_MC[i]=fMC_startXYZT[i][1];
    trkstartz_MC[i]=fMC_startXYZT[i][2];
    trkendx_MC[i]=fMC_endXYZT[i][0];
    trkendy_MC[i]=fMC_endXYZT[i][1];
    trkendz_MC[i]=fMC_endXYZT[i][2];
    tmpVec= positionEnd-positionStart;
    trklen_MC[i]=(positionEnd-positionStart).Rho();
    StartTime_MC[i] = particle->T();                                 // nsec
    double mcdx = StartTime_MC[i] * 1.e-3 * larprop->DriftVelocity();   // cm
    // Calculate the points where this mc particle enters and leaves the
    // fiducial volume, and the length in the fiducial volume.
    TVector3 mcstart;
    TVector3 mcend;
    TVector3 mcstartmom;
    TVector3 mcendmom;
    double plen = length(*particle, mcdx, mcstart, mcend, mcstartmom, mcendmom);
    trklen_cut_MC[i]=plen;
    //trkdedx_MC[i]=trkde/trklen_cut_MC[i];
    //        momentumStart.GetXYZT(fMC_startMomentum[i]);
    //        momentumEnd.GetXYZT(fMC_endMomentum[i]);
    
    // Calculate the points where this mc particle enters and leaves the
    // fiducial volume, and the length in the fiducial volume.
    
    // Get the displacement of this mc particle in the global coordinate system.
    //  TVector3 mcpos = mcstart - pos;
    //TVector3 mcpos = pos -mcstart ;
    // Rotate the momentum and position to the
    // track-local coordinate system.
    //TVector3 mcmoml = rot * mcstartmom;
    //TVector3 mcposl = rot * mcpos;
    trktheta_xz_MC[i] = std::atan2(mcstartmom.X(), mcstartmom.Z());
    trktheta_yz_MC[i] = std::atan2(mcstartmom.Y(), mcstartmom.Z());
    trketa_xy_MC[i] = std::atan2(mcstartmom.X(), mcstartmom.Y());
    trketa_zy_MC[i] = std::atan2(mcstartmom.Z(), mcstartmom.Y());
    
    
    
    trktheta_MC[i]=mcstartmom.Theta();
    trkphi_MC[i]=mcstartmom.Phi();
    
    /*
      trkcolinearity[i] = mcmoml.Z() / mcmoml.Mag();
      double u = mcposl.X();
      double v = mcposl.Y();
      double w = mcposl.Z();
      trkwmatchdisp[i]=w;
      std::cout << "++++++" << std::endl;
      std::cout << "w " << w << std::endl;
      std::cout << "trkcolinearity  " << trkcolinearity[i] << std::endl;
      std::cout << "plen  " << plen << std::endl;
      std::cout << "mcstartmom mag  " << mcstartmom.Mag() << std::endl;
      
      std::cout << "++++++" << std::endl;
      
      double pu = mcmoml.X();
      double pv = mcmoml.Y();
      double pw = mcmoml.Z();
      double dudw = pu / pw;
      double dvdw = pv / pw;
      std::cout << "pu  "<<pu << "pv  " << pv << std::endl;
      std::cout << "pw  "<<pw  << std::endl;
      
      std::cout << "u  "<<u << "v  " << v << std::endl; 
      
      double u0 = u - w * dudw;
      double v0 = v - w * dvdw;
      trkmatchdisp[i]=abs( std::sqrt(u0*u0 + v0*v0));
    */
    
    ++i;
    if (i == kMaxTrack) break;
    //} // particle loop done 
   
    // **********************
    //  Histograms:
    // *********************
    
    // 
    //  Not any more
    // *********************
    
    /*  art::ServiceHandle<cheat::BackTracker> bt2;
	const   sim::ParticleList& pplist=bt2->ParticleList();
	std::vector<const simb::MCParticle*> plist2;
	plist2.reserve(pplist.size());*/
    
    //    art::Handle< std::vector<recob::Track> > trackh;
    //    evt.getByLabel(fTrackModuleLabel, trackh);
    /* 
       if(!trackh.isValid()) continue;
       unsigned int ntrack = trackh->size();
       for(unsigned int i = 0; i < ntrack; ++i) {
       
       art::Ptr<recob::Track> ptrack(trackh, i);
       art::FindMany<recob::Hit>       fmhit(trackListHandle, evt, fTrackModuleLabel);
       std::vector<const recob::Hit*> hits = fmhit.at(i);
       
       //
       // Trick learned from the newest TrackAna
       // Extract hits associated with this track.
       art::FindManyP<recob::Hit> tkhit_find(trackh, evt, fTrackModuleLabel);
       std::vector<art::Ptr<recob::Hit> > trackhits;
       tkhit_find.get(i, trackhits);
       //
       //
       
       const recob::Track& track = *ptrack;
    */
    /*auto pcoll{ptrack};
      art::FindManyP<recob::SpacePoint> fs(pcoll, evt, fTrkSpptAssocModuleLabel);
      auto sppt = fs.at(0);
      art::FindManyP<recob::Hit> fh(sppt, evt, fHitSpptAssocModuleLabel);*/
    ////
    ///              figuring out which TPC
    ///
    ///
    //
    //        auto pcoll { ptrack };
    //art::FindManyP<recob::SpacePoint> fs( pcoll, evt, fTrkSpptAssocModuleLabel);
    //        auto sppt = fs.at(0);//.at(is);
    //        art::FindManyP<recob::Hit> fh( sppt, evt, fHitSpptAssocModuleLabel);
    //	auto hit = fh.at(0).at(0);
    //auto hit = fmhit.at(0).at(0);
    /*	for(int ii=0;ii<hitlist->size())
	{
	
	}
	hit_trkid[i] = fmtk.at(i)[0]->ID();
	int hit_ttpc=-1;
	if(hits.size()!=0)
	{
	geo::WireID tmpWireid=hits.at(0)->WireID();
	hit_ttpc=tmpWireid.TPC;
	}
	else hit_ttpc=1; */
    
    /*  art::Handle< std::vector<recob::Hit> > hitListHandle;
	std::vector<art::Ptr<recob::Hit> > hitlist;
	if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
	art::fill_ptr_vector(hitlist, hitListHandle);
    */
    // Calculate the x offset due to nonzero reconstructed time.
    //double recotime = track.Time() * detprop->SamplingRate();       // nsec
    
    //}
    //}
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
  fTree->Branch("efield",efield,"efield[3]/D");
  fTree->Branch("t0",&t0,"t0/I");
  fTree->Branch("trigtime",trigtime,"trigtime[16]/I");

  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("ntrkhits",ntrkhits,"ntrkhits[ntracks_reco]/I");
  fTree->Branch("trkid",trkid,"trkid[ntracks_reco]/I");  
  fTree->Branch("trkstartx",trkstartx,"trkstartx[ntracks_reco]/D");
  fTree->Branch("trkstarty",trkstarty,"trkstarty[ntracks_reco]/D");
  fTree->Branch("trkstartz",trkstartz,"trkstartz[ntracks_reco]/D");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/D");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/D");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/D");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/D");
  fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/D");
  fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/D");
  fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/D");
  fTree->Branch("trktheta_xz",trktheta_xz,"trktheta_xz[ntracks_reco]/D");
  fTree->Branch("trktheta_yz",trktheta_yz,"trktheta_yz[ntracks_reco]/D");
  fTree->Branch("trketa_xy",trketa_xy,"trketa_xy[ntracks_reco]/D");
  fTree->Branch("trketa_zy",trketa_zy,"trketa_zy[ntracks_reco]/D");
  fTree->Branch("trktheta",trktheta,"trktheta[ntracks_reco]/D");
  fTree->Branch("trkphi",trkphi,"trkphi[ntracks_reco]/D");
  fTree->Branch("trkang",trkang,"trkang[ntracks_reco]/D"); 
  fTree->Branch("trkd2",trkd2,"trkd2[ntracks_reco]/D");
  fTree->Branch("trkcolin",trkcolin,"trkcolin[ntracks_reco]/D");
  fTree->Branch("trkdedx",trkdedx,"trkdedx[ntracks_reco]/D");
  fTree->Branch("trkdedx2",trkdedx2,"trkdedx2[ntracks_reco][3][1000]/D");
  fTree->Branch("trkdqdx",trkdqdx,"trkdqdx[ntracks_reco][3][1000]/D");
  fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][3]/D");
  fTree->Branch("trkpitchHit",trkpitchHit,"trkpitchHit[ntracks_reco][3][1000]/D"); 
  fTree->Branch("trkkinE",trkkinE,"trkkinE[ntracks_reco][3]/D"); 
  fTree->Branch("trkrange",trkrange,"trkrange[ntracks_reco][3]/D"); 
  fTree->Branch("trkTPC",trkTPC,"trkTPC[ntracks_reco][3][1000]/D");
  fTree->Branch("trkplaneid",trkplaneid,"trkplaneid[ntracks_reco][3][1000]/D");
  fTree->Branch("trkresrg",trkresrg,"trkresrg[ntracks_reco][3][1000]/D");
  fTree->Branch("trkPosx",trkPosx,"trkPosx[ntracks_reco][3][1000]/D");
  fTree->Branch("trkPosy",trkPosy,"trkPosy[ntracks_reco][3][1000]/D");
  fTree->Branch("trkPosz",trkPosz,"trkPosz[ntracks_reco][3][1000]/D");
  fTree->Branch("trklen",trklen,"trklen[ntracks_reco]/D");
  fTree->Branch("trklen_L",trklen_L,"trklen_L[ntracks_reco]/D");
  fTree->Branch("trkcolinearity",trkcolinearity,"trkcolinearity[ntracks_reco]/D");
  fTree->Branch("trkmatchdisp",trkmatchdisp,"trkmatchdisp[ntracks_reco]/D");
  fTree->Branch("trkwmatchdisp",trkwmatchdisp,"trkwmatchdisp[ntracks_reco]/D");
  fTree->Branch("trklenratio",trklenratio,"trklenratio[ntracks_reco]/D");
  fTree->Branch("trkdQdxSum",trkdQdxSum,"trkdQdxSum[ntracks_reco]/D");
  fTree->Branch("trkdQdxAverage",trkdQdxAverage,"trkdQdxAverage[ntracks_reco]/D");
  fTree->Branch("trkdEdxSum",trkdEdxSum,"trkdEdxSum[ntracks_reco]/D");
  fTree->Branch("trkdEdxAverage",trkdEdxAverage,"trkdEdxAverage[ntracks_reco]/D");
  fTree->Branch("trkMCTruthT0",trkMCTruthT0,"trkMCTruthT0[ntracks_reco]/D");
  fTree->Branch("trkMCTruthPdG",trkMCTruthPdG,"trkMCTruthPdG[ntracks_reco]/D");
  
  fTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
  fTree->Branch("trkid_MC",trkid_MC,"trkid_MC[nMCParticles]/I");
  fTree->Branch("trkpdg_MC",trkpdg_MC,"trkpdg_MC[nMCParticles]/I");
  fTree->Branch("nCryoHits_MC",&nCryoHits_MC,"nCryoHits_MC[nMCParticles]/I");
  fTree->Branch("trkMother_MC",trkMother_MC,"trkMother_MC[nMCParticles]/I");
  fTree->Branch("trkNumDaughters_MC",trkNumDaughters_MC,"trkNumDaughters_MC[nMCParticles]/I");
  fTree->Branch("trkFirstDaughter_MC",trkFirstDaughter_MC,"trkFirstDaughter_MC[nMCParticles]/I");
  fTree->Branch("trkLastDaughter_MC",trkLastDaughter_MC,"trkLastDaughter_MC[nMCParticles]/I");
  fTree->Branch("trkPrimary_MC",trkPrimary_MC,"trkPrimarys_MC[nMCParticles]/I");
  fTree->Branch("StartTime_MC",StartTime_MC,"StartTime_MC[nMCParticles]/D");
  fTree->Branch("trkstartx_MC",trkstartx_MC,"trkstartx_MC[nMCParticles]/D");
  fTree->Branch("trkstarty_MC",trkstarty_MC,"trkstarty_MC[nMCParticles]/D");
  fTree->Branch("trkstartz_MC",trkstartz_MC,"trkstartz_MC[nMCParticles]/D");
  fTree->Branch("trkendx_MC",trkendx_MC,"trkendx_MC[nMCParticles]/D");
  fTree->Branch("trkendy_MC",trkendy_MC,"trkendy_MC[nMCParticles]/D");
  fTree->Branch("trkendz_MC",trkendz_MC,"trkendz_MC[nMCParticles]/D");
  fTree->Branch("Hits_posx_MC",Hits_posx_MC,"Hits_posx_MC[nMCParticles][1000]/D");
  fTree->Branch("Hits_posy_MC",Hits_posy_MC,"Hits_posy_MC[nMCParticles][1000]/D");
  fTree->Branch("Hits_posz_MC",Hits_posz_MC,"Hits_posz_MC[nMCParticles][1000]/D");
  fTree->Branch("Hits_mom_MC",Hits_mom_MC,"Hits_mom_MC[nMCParticles][1000]/D");
  fTree->Branch("Hits_E_MC",Hits_E_MC,"Hits_E_MC[nMCParticles][1000]/D");
  fTree->Branch("trkenergy_MC",trkenergy_MC,"trkenergy_MC[nMCParticles]/D");
  fTree->Branch("trkmom_MC",trkmom_MC,"trkmom_MC[nMCParticles]/D");
  fTree->Branch("trkmom_XMC",trkmom_XMC,"trkmom_XMC[nMCParticles]/D");
  fTree->Branch("trkmom_YMC",trkmom_YMC,"trkmom_YMC[nMCParticles]/D");
  fTree->Branch("trkmom_ZMC",trkmom_ZMC,"trkmom_ZMC[nMCParticles]/D");
  fTree->Branch("trkstartdoc_XMC",trkstartdoc_XMC,"trkstartdoc_XMC[nMCParticles]/D");
  fTree->Branch("trkstartdoc_YMC",trkstartdoc_YMC,"trkstartdoc_YMC[nMCParticles]/D");
  fTree->Branch("trkstartdoc_ZMC",trkstartdoc_ZMC,"trkstartdoc_ZMC[nMCParticles]/D");
  fTree->Branch("mcpos_x",&mcpos_x,"mcpos_x[nMCParticles]/D");
  fTree->Branch("mcpos_y",&mcpos_y,"mcpos_y[nMCParticles]/D");
  fTree->Branch("mcpos_z",&mcpos_z,"mcpos_z[nMCParticles]/D"); 
  fTree->Branch("mcang_x",&mcang_x,"mcang_x[nMCParticles]/D");
  fTree->Branch("mcang_y",&mcang_y,"mcang_y[nMCParticles]/D");
  fTree->Branch("mcang_z",&mcang_z,"mcang_z[nMCParticles]/D");
  fTree->Branch("trktheta_xz_MC",trktheta_xz_MC,"trktheta_xz_MC[nMCParticles]/D");
  fTree->Branch("trktheta_yz_MC",trktheta_yz_MC,"trktheta_yz_MC[nMCParticles]/D");
  fTree->Branch("trktheta_MC",trktheta_MC,"trktheta_MC[nMCParticles]/D");
  fTree->Branch("trkphi_MC",trkphi_MC,"trkphi_MC[nMCParticles]/D");
  fTree->Branch("trketa_xy_MC",trketa_xy_MC,"trketa_xy_MC[nMCParticles]/D");
  fTree->Branch("trketa_zy_MC",trketa_zy_MC,"trketa_zy_MC[nMCParticles]/D");
  fTree->Branch("trkdedx_MC",trkdedx_MC,"trkdedx_MC[nMCParticles]/D");
  fTree->Branch("trkdq_MC",trkdq_MC,"trkdq_MC[nMCParticles]/D");
  fTree->Branch("trklen_MC",trklen_MC,"trklen_MC[nMCParticles]/D");
  fTree->Branch("trklen_cut_MC",trklen_cut_MC,"trklen_cut_MC[nMCParticles]/D");
 
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("nclust",&nclust,"nclust/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I");
  fTree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/D");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/D");
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I");

  //  art::ServiceHandle<sim::LArG4Parameters> larParameters;
  //  fElectronsToGeV = 1./larParameters->GeVToElectrons();


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
    trklen_MC[i] = -99999;
    trklen_cut_MC[i] = -99999;
    trkenergy_MC[i] = -99999;
    trkmom_MC[i] = -99999;
    trkmom_XMC[i] = -99999;
    trkmom_YMC[i] = -99999;
    trkmom_ZMC[i] = -99999;
    trkstartdoc_XMC[i] = -99999;
    trkstartdoc_YMC[i] = -99999;
    trkstartdoc_ZMC[i] = -99999;
    trkid_MC[i] = -99999;
    trkpdg_MC[i] = -99999;
    trkMother_MC[i] = -99999;
    trkNumDaughters_MC[i] = -99999;
    trkFirstDaughter_MC[i] = -99999;
    trkLastDaughter_MC[i] = -99999;
    trkPrimary_MC[i] = 0;
    trkd2[i] = -99999;
    trkcolin[i] = -99999;
    trktheta_xz_MC[i] = -99999;
    trktheta_yz_MC[i] = -99999;
    trktheta_MC[i] = -99999;
    trkphi_MC[i] = -99999;
    trkdedx[i] = -99999;
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
    nCryoHits_MC[i]=0;
    StartTime_MC[i]=0;
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
    trkdedx_MC[i] = -99999;
    trkdq_MC[i] = -99999;

    trkstartdcosx[i] = -99999;
    trkstartdcosy[i] = -99999;
    trkstartdcosz[i] = -99999;
    trklen[i] = -99999;
    trklen_L[i] = -99999;
    trkid[i] = -99999;
    trkang[i] = -99999;
    trkcolinearity[i] = -99999;
    trkmatchdisp[i] = -99999;
    trkwmatchdisp[i] = -99999;
    trklenratio[i] = -99999;
    trkMCTruthT0[i] = -99999;
    trkMCTruthPdG[i] = -99999;

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
      Hits_mom_MC[i][j] = -99999;
      Hits_E_MC[i][j] = -99999;
    }
    for (int j = 0; j<3; ++j){
      trkpitch[i][j] = -99999;
    }
  }
  nhits = -99999;

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

}

void AnaTree::AnaTree::endJob()
{
}


DEFINE_ART_MODULE(AnaTree::AnaTree)
