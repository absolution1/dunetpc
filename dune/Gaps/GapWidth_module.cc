////////////////////////////////////////////////////////////////////////
// Class:      GapWidth
// Module Type: analyzer
// File:        GapWidth_module.cc
// Author:      Tristan Blackburn (Sussex) t.blackburn@sussex.ac.uk
//
// Module skeleton is AnaTree - K. Warburton (Sheffield)                   
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
//Commenting this out since it doesn't appear to be used.
//#include "larsim/MCCheater/BackTrackerService.h"
//#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TLorentzVector.h"
#include "TH2F.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TPaveStats.h"
#include "TROOT.h"

//standard library includes
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits> // std::numeric_limits<>
#include <vector>
#include <algorithm>
#include <array>
#include <iterator>
#include <cmath>

const int kMaxTrack      =  1000; //maximum number of tracks
const int kMaxHits       = 10000; //maximum number of hits
//const int kMaxClust      = 10000; //maximum number of clusters // unused
const int kMaxTrackHits  =  1000; //maximum number of space points
const int kMaxEvent      = 10000; //maximum number of events

//Some loop parameters
	    
double track1x = 0;
double track1y = 0;
double track1z = 0;
  
double track2x = 0;
double track2y = 0;
double track2z = 0;
  
double track1grad = 0; 
double track2grad = 0;

double  minx = 0;
double  miny = 0;
double  minz = 0;

int posminx; 
int posminy; 
int posminz;

double zdiff = 0;
double xdiff = 0;
double ydiff = 0;

double z = 0;
double y = 0;
double translation = 0;


namespace GapWidth {
  class GapWidth;

  double xextrap(double z) {
    if (track1grad <= 0) return fabs((track1x - (track1grad)*(minz + z)) - track2x);
    else return fabs((track1x + (track1grad)*(minz + z)) - track2x);
  }
  
  void minuitFunctionx(int& nDim, double* gout, double& result, double par[], int flg){
    result = xextrap(par[0]);
  }

  double zextrap(double y) {
    if (track1grad <= 0) return fabs((track1z - (track1grad)*(miny + y)) - track2z);
    else return fabs((track1z + (track1grad)*(miny + y)) - track2z);
  }

  void minuitFunctionz(int& nDim, double* gout, double& result, double par[], int flg){
    result = zextrap(par[0]);
  }
}


class GapWidth::GapWidth : public art::EDAnalyzer {
public:
  explicit GapWidth(fhicl::ParameterSet const & p);
  virtual ~GapWidth();

  void analyze(art::Event const & e) override;
  void beginJob() override;
  void endJob() override;


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

  int ntracks_reco;         //number of reconstructed tracks
  double trkstartx[kMaxTrack];
  double trkstarty[kMaxTrack];
  double trkstartz[kMaxTrack];
  double trkendx[kMaxTrack];
  double trkendy[kMaxTrack];
  double trkendz[kMaxTrack];

  double trkd2[kMaxTrack];
  double trklen[kMaxTrack];
  double trklen_L[kMaxTrack];
  int    trkid[kMaxTrack];
  
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
  int    nhits2;
  int    nclust;

  int    hit_tpc[kMaxHits];
  int    hit_plane[kMaxHits];
  int    hit_wire[kMaxHits];
  int    hit_channel[kMaxHits];
  double hit_peakT[kMaxHits];
  double hit_charge[kMaxHits];
  double hit_ph[kMaxHits];
  int    hit_trkid[kMaxHits];

  double trkMCTruthT0[kMaxTrack];

  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fClusterModuleLabel;
  std::string fTrkSpptAssocModuleLabel;
  std::string fHitSpptAssocModuleLabel;
  std::string fSimulationProducerLabel; 
  std::string fCalorimetryModuleLabel; 
  std::string fMCTruthT0ModuleLabel;

  //double fElectronsToGeV; // conversion factor // unused
  art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

  double tracklengthXZ[kMaxTrack];
  double tracklengthYZ[kMaxTrack];
  double trackgradientXZ[kMaxTrack];
  double trackgradientYZ[kMaxTrack];

  //vectors of possible minimised distances from possible start and end points between two adjacent TPC tracks
  std::vector<double> possiblex; 
  std::vector<double> possibley;
  std::vector<double> possiblez;

  //Arrays for computing the track TPC using midpoints
  double trackmidx[kMaxTrack];
  double trackmidy[kMaxTrack];
  double trackmidz[kMaxTrack];
  int     tracktpc[kMaxTrack];

  //Array to store output
  double output[kMaxEvent][5][5];

  // double xyz[3];                                                                                                                                   
  // double abc[3];                                                                 
  // int chan;
  // int cryo = fGeometry->Ncryostats();                                                             
  // for (int c = 0; c<cryo; ++c){                   
  //   int tpc =fGeometry->NTPC(c);                                               
  //   for (int t=0; t<tpc; ++t){                                            
  //     int Nplanes=fGeometry->Nplanes(t,c);                                      
  //     for (int p=0; p<Nplanes; ++p) {                                        
  // 	int Nwires = fGeometry->Nwires(p,t,c);                                  
  // 	for (int w=0; w<Nwires; ++w){
  // 	  fGeometry->WireEndPoints(c,t,p,w,xyz,abc);
  // 	  chan=fGeometry->PlaneWireToChannel(p,w,t,c);  
	  
	  
  // 	  if (t % 2 != 0 && p == 2){
  //           std::cout << "FLAG " << chan << " " << c << " " << t << " " << p << " " << w << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] <<  " " << abc[0] << " " << abc[1] << " " << abc[2] << std::endl;                                         
  // 	  }
  // 	  //c=cryo, t=tpc, p=plane, w=wire, start (x,y,z), end (x,y,z)
 
  // 	}
  //     }
  //   }  
  // }

  //hard code (need to soft code) some detector YZ detector boundaries using geometry v5
  double upper1 = 113.142; //Upper & lower refer to Y
  double lower1 = -82.308;
  double left1  = 0.29937; //Left and right refer to Z
  double right1 = 50.1445;

  double upper3 = -1.44705;
  double lower3 = -84.4553;
  double left3  =  52.6724;
  double right3 =  102.518;
  
  double upper5 = 113.142;
  double lower5 = 1.46205;
  double left5  = 52.2234;
  double right5 = 102.068;

  double upper7 = 113.142;
  double lower7 = -82.308;
  double left7  = 104.147;
  double right7 = 153.992;

};

GapWidth::GapWidth::GapWidth(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fHitsModuleLabel         ( pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel        ( pset.get< std::string >("TrackModuleLabel"))
  , fClusterModuleLabel      ( pset.get< std::string >("ClusterModuleLabel"))
  , fTrkSpptAssocModuleLabel ( pset.get< std::string >("TrkSpptAssocModuleLabel"))
  , fHitSpptAssocModuleLabel ( pset.get< std::string >("HitSpptAssocModuleLabel"))
  , fSimulationProducerLabel ( pset.get< std::string >("SimulationLabel"))
  , fCalorimetryModuleLabel  ( pset.get< std::string >("CalorimetryModuleLabel"))
  , fMCTruthT0ModuleLabel    ( pset.get< std::string >("MCTruthT0ModuleLabel"))
{
}

GapWidth::GapWidth::~GapWidth()
{
  // Clean up dynamic memory and other resources here.
}

void GapWidth::GapWidth::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  ResetVars();

  art::ServiceHandle<geo::Geometry> geom;
  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
//  art::ServiceHandle<cheat::BackTrackerService> bktrk;
//  Commenting this out since it does not appear to be used anywhere.

  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory topdir = tfs->mkdir("trkgaps", "Gap histograms");
  art::TFileDirectory gap1 = topdir.mkdir("Gap 1");
  art::TFileDirectory gap2 = topdir.mkdir("Gap 2");
  art::TFileDirectory gap3 = topdir.mkdir("Gap 3");
  art::TFileDirectory gap4 = topdir.mkdir("Gap 4");
  art::TFileDirectory gap5 = topdir.mkdir("Gap 5");
  
  TH1D * gapit1[kMaxEvent];
  TH1D * gapdif1[kMaxEvent];

  TH1D * gapit2[kMaxEvent];
  TH1D * gapdif2[kMaxEvent];

  TH1D * gapit3[kMaxEvent];
  TH1D * gapdif3[kMaxEvent];

  TH1D * gapit4[kMaxEvent];
  TH1D * gapdif4[kMaxEvent];

  TH1D * gapit5[kMaxEvent];
  TH1D * gapdif5[kMaxEvent];

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

  art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
  evt.getByLabel(fSimulationProducerLabel, simChannelHandle);

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

  TVector3 larStart;
  TVector3 larEnd;
  
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
  //  NOW DO THE TRACK LIST STUFF 
  // ------------------------------------
  
  if ( trackListHandle.isValid() ) {
    art::FindManyP<recob::SpacePoint> fmsp  (trackListHandle, evt, fTrackModuleLabel);
    art::FindManyP<recob::Hit>        fmht  (trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal (trackListHandle, evt, fCalorimetryModuleLabel);
    art::FindMany<anab::T0>           fmt0 (trackListHandle, evt,  fMCTruthT0ModuleLabel);
     
    for(int i=0; i<std::min(int(tracklist.size()),kMaxTrack);++i){


      //***************************************************
      //   T0 stuff - So can correct X Start/End positions
      //***************************************************
      double TickT0 = 0;
      if ( fmt0.isValid() ) {
	std::vector<const anab::T0*> T0s = fmt0.at(i);
	for (size_t t0size =0; t0size < T0s.size(); t0size++) {
	  trkMCTruthT0[i]  = T0s[t0size]->Time();
	  TickT0 = trkMCTruthT0[i] / detprop->SamplingRate();
	} // T0 size
      } // T0 valid
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
      recob::Track::Point_t trackStart, trackEnd;
      std::tie(trackStart, trackEnd) = tracklist[i]->Extent(); 
      larStart = tracklist[i]->VertexDirection();
      larEnd = tracklist[i]->EndDirection();

      trkstartx[i]      = trackStart.X() - detprop->ConvertTicksToX( TickT0, allHits[Hit_Size-1]->WireID().Plane, allHits[Hit_Size-1]->WireID().TPC, allHits[Hit_Size-1]->WireID().Cryostat ); // Correct X, last entry is first 'hit'
      trkstarty[i]      = trackStart.Y();
      trkstartz[i]      = trackStart.Z();
      trkendx[i]        = trackEnd.X() - detprop->ConvertTicksToX( TickT0, allHits[0]->WireID().Plane, allHits[0]->WireID().TPC, allHits[0]->WireID().Cryostat ); // Correct X, first entry is last 'hit'
      trkendy[i]        = trackEnd.Y();
      trkendz[i]        = trackEnd.Z();
      trkstartdcosx[i]  = larStart[0];
      trkstartdcosy[i]  = larStart[1];
      trkstartdcosz[i]  = larStart[2];
      trkenddcosx[i]    = larEnd[0];
      trkenddcosy[i]    = larEnd[1];
      trkenddcosz[i]    = larEnd[2];
      TLorentzVector v1(trackStart.X(),trackStart.Y(),trackStart.Z(),0);
      TLorentzVector v2(trackEnd.X(),trackEnd.Y(),trackEnd.Z(),0);
      trklen[i]=(v2-v1).Rho();

      art::Ptr<recob::Track> ptrack(trackh, i);
      const recob::Track& track = *ptrack;
      trklen_L[i]=track.Length();

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
      TVector3 V1(trackStart.X(),trackStart.Y(),trackStart.Z());
      TVector3 V2(trackEnd.X(),trackEnd.Y(),trackEnd.Z());
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
	  
	  if (spts[j]->XYZ()[0] < -1000) continue;
	  else trkx[i][j] = spts[j]->XYZ()[0];
	  if (spts[j]->XYZ()[1] < -1000) continue;
	  else trky[i][j] = spts[j]->XYZ()[1];
	  if (spts[j]->XYZ()[2] < -1000) continue;
	  else trky[i][j] = spts[j]->XYZ()[2];

	  std::cout << "X, Y, Z: " << spts[j]->XYZ()[0] << ", " << spts[j]->XYZ()[1] << ", " << spts[j]->XYZ()[2] << std::endl;
	 	 
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
	    trkdedx2[i][jj][k]     = (calos[jj] -> dEdx())[k];
	    trkdqdx[i][jj][k]      = (calos[jj] -> dQdx())[k];
	    trkpitchHit[i][jj][k]  = (calos[jj] -> TrkPitchVec())[k];
	    trkresrg[i][jj][k]     = (calos[jj] -> ResidualRange())[k];
	    trkTPC[i][jj][k]       = (calos[jj]->PlaneID()).TPC;
	    trkplaneid[i][jj][k]   = (calos[jj]->PlaneID()).Plane;
	    
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
	    trkpitch[i][j] = lar::util::TrackPitchInView(*(tracklist[i]), geo::kU);
	  else if (j==1)
	    trkpitch[i][j] = lar::util::TrackPitchInView(*(tracklist[i]), geo::kV);
	  else if (j==2)
	    trkpitch[i][j] = lar::util::TrackPitchInView(*(tracklist[i]), geo::kZ);
	}
	catch( cet::exception &e) {
	  mf::LogWarning("GapWidth")<<"caught exception "<<e<<"\n setting pitch to 0";
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
      if (fmtk.at(i).size()!=0 && fmtk.at(i)[0]->ID() <= 7){
       	std::cout << "Hit TrackID: " << fmtk.at(i)[0]->ID() << std::endl;
      	hit_trkid[i] = fmtk.at(i)[0]->ID();
      }
      if ( i == kMaxHits ) break;
    }
  }

  gapit1[event]  = gap1.make<TH1D>(Form("gapit1, from event %d",event),Form("Gap Width Against Translation, Event %d",event), 101, -5.05, 5.05);
  gapdif1[event] = gap1.make<TH1D>(Form("gapdif1 from event %d",event),Form("Gap Difference from Nominal, Event %d",event), 101, -5.05, 5.05);

  gapit2[event]  = gap2.make<TH1D>(Form("gapit2, from event %d",event),Form("Gap Width Against Translation, Event %d",event), 101, -5.05, 5.05);
  gapdif2[event] = gap2.make<TH1D>(Form("gapdif2 from event %d",event),Form("Gap Difference from Nominal, Event %d",event), 101, -5.05, 5.05);

  gapit3[event]  = gap3.make<TH1D>(Form("gapit3, from event %d",event),Form("Gap Width Against Translation, Event %d",event), 101, -5.05, 5.05);
  gapdif3[event] = gap3.make<TH1D>(Form("gapdif3 from event %d",event),Form("Gap Difference from Nominal, Event %d",event), 101, -5.05, 5.05);

  gapit4[event]  = gap4.make<TH1D>(Form("gapit4, from event %d",event),Form("Gap Width Against Translation, Event %d",event), 101, -5.05, 5.05);
  gapdif4[event] = gap4.make<TH1D>(Form("gapdif4 from event %d",event),Form("Gap Difference from Nominal, Event %d",event), 101, -5.05, 5.05);

  gapit5[event]  = gap5.make<TH1D>(Form("gapit5, from event %d",event),Form("Gap Width Against Translation, Event %d",event), 101, -5.05, 5.05);
  gapdif5[event] = gap5.make<TH1D>(Form("gapdif5 from event %d",event),Form("Gap Difference from Nominal, Event %d",event), 101, -5.05, 5.05);
       
    //First loop assigns tracks TPCs and limits tracks such that they don't extend beyond TPC's they occur in
    for (int i = 0; i < ntracks_reco; i++) { 

      //First, find track TPC's using mid points

      trackmidx[i] = (trkendx[i] + trkstartx[i])/2.;
      trackmidy[i] = (trkendy[i] + trkstarty[i])/2.;
      trackmidz[i] = (trkendz[i] + trkstartz[i])/2.;

      if (trackmidx[i] > 0.) {
	if (trackmidz[i] > 0. && trackmidz[i] < 51.) tracktpc[i]=1;
	if (trackmidz[i] > 51. && trackmidz[i] < 103. && trackmidy[i] < 0.) tracktpc[i]=3;
	if (trackmidz[i] > 51. && trackmidz[i] < 103. && trackmidy[i] > 0.) tracktpc[i]=5;
	if (trackmidz[i] > 103. && trackmidz[i] < 155.) tracktpc[i]=7;
      }

      //Floor boundaries to coincide with wire ends
      
      if (tracktpc[i] == 1) {
	if (trkendy[i]   > upper1) trkendy[i] = upper1;
	if (trkstarty[i] > upper1) trkstarty[i] = upper1;
	if (trkendy[i]   < lower1) trkendy[i] = lower1;
	if (trkstarty[i] < lower1) trkstarty[i] = lower1;
	if (trkendz[i]   < left1 ) trkendz[i] = left1;
	if (trkstartz[i] < left1 ) trkstartz[i] = left1;
	if (trkendz[i]   > right1) trkendz[i] = right1;
	if (trkstartz[i] > right1) trkstartz[i] = right1;	  
      }

      if (tracktpc[i] == 3) {
	if (trkendy[i]   > upper3) trkendy[i] = upper3;
	if (trkstarty[i] > upper3) trkstarty[i] = upper3;
	if (trkendy[i]   < lower3) trkendy[i] = lower3;
	if (trkstarty[i] < lower3) trkstarty[i] = lower3;
	if (trkendz[i]   < left3 ) trkendz[i] = left3;
	if (trkstartz[i] < left3 ) trkstartz[i] = left3;
	if (trkendz[i]   > right3) trkendz[i] = right3;
	if (trkstartz[i] > right3) trkstartz[i] = right3;	
      }

      if (tracktpc[i] == 5) {
	if (trkendy[i]   > upper5) trkendy[i] = upper5;
	if (trkstarty[i] > upper5) trkstarty[i] = upper5;
	if (trkendy[i]   < lower5) trkendy[i] = lower5;
	if (trkstarty[i] < lower5) trkstarty[i] = lower5;
	if (trkendz[i]   < left5 ) trkendz[i] = left5;
	if (trkstartz[i] < left5 ) trkstartz[i] = left5;
	if (trkendz[i]   > right5) trkendz[i] = right5;
	if (trkstartz[i] > right5) trkstartz[i] = right5;	
      }

      if (tracktpc[i] == 7) {
	if (trkendy[i]   > upper7) trkendy[i] = upper7;
	if (trkstarty[i] > upper7) trkstarty[i] = upper7;
	if (trkendy[i]   < lower7) trkendy[i] = lower7;
	if (trkstarty[i] < lower7) trkstarty[i] = lower7;
	if (trkendz[i]   < left7 ) trkendz[i] = left7;
	if (trkstartz[i] < left7 ) trkstartz[i] = left7;
	if (trkendz[i]   > right7) trkendz[i] = right7;
	if (trkstartz[i] > right7) trkstartz[i] = right7;	
      }

      tracklengthXZ[i] = TMath::Power((TMath::Power(std::fabs(trkendz[i] - trkstartz[i]), 2) + TMath::Power(std::fabs(trkendx[i] - trkstartx[i]), 2)), 0.5);
      tracklengthYZ[i] = TMath::Power((TMath::Power(std::fabs(trkendz[i] - trkstartz[i]), 2) + TMath::Power(std::fabs(trkendy[i] - trkstarty[i]), 2)), 0.5);
      
      if (trkendz[i]>trkstartz[i]) xdiff = (trkendx[i] - trkstartx[i]);
      else                         xdiff = (trkstartx[i] - trkendx[i]);
      
      if (trkendz[i]>trkstartz[i]) ydiff = (trkendy[i] - trkstarty[i]);
      else                         ydiff = (trkstarty[i] - trkendy[i]);
      
      if (trkendz[i]>trkstartz[i]) zdiff = fabs(trkendz[i] - trkstartz[i]);
      else                         zdiff = fabs(trkstartz[i] - trkendz[i]);
      
      trackgradientXZ[i] = zdiff/xdiff;
      trackgradientYZ[i] = zdiff/ydiff;
    }

    for (int translationsteps = 0; translationsteps < 101; translationsteps ++){
      translation = (double(translationsteps)-50)/10.;

      //Start of loop for track comparison
      for (int i = 0; i < ntracks_reco; i++){
	
	//need comparison track
	for (int k = 0; k < ntracks_reco; k++){
 
	  possiblex.clear();
	  possibley.clear();
	  possiblez.clear();

	  track1x = 0;
	  track1y = 0;
	  track1z = 0;
	  track2x = 0;
	  track2y = 0;
	  track2z = 0;
	  track1grad = 0;
	  track2grad = 0;
	  posminx = 0;
	  posminy = 0;
	  posminz = 0;

	  //force tracks to be adjacent in tpcs with similar gradients
	  if ((tracktpc[i] != tracktpc[k]) && (abs(tracktpc[i] - tracktpc[k]) <= 4) && (fabs(atan(trackgradientXZ[k]))*TMath::Pi()/180 - (atan(trackgradientXZ[i]))*TMath::Pi()/180 < 0.1) && (tracktpc[i] < tracktpc[k])) { 

	    possiblex.push_back(fabs(trkendx[i] - trkendx[k]));
	    possiblex.push_back(fabs(trkendx[i] - trkstartx[k]));
	    possiblex.push_back(fabs(trkstartx[i] - trkendx[k]));
	    possiblex.push_back(fabs(trkstartx[i] - trkstartx[k]));

	    possibley.push_back(fabs(trkendy[i] - trkendy[k]));
	    possibley.push_back(fabs(trkendy[i] - trkstarty[k]));
	    possibley.push_back(fabs(trkstarty[i] - trkendy[k]));
	    possibley.push_back(fabs(trkstarty[i] - trkstarty[k]));

	    possiblez.push_back(fabs(trkendz[i] - trkendz[k]));
	    possiblez.push_back(fabs(trkendz[i] - trkstartz[k]));
	    possiblez.push_back(fabs(trkstartz[i] - trkendz[k]));
	    possiblez.push_back(fabs(trkstartz[i] - trkstartz[k]));	  	

	    minx = * std::min_element(possiblex.begin(), possiblex.end());
	    miny = * std::min_element(possibley.begin(), possibley.end());
	    minz = * std::min_element(possiblez.begin(), possiblez.end());
	 
	    for (size_t g = 0; g < possiblex.size(); g++) {
	      if (possiblex[g] == minx) {
		posminx = g;
	      }	    
	    }
	    for (size_t g = 0; g < possibley.size(); g++) {
	      if (possibley[g] == miny) {
		posminy = g;
	      }
	    }
	    for (size_t g = 0; g < possiblez.size(); g++) {
	      
	      if (possiblez[g] == minz) { 
		posminz = g;
	      }
	    }  


	    //GAP1
	    if ((tracktpc[i] == 1 && tracktpc[k] == 5) || (tracktpc[i] == 5 && tracktpc[k] == 1)) {	  
	      if (minz <= 6.) {
		//std::cout << "From the comparison between tracks with IDs: " << trkid[i] << " and " << trkid[k];	     
		//std::cout << ", the size of Gap 1 is: " << minz << " at position " << posminz << std::endl; 


		if (posminz <= 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkstartz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkstartx[k] + translation;
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkendx[i] + translation;
		    track2y =  trkendy[i];
		    track2z =  trkendz[i];
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientXZ[k];
		  }	
		}

		if (posminz <= 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkendz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkendx[k] + translation;
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkendx[i] + translation;
		    track2y =  trkendy[i];
		    track2z =  trkendz[i];
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
	      
		if (posminz  > 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
   	
		  if (trkstartz[i] < trkstartz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkstartx[k] + translation;
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkstartx[i] + translation;
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i];
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
		if (posminz  > 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkstartz[i] < trkendz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkendx[k] + translation;
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkstartx[i] + translation;
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i];
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
	      
	    	double error = 0.01;
		TMinuit* mingen = new TMinuit(1);
		mingen->SetFCN(minuitFunctionx);
		mingen->SetPrintLevel(-1);
		mingen->DefineParameter(0,"gap extension", -20, 0.001, -100, 100);
		mingen->Migrad();
		mingen->GetParameter(0, z, error);
	
		output[event][0][0] = double(trkid[i]);
		output[event][0][1] = double(trkid[k]);
		output[event][0][2] = double(tracktpc[i]);
		output[event][0][3] = double(tracktpc[k]);
		output[event][0][4] = (minz + z);

	      }
	    } // End Gap 1

	    //GAP2
	    if ((tracktpc[i] == 5 && tracktpc[k] == 7) || (tracktpc[i] == 7 && tracktpc[k] == 5)) {	  
	      if (minz <= 6.) {
  
		if (posminz <= 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkstartz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkstartx[k] + translation;
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkendx[i] + translation;
		    track2y =  trkendy[i];
		    track2z =  trkendz[i];
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientXZ[k];
		  }	
		}

		if (posminz <= 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkendz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkendx[k] + translation;
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkendx[i] + translation;
		    track2y =  trkendy[i];
		    track2z =  trkendz[i];
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
	      
		if (posminz  > 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
   	
		  if (trkstartz[i] < trkstartz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkstartx[k] + translation;
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkstartx[i] + translation;
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i];
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
		if (posminz  > 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkstartz[i] < trkendz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkendx[k] + translation;
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkstartx[i] + translation;
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i];
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
	      
	    	double error = 0.01;
		TMinuit* mingen = new TMinuit(1);
		mingen->SetFCN(minuitFunctionx);
		mingen->SetPrintLevel(-1);
		mingen->DefineParameter(0,"gap extension", -20, 0.001, -100, 100);
		mingen->Migrad();
		mingen->GetParameter(0, z, error);
	
		output[event][1][0] = double(trkid[i]);
		output[event][1][1] = double(trkid[k]);
		output[event][1][2] = double(tracktpc[i]);
		output[event][1][3] = double(tracktpc[k]);
		output[event][1][4] = (minz + z);

	      }
	    } // End Gap 2

	    //GAP3
	    if ((tracktpc[i] == 1 && tracktpc[k] == 3) || (tracktpc[i] == 3 && tracktpc[k] == 1)) {	  
	      if (minz <= 6.) {

		if (posminz <= 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkstartz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkstartx[k] + translation;
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkendx[i] + translation;
		    track2y =  trkendy[i];
		    track2z =  trkendz[i];
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientXZ[k];
		  }	
		}

		if (posminz <= 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkendz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkendx[k] + translation;
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkendx[i] + translation;
		    track2y =  trkendy[i];
		    track2z =  trkendz[i];
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
	      
		if (posminz  > 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
   	
		  if (trkstartz[i] < trkstartz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkstartx[k] + translation;
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkstartx[i] + translation;
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i];
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
		if (posminz  > 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkstartz[i] < trkendz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkendx[k] + translation;
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkstartx[i] + translation;
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i];
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
	      
	    	double error = 0.01;
		TMinuit* mingen = new TMinuit(1);
		mingen->SetFCN(minuitFunctionx);
		mingen->SetPrintLevel(-1);
		mingen->DefineParameter(0,"gap extension", -20, 0.001, -100, 100);
		mingen->Migrad();
		mingen->GetParameter(0, z, error);
	
		output[event][2][0] = double(trkid[i]);
		output[event][2][1] = double(trkid[k]);
		output[event][2][2] = double(tracktpc[i]);
		output[event][2][3] = double(tracktpc[k]);
		output[event][2][4] = (minz + z);

	      }
	    } // End Gap 3

	    //GAP4
	    if ((tracktpc[i] == 3 && tracktpc[k] == 7) || (tracktpc[i] == 7 && tracktpc[k] == 3)) {	  
	      if (minz <= 6.) {

		if (posminz <= 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkstartz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkstartx[k] + translation;
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkendx[i] + translation;
		    track2y =  trkendy[i];
		    track2z =  trkendz[i];
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientXZ[k];
		  }	
		}

		if (posminz <= 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkendz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkendx[k] + translation;
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkendx[i] + translation;
		    track2y =  trkendy[i];
		    track2z =  trkendz[i];
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
	      
		if (posminz  > 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
   	
		  if (trkstartz[i] < trkstartz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkstartx[k] + translation;
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkstartx[i] + translation;
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i];
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
		if (posminz  > 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkstartz[i] < trkendz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkendx[k] + translation;
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k];
		    track1grad = trackgradientXZ[i];
		  }
		  else {
		    track2x =  trkstartx[i] + translation;
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i];
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientXZ[k];
		  }
		}
	      
	    	double error = 0.01;
		TMinuit* mingen = new TMinuit(1);
		mingen->SetFCN(minuitFunctionx);
		mingen->SetPrintLevel(-1);
		mingen->DefineParameter(0,"gap extension", -20, 0.001, -100, 100);
		mingen->Migrad();
		mingen->GetParameter(0, z, error);
	
		output[event][3][0] = double(trkid[i]);
		output[event][3][1] = double(trkid[k]);
		output[event][3][2] = double(tracktpc[i]);
		output[event][3][3] = double(tracktpc[k]);
		output[event][3][4] = (minz + z);

	      }
	    } // End Gap 4

	    //GAP 5
	    if ((tracktpc[i] == 3 && tracktpc[k] == 5) || (tracktpc[i] == 5 && tracktpc[k] == 3)) {
	      if (miny <= 6.) {

		if (posminz <= 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkstartz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkstartx[k];
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k] + translation;
		    track1grad = trackgradientYZ[i];
		    track2grad = trackgradientYZ[k];
		  }
		  else {
		    track2x =  trkendx[i];
		    track2y =  trkendy[i];
		    track2z =  trkendz[i] + translation;
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientYZ[k];
		    track2grad = trackgradientYZ[i];
		  }	
		}

		if (posminz <= 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching end of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkendz[i] < trkendz[k]){
		    track1x =  trkendx[i];
		    track1y =  trkendy[i];
		    track1z =  trkendz[i];
		    track2x =  trkendx[k];
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k] + translation;
		    track1grad = trackgradientYZ[i];
		    track2grad = trackgradientYZ[k];
		  }
		  else {
		    track2x =  trkendx[i];
		    track2y =  trkendy[i];
		    track2z =  trkendz[i] + translation;
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientYZ[k];
		    track2grad = trackgradientYZ[i];
		  }
		}
	      
		if (posminz  > 1 && posminz % 2 ==1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the start of track: " << trkid[k] << std::endl;
   	
		  if (trkstartz[i] < trkstartz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkstartx[k];
		    track2y =  trkstarty[k];		  
		    track2z =  trkstartz[k] + translation;
		    track1grad = trackgradientYZ[i];
		    track2grad = trackgradientYZ[k];
		  }
		  else {
		    track2x =  trkstartx[i];
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i] + translation;
		    track1x =  trkstartx[k];
		    track1y =  trkstarty[k];		  
		    track1z =  trkstartz[k];
		    track1grad = trackgradientYZ[k];
		    track2grad = trackgradientYZ[i];
		  }
		}
		if (posminz  > 1 && posminz % 2 !=1){
		  //std::cout << "Z minimised matching start of track: " << trkid[i] << " and the end of track: " << trkid[k] << std::endl;
		
		  if (trkstartz[i] < trkendz[k]){
		    track1x =  trkstartx[i];
		    track1y =  trkstarty[i];
		    track1z =  trkstartz[i];
		    track2x =  trkendx[k];
		    track2y =  trkendy[k];		  
		    track2z =  trkendz[k] + translation;
		    track1grad = trackgradientYZ[i];
		    track2grad = trackgradientYZ[k];
		  }
		  else {
		    track2x =  trkstartx[i];
		    track2y =  trkstarty[i];
		    track2z =  trkstartz[i] + translation;
		    track1x =  trkendx[k];
		    track1y =  trkendy[k];		  
		    track1z =  trkendz[k];
		    track1grad = trackgradientYZ[k];
		    track2grad = trackgradientYZ[i];
		  }
		}	       

		double error = 0.01;
		TMinuit* min5 = new TMinuit(1);
		min5->SetFCN(minuitFunctionz);
		min5->SetPrintLevel(-1);
		min5->DefineParameter(0,"gap extension", -20, 0.001, -100, 100);
		min5->Migrad();
		min5->GetParameter(0, y, error);

		output[event][4][0] = double(trkid[i]);
		output[event][4][1] = double(trkid[k]);
		output[event][4][2] = double(tracktpc[i]);
		output[event][4][3] = double(tracktpc[k]);
		output[event][4][4] = (miny + y);		
		
	      }
	    }//End Gap 5 




	  }//End of check that TPCs are adjacent and that track ID's are not identical
	}//End of loop over tracks[k]
      } //End of loop over tracks[i]           
      
      gapit1[event]->SetBinContent(translationsteps+1, output[event][0][4]);
      gapdif1[event]->SetBinContent(translationsteps+1, fabs(output[event][0][4] - 2.0789));

      gapit2[event]->SetBinContent(translationsteps+1, output[event][1][4]);
      gapdif2[event]->SetBinContent(translationsteps+1, fabs(output[event][1][4] - 2.079));

      gapit3[event]->SetBinContent(translationsteps+1, output[event][2][4]);
      gapdif3[event]->SetBinContent(translationsteps+1, fabs(output[event][2][4] - 2.5279));

      gapit4[event]->SetBinContent(translationsteps+1, output[event][3][4]);
      gapdif4[event]->SetBinContent(translationsteps+1, fabs(output[event][3][4] - 1.629));

      gapit5[event]->SetBinContent(translationsteps+1, output[event][4][4]);
      gapdif5[event]->SetBinContent(translationsteps+1, fabs(output[event][4][4] - 2.9091));
      
    }
    
  // -------------------------------------------
  // FINALLLY Fill Tree:
  // -------------------------------------------
  fTree->Fill();
  
}


void GapWidth::GapWidth::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs2;
  fTree = tfs2->make<TTree>("gapwidth","analysis tree");

  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("efield",efield,"efield[3]/D");
  fTree->Branch("t0",&t0,"t0/I");

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
  
  //Array branches start
  fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/D");
  fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/D");
  fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/D");
 
  //Array branches end
  
  fTree->Branch("trktheta_xz",trktheta_xz,"trktheta_xz[ntracks_reco]/D");
  fTree->Branch("trktheta_yz",trktheta_yz,"trktheta_yz[ntracks_reco]/D");
  fTree->Branch("trketa_xy",trketa_xy,"trketa_xy[ntracks_reco]/D");
  fTree->Branch("trketa_zy",trketa_zy,"trketa_zy[ntracks_reco]/D");
  fTree->Branch("trktheta",trktheta,"trktheta[ntracks_reco]/D");
  fTree->Branch("trkphi",trkphi,"trkphi[ntracks_reco]/D");
  fTree->Branch("trkd2",trkd2,"trkd2[ntracks_reco]/D");

  //Array Branches
  fTree->Branch("trkdedx2",trkdedx2,"trkdedx2[ntracks_reco][3][1000]/D");
  fTree->Branch("trkdqdx",trkdqdx,"trkdqdx[ntracks_reco][3][1000]/D");

  fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][3]/D");

  //Array
  fTree->Branch("trkpitchHit",trkpitchHit,"trkpitchHit[ntracks_reco][3][1000]/D");
  fTree->Branch("trkkinE",trkkinE,"trkkinE[ntracks_reco][3]/D"); 
  fTree->Branch("trkrange",trkrange,"trkrange[ntracks_reco][3]/D");
  
  //Array branches 
  fTree->Branch("trkTPC",trkTPC,"trkTPC[ntracks_reco][3][1000]/D");
  fTree->Branch("trkplaneid",trkplaneid,"trkplaneid[ntracks_reco][3][1000]/D");
  fTree->Branch("trkresrg",trkresrg,"trkresrg[ntracks_reco][3][1000]/D");
  fTree->Branch("trkPosx",trkPosx,"trkPosx[ntracks_reco][3][1000]/D");
  fTree->Branch("trkPosy",trkPosy,"trkPosy[ntracks_reco][3][1000]/D");
  fTree->Branch("trkPosz",trkPosz,"trkPosz[ntracks_reco][3][1000]/D");
  //Array branches end

  fTree->Branch("trklen",trklen,"trklen[ntracks_reco]/D");
  fTree->Branch("trklen_L",trklen_L,"trklen_L[ntracks_reco]/D");
  fTree->Branch("trkdQdxSum",trkdQdxSum,"trkdQdxSum[ntracks_reco]/D");
  fTree->Branch("trkdQdxAverage",trkdQdxAverage,"trkdQdxAverage[ntracks_reco]/D");
  fTree->Branch("trkdEdxSum",trkdEdxSum,"trkdEdxSum[ntracks_reco]/D");
  fTree->Branch("trkdEdxAverage",trkdEdxAverage,"trkdEdxAverage[ntracks_reco]/D");

  fTree->Branch("trkMCTruthT0",trkMCTruthT0,"trkMCTruthT0[ntracks_reco]/D");
 
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("nhits2",&nhits2,"nhits2/I");
  fTree->Branch("nclust",&nclust,"nclust/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits2]/I");
  fTree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits2]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits2]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits2]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits2]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits2]/D");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits2]/D");  
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits2]/I");

}

//void GapWidth::GapWidth::reconfigure(fhicl::ParameterSet const & p)
//{
//  // Implementation of optional member function here.
//}
void GapWidth::GapWidth::ResetVars(){

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  for (int i = 0; i<3; ++i){
    efield[i] = -99999;
  }
  t0 = -99999;
  ntracks_reco = -99999;
  
  for (int i = 0; i < kMaxTrack; ++i){
    trkstartx[i] = -99999;
    trkstarty[i] = -99999;
    trkstartz[i] = -99999;
    trkendx[i] = -99999;
    trkendy[i] = -99999;
    trkendz[i] = -99999;
    
    trkd2[i] = -99999;
    trktheta_xz[i] = -99999;
    trktheta_yz[i] = -99999;
    trketa_xy[i] = -99999;
    trketa_zy[i] = -99999;
    trktheta[i] = -99999;
    trkphi[i] = -99999;
  
    trkdQdxSum[i] = 0;
    trkdEdxSum[i] = 0;
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

    trkenddcosx[i] = -99999;
    trkenddcosy[i] = -99999;
    trkenddcosz[i] = -99999;
    ntrkhits[i] = -99999;
    for (int j = 0; j<kMaxTrackHits; ++j){
      trkx[i][j] = -99999;
      trky[i][j] = -99999;
      trkz[i][j] = -99999;
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
}

void GapWidth::GapWidth::endJob()
{
}

DEFINE_ART_MODULE(GapWidth::GapWidth)
