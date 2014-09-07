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


#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TLorentzVector.h"

const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxHits       = 10000; //maximum number of hits
const int kMaxClust       = 10000; //maximum number of clusters
const int kMaxTrackHits  = 1000;  //maximum number of space points

namespace AnaTree {
  class AnaTree;
}

class AnaTree::AnaTree : public art::EDAnalyzer {
public:
  explicit AnaTree(fhicl::ParameterSet const & p);
  virtual ~AnaTree();

  void analyze(art::Event const & e) override;

  void beginJob() override;
  //void reconfigure(fhicl::ParameterSet const & p) override;

private:

  void ResetVars();
  
  TTree* fTree;
  //run information
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

 double trkstartx_MC[kMaxTrack];
  double trkstarty_MC[kMaxTrack];
  double trkstartz_MC[kMaxTrack];
  double trkendx_MC[kMaxTrack];
  double trkendy_MC[kMaxTrack];
  double trkendz_MC[kMaxTrack];
  double trklen_MC[kMaxTrack];

  double trkmom[kMaxTrack];
  double trkd2[kMaxTrack];
  double trkcolin[kMaxTrack];
  double trklen[kMaxTrack];
  double trkid[kMaxTrack];
  double mcang_x;
  double mcang_y;
  double mcang_z;
  double mcpos_x;
  double mcpos_y;
  double mcpos_z;

  double trkang[kMaxTrack];
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
  int  hit_tpc[kMaxHits];
  int    hit_plane[kMaxHits];
  int    hit_wire[kMaxHits];
  int    hit_channel[kMaxHits];
  double hit_peakT[kMaxHits];
  double hit_charge[kMaxHits];
  double hit_ph[kMaxHits];
  int    hit_trkid[kMaxHits];
  int    nclust;
   int  clust_startw[kMaxClust];
  int    clust_endw[kMaxClust];
  double    clust_startt[kMaxClust];
  double   clust_endt[kMaxClust];
  int   clust_view[kMaxClust];


  std::string fTrigModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fClusterModuleLabel;


};


AnaTree::AnaTree::AnaTree(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fTrigModuleLabel       (pset.get< std::string >("TrigModuleLabel"))
  , fHitsModuleLabel       (pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel       (pset.get< std::string >("TrackModuleLabel"))
  , fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel"))
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

  //  Event.getByLabel(fSimulationProducerLabel, particleHandle);
  //  art::Handle< std::vector<simb::MCParticle> > particleHandle;

  art::ServiceHandle<cheat::BackTracker> bt;
  //  int trkid=1;
  const sim::ParticleList& plist = bt->ParticleList();

  simb::MCParticle *particle=0;
  //  const simb::MCParticle *particle = bt->TrackIDToParticle(trkid);
  for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    particle = ipar->second;
  }

  //  size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
  //  int last = numberTrajectoryPoints - 1;
  //  const TLorentzVector& momentumStart = particle->Momentum(0);
  //  const TLorentzVector& momentumEnd   = particle->Momentum(last);
  //    TLorentzVector tmpVec= momentumEnd-momentumStart;

  //  TLorentzVector tmpVec= momentumStart;
  
    TLorentzVector tmpVec;
    /*
  for ( auto const& particle : (*particleHandle) )
   {
     // we know it is only one MC Particle for now
     size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();
     int last = numberTrajectoryPoints - 1;
     const TLorentzVector& momentumStart = particle.Momentum(0);
     const TLorentzVector& momentumEnd   = particle.Momentum(last);
     tmpVec= momentumEnd-momentumStart;
   }
  */
  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  art::FindMany<recob::Track>       fmtk(hitListHandle, evt, fTrackModuleLabel);

 


  //track information
  ntracks_reco=tracklist.size();

  TLorentzVector vectx(1,0,0,0);
  TLorentzVector vecty(0,1,0,0);
  TLorentzVector vectz(0,0,1,0);
  //mcang_x=(1./TMath::DegToRad())*tmpVec.Angle(vectx.Vect());
  //mcang_x=(1./TMath::DegToRad())*(1.-tmpVec.Phi());
  //  mcang_y=(1./TMath::DegToRad())*tmpVec.Phi();
  //  mcang_z=(1./TMath::DegToRad())*tmpVec.Theta();
  mcang_x=particle->Px();
  mcang_y=particle->Py();
  mcang_z=particle->Pz();


  //  mcpos_x=particle->Vx();
  //  mcpos_y=particle->Vy();
  //  mcpos_z=particle->Vz();
  tmpVec=particle->Position();
  mcpos_x=(1./TMath::DegToRad())*tmpVec.Angle(vectx.Vect());
  mcpos_y=(1./TMath::DegToRad())*tmpVec.Angle(vecty.Vect());
  mcpos_z=(1./TMath::DegToRad())*tmpVec.Angle(vectz.Vect());

  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;

  /////////////////////
  /////////////  Extra Test
  /////////
  ////////////////////

  art::Handle< std::vector<simb::MCParticle> > particleHandle2;
    evt.getByLabel("largeant", particleHandle2);

    // put it in a more easily usable form
    std::vector< art::Ptr<simb::MCParticle> > particles2;
    art::fill_ptr_vector(particles2, particleHandle2);

    art::Handle< std::vector<sim::SimChannel> > simChannelHandle2;
    evt.getByLabel("largeant", simChannelHandle2);    

    //    int fMC_Ntrack = particles2.size();

    float fMC_startXYZT[1000][4];
    float fMC_endXYZT[1000][4];

    int i=0; // track index
    for (auto const& particle: particles2 ) {
      //        int Ndaughters = particle->NumberDaughters();
      //        vector<int> daughters;
      //        for (int i=0; i<Ndaughters; i++) {
      //            daughters.push_back(particle->Daughter(i));
      //        }
      //        fMC_daughters.push_back(daughters);
      size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
      if(!(particle->Process()=="primary" && particle->PdgCode()==13))
	continue;
      double origin[3] = {0.};
      double world[3] = {0.};
      double xyztArray[4];
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


      const TLorentzVector& positionStart = particle->Position(0);
      int last = numberTrajectoryPoints - 1;
      TLorentzVector& positionEnd  =( TLorentzVector&)particle->Position(last);

      bool insideActiveVolume=true;
      for(size_t ii=0;ii<numberTrajectoryPoints;ii++)
	  {
	    const TLorentzVector& tmpPosition=particle->Position(ii);
	    tmpPosition.GetXYZT(xyztArray);
	    for(int p=0;p<3;p++)
	      {
		if((xyztArray[p]>cryoBound_pos[p]) || (xyztArray[p]<cryoBound_neg[p]))
		  {
		    insideActiveVolume=false;
		    break;
		  }
	      }
	    
	    if(!insideActiveVolume) 
	      {
		positionEnd=(const TLorentzVector&)particle->Position(ii-1);
		break;
	      }
	  }
     
	//        const TLorentzVector& momentumStart = particle->Momentum(0);
	//        const TLorentzVector& momentumEnd   = particle->Momentum(last);
      TLorentzVector& momentumStart  =( TLorentzVector&)particle->Momentum(0);
      trkmom[i]=momentumStart.P();
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
	//        momentumStart.GetXYZT(fMC_startMomentum[i]);
	//        momentumEnd.GetXYZT(fMC_endMomentum[i]);
        i++;
    } // particle loop done 
    //////////////////////
    /////////////  End of Extra Test
    //////
    ////////////////////



    // **********************
    // **********************
    //
    //  Tracks:
    //
    //
    //
    // *********************
    // *********************

    
  for(size_t i=0; i<tracklist.size();++i){
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
    //    trkang[i]=TMath::Cos((v2-v1).Angle(tmpVec.Vect()));
    trkang[i]=TMath::Cos((v2-v1).Angle(tmpVec.Vect()));
    //    trklen[i]=v.Mag();
    trkenddcosx[i]    = larEnd[0];
    trkenddcosy[i]    = larEnd[1];
    trkenddcosz[i]    = larEnd[2];
    ntrkhits[i] = fmsp.at(i).size();
    std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
    TVector3 V1(trackStart[0],trackStart[1],trackStart[2]);
    TVector3 V2(trackEnd[0],trackEnd[1],trackEnd[2]);
    TVector3 vOrth=(V2-V1).Orthogonal();

    TVector3 pointVector=V1;
    double distance_squared=0;
    double distance=0;

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

    for (int j = 0; j<3; ++j){
      try{
	if (j==0)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kU);
	else if (j==1)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kV);
	else if (j==2)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kZ);
      }
      catch( cet::exception &e){
	mf::LogWarning("AnaTree")<<"caught exeption "<<e<<"\n setting pitch to 0";
	trkpitch[i][j] = 0;
      }
    }
    TMatrixD rot(3,3);
    int start_point =0;
    simb::MCParticle* particle=0;
    tracklist[i]->GlobalToLocalRotationAtPoint(start_point, rot);
    for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
      particle = ipar->second;
    }
    int pdg = particle->PdgCode();
    if (pdg!=13) continue;
    TVector3 startmom;
    startmom=particle->Momentum(0).Vect();
    TVector3 mcmoml=rot * startmom;
    trkcolin[i]=mcmoml.Z()/mcmoml.Mag();

    
  }
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);
  
  // **********************
  // **********************
  //
  //  Clusters:
  //
  //
  //
  // *********************
  // *********************
  
  nclust=clusterlist.size();
  for (size_t iclu = 0; iclu<clusterlist.size(); ++iclu){
    double w0 = clusterlist[iclu]->StartPos()[0];
    double w1 = clusterlist[iclu]->EndPos()[0];
    double t0 = clusterlist[iclu]->StartPos()[1];
    double t1 = clusterlist[iclu]->EndPos()[1];
    int tmpview = clusterlist[iclu]->View();
    
    clust_startw[iclu]=w0;
    clust_endw[iclu]=w1;
    clust_startt[iclu]=t0;
    clust_endt[iclu]=t1;
    clust_view[iclu]=tmpview;

    }

  
  // **********************
  // **********************
  //
  //  Hits:
  //
  //
  //
  // *********************
  // *********************
  
  
  nhits = hitlist.size();
  for (size_t i = 0; i<hitlist.size(); ++i){
    unsigned int channel = hitlist[i]->Channel();
    geo::WireID wireid = hitlist[i]->WireID();
    hit_tpc[i]     =wireid.TPC;
    hit_plane[i]   = wireid.Plane;
    hit_wire[i]    = wireid.Wire;
    hit_channel[i] = channel;
    hit_peakT[i]   = hitlist[i]->PeakTime();
    hit_charge[i]  = hitlist[i]->Charge();
    hit_ph[i]      = hitlist[i]->Charge(true);
    if (fmtk.at(i).size()!=0){
      hit_trkid[i] = fmtk.at(i)[0]->ID();
    }
  }
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
  fTree->Branch("trkstartx",trkstartx,"trkstartx[ntracks_reco]/D");
  fTree->Branch("trkstarty",trkstarty,"trkstarty[ntracks_reco]/D");
  fTree->Branch("trkstartz",trkstartz,"trkstartz[ntracks_reco]/D");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/D");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/D");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/D");
fTree->Branch("trkstartx_MC",trkstartx_MC,"trkstartx_MC[ntracks_reco]/D");
  fTree->Branch("trkstarty_MC",trkstarty_MC,"trkstarty_MC[ntracks_reco]/D");
  fTree->Branch("trkstartz_MC",trkstartz_MC,"trkstartz_MC[ntracks_reco]/D");
  fTree->Branch("trkendx_MC",trkendx_MC,"trkendx_MC[ntracks_reco]/D");
  fTree->Branch("trkendy_MC",trkendy_MC,"trkendy_MC[ntracks_reco]/D");
  fTree->Branch("trkendz_MC",trkendz_MC,"trkendz_MC[ntracks_reco]/D");
  fTree->Branch("trklen_MC",trklen_MC,"trklen_MC[ntracks_reco]/D");
  fTree->Branch("trkmom",trkmom,"trkmom[ntracks_reco]/D");
  fTree->Branch("trkd2",trkd2,"trkd2[ntracks_reco]/D");
  fTree->Branch("trkcolin",trkcolin,"trkcolin[ntracks_reco]/D");

  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/D");
  fTree->Branch("trklen",trklen,"trklen[ntracks_reco]/D");
  fTree->Branch("trkid",trkid,"trkid[ntracks_reco]/D");
  fTree->Branch("mcang_x",&mcang_x,"mcang_x/D");
  fTree->Branch("mcang_y",&mcang_y,"mcang_y/D");
  fTree->Branch("mcang_z",&mcang_z,"mcang_z/D");

  fTree->Branch("mcpos_x",&mcpos_x,"mcpos_x/D");
  fTree->Branch("mcpos_y",&mcpos_y,"mcpos_y/D");
  fTree->Branch("mcpos_z",&mcpos_z,"mcpos_z/D");

  fTree->Branch("trkang",trkang,"trkang[ntracks_reco]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/D");
  fTree->Branch("ntrkhits",ntrkhits,"ntrkhits[ntracks_reco]/I");
  fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/D");
  fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/D");
  fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/D");
  fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][3]/D");
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I");
  fTree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/D");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/D");
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I");

 fTree->Branch("nclust",&nclust,"nclust/I");
  fTree->Branch("clust_startw",clust_startw,"clust_startw[nclust]/I");
  fTree->Branch("clust_endw",clust_endw,"clust_endw[nclust]/I");
  fTree->Branch("clust_view",clust_view,"clust_view[nclust]/I");
  fTree->Branch("clust_startt",clust_startt,"clust_startt[nclust]/D");
  fTree->Branch("clust_endt",clust_endt,"hclust_endt[nclust]/D");
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
  mcang_x = -99999;
  mcang_y = -99999;
  mcang_z = -99999;

  mcpos_x = -99999;
  mcpos_y = -99999;
  mcpos_z = -99999;


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
    trkmom[i] = -99999;
    trkd2[i] = -99999;
    trkcolin[i] = -99999;

    trkstartdcosx[i] = -99999;
    trkstartdcosy[i] = -99999;
    trkstartdcosz[i] = -99999;
    trklen[i] = -99999;
    trkid[i] = -99999;
    trkang[i] = -99999;
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
   for (int i = 0; i<kMaxClust; ++i){
    clust_startw[i] = -99999;
    clust_endw[i] = -99999;
    clust_startt[i] = -99999;
    clust_endt[i] = -99999;
    clust_view[i] = -99999;

  }
}

DEFINE_ART_MODULE(AnaTree::AnaTree)
