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
const int kMaxClust       = 10000; //maximum number of clusters
const int kMaxTrackHits  = 1000;  //maximum number of space points

namespace AnaTree {
  class AnaTree;
}
namespace {

  // Local functions.

  // Calculate distance to boundary.

  //----------------------------------------------------------------------------

  double bdist(const TVector3& pos, unsigned int tpc = 0, unsigned int /*cstat*/ = 0)

  {

    // Get geometry.

    double d3,d4;

    art::ServiceHandle<geo::Geometry> geom;

      

    if(tpc==2 || tpc==3 || tpc==4 || tpc==5)

      {

        d3 = pos.Y() - 1.25;     // 1.25cm APA 2/3 distance to horizontal.

        //    double d3 = pos.Y() + 85.0;     // Distance to bottom.

        d4 = 113.0 - pos.Y();     // Distance to top.

        //    double d4 = 113.0 - pos.Y();     // Distance to top.

      }

    else  //tpc 0 1  6  7

      {

        d3 = pos.Y() + geom->DetHalfHeight(tpc)-15.0;     // Distance to bottom.

        //    double d3 = pos.Y() + 85.0;     // Distance to bottom.

        d4 = geom->DetHalfHeight(tpc)+15.0 - pos.Y();     // Distance to top.

        //    double d4 = 113.0 - pos.Y();     // Distance to top.

      }

    

    //    mf::LogVerbatim("output") <<"d3" << d3;

    //    mf::LogVerbatim("output") <<"d4" << d4;

    double d1 = abs(pos.X());          // Distance to right side (wires).

    double d2=2.*geom->DetHalfWidth(tpc)- abs(pos.X());

    //    mf::LogVerbatim("output") <<"d1" << d1;

    //    mf::LogVerbatim("output") <<"d2" << d2;

    //    double d2 = 226.539 - pos.X();   // Distance to left side (cathode).

    double d5 = 0.;

    double d6 = 0.;

    

    if(tpc==0 || tpc==1)

      {

        d5 = pos.Z()+1.0;                             // Distance to front.

        d6 = geom->DetLength(tpc) -1.0- pos.Z();         // Distance to back.

      }

    else if (tpc==2||tpc==3 || tpc==4 || tpc==5)

      {

        d5 = pos.Z()-51.0;                             // Distance to front.     

        d6 = geom->DetLength(tpc) +51.0- pos.Z();         // Distance to back.   

      }

    else if (tpc==6 || tpc==7)

      {

        d5 = pos.Z()-103.0;                             // Distance to front.     

        d6 = geom->DetLength(tpc) +103.0- pos.Z();         // Distance to back.   

        

      }

    if(d6<0){

      //      mf::LogVerbatim("output")<< "z"  <<pos.Z();

      //      mf::LogVerbatim("output")<< "Tpc" <<tpc;

      //      mf::LogVerbatim("output")<< "DetLength" <<geom->DetLength(tpc);

      

    }

    //    mf::LogVerbatim("output") <<"d5" << d5;

    //    mf::LogVerbatim("output") <<"d6" << d6;

    double result = std::min(std::min(std::min(std::min(std::min(d1, d2), d3), d4), d5), d6);

    //    mf::LogVerbatim("output")<< "bdist" << result;

    //    mf::LogVerbatim("output")<< "Height" << geom->DetHalfHeight(tpc);

    //    mf::LogVerbatim("output")<< "Width" << geom->DetHalfWidth(tpc);

    if(result<0) result=0;

    return result;
  }
  // Length of reconstructed track.
  //----------------------------------------------------------------------------
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

  void effcalc(const TH1* hnum, const TH1* hden, TH1* heff)

  {

    int nbins = hnum->GetNbinsX();

    if (nbins != hden->GetNbinsX())

      throw cet::exception("TrackAnaCT") << "effcalc[" __FILE__ "]: incompatible histograms (I)\n";

    if (nbins != heff->GetNbinsX())

      throw cet::exception("TrackAnaCT") << "effcalc[" __FILE__ "]: incompatible histograms (II)\n";

    // Loop over bins, including underflow and overflow.

    for(int ibin = 0; ibin <= nbins+1; ++ibin) {

      double num = hnum->GetBinContent(ibin);

      double den = hden->GetBinContent(ibin);

      if(den == 0.) {

        heff->SetBinContent(ibin, 0.);

        heff->SetBinError(ibin, 0.);

      }

      else {

        double eff = num / den;

        if(eff < 0.)

          eff = 0.;

        if(eff > 1.)

          eff = 1.;

        double err = std::sqrt(eff * (1.-eff) / den);

        heff->SetBinContent(ibin, eff);

        heff->SetBinError(ibin, err);

      }

    }

    heff->SetMinimum(0.);

    heff->SetMaximum(1.05);

    heff->SetMarkerStyle(20);

  }

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
  double trklen_cut_MC[kMaxTrack];

  double trkmom_MC[kMaxTrack];
  double trkmom_XMC[kMaxTrack];
  double trkmom_YMC[kMaxTrack];
  double trkmom_ZMC[kMaxTrack];
  double trkstartdoc_XMC[kMaxTrack];
  double trkstartdoc_YMC[kMaxTrack];
  double trkstartdoc_ZMC[kMaxTrack];
  double trkpdg_MC[kMaxTrack];
  double trkd2[kMaxTrack];
  double trkcolin[kMaxTrack];
  double trklen[kMaxTrack];
  double trklen_L[kMaxTrack];
  double trkid[kMaxTrack];
  double trktheta_xz_MC[kMaxTrack];
  double trktheta_yz_MC[kMaxTrack];
  double trketa_xy_MC[kMaxTrack];
  double trketa_zy_MC[kMaxTrack];
  double trktheta_MC[kMaxTrack];
  double trkphi_MC[kMaxTrack];

  double trktheta_xz[kMaxTrack];
  double trketa_xy[kMaxTrack];
  double trktheta_yz[kMaxTrack];
  double trketa_zy[kMaxTrack];
  double trktheta[kMaxTrack];
  double trkphi[kMaxTrack];

  double trkdedx[kMaxTrack];
  double trkdedx2[kMaxTrack][3][1000];
  double trkdqdx[kMaxTrack][3][1000];
  double trkpitchHit[kMaxTrack][3][1000];
  double trkkinE[kMaxTrack][3];
  double trkplaneid[kMaxTrack][3][1000];
  double trkresrg[kMaxTrack][3][1000];
  double trkdedx_MC[kMaxTrack];
  double trkdq_MC[kMaxTrack];
  double mcang_x;
  double mcang_y;
  double mcang_z;
  double mcpos_x;
  double mcpos_y;
  double mcpos_z;

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

  struct RecoHists
    {
      // Constructors.
      RecoHists();
      //      ~RecoHists();
      RecoHists(const std::string& subdir);
      // Pure reco track histograms.
      TH1F* fHstartx;      // Starting x position.
      TH1F* fHstarty;      // Starting y position.
      TH1F* fHstartz;      // Starting z position.
      TH1F* fHstartd;      // Starting distance to boundary.
      TH1F* fHendx;        // Ending x position.
      TH1F* fHendy;        // Ending y position.
      TH1F* fHendz;        // Ending z position.
      TH1F* fHendd;        // Ending distance to boundary.
      TH1F* fHtheta;       // Theta.
      TH1F* fHphi;         // Phi.
      TH1F* fHtheta_xz;    // Theta_xz.
      TH1F* fHtheta_yz;    // Theta_yz.
      TH1F* fHmom;         // Momentum.
      TH1F* fHlen;         // Length.
      // Histograms for the consituent Hits
      TH1F* fHHitChg;       // hit charge
      TH1F* fHHitWidth;     // hit width
      TH1F* fHHitPdg;       // Pdg primarily responsible.
      TH1F* fHHitTrkId;     // TrkId
      TH1F* fModeFrac;       // mode fraction
      TH1F* fNTrkIdTrks;    // # of stitched tracks in which unique TrkId appears
      TH2F* fNTrkIdTrks2;   
      TH2F* fNTrkIdTrks3;   
    };

    // Struct for mc particles and mc-matched tracks.
    struct MCHists
    {
      // Constructors.
      MCHists();
      MCHists(const std::string& subdir);
      // Reco-MC matching.
      TH2F* fHduvcosth;    // 2D mc vs. data matching, duv vs. cos(theta).
      TH1F* fHcosth;       // 1D direction matching, cos(theta).
      TH1F* fHmcu;         // 1D endpoint truth u.
      TH1F* fHmcv;         // 1D endpoint truth v.
      TH1F* fHmcw;         // 1D endpoint truth w.
      TH1F* fHupull;       // 1D endpoint u pull.
      TH1F* fHvpull;       // 1D endpoint v pull.
      TH1F* fHmcdudw;      // Truth du/dw.
      TH1F* fHmcdvdw;      // Truth dv/dw.
      TH1F* fHdudwpull;    // du/dw pull.
      TH1F* fHdvdwpull;    // dv/dw pull.
      TH1F* fHHitEff;      // Hit efficiency.
      TH1F* fHHitPurity;   // Hit purity.

      // Histograms for matched tracks.
      TH1F* fHstartdx;     // Start dx.
      TH1F* fHstartdy;     // Start dy.
      TH1F* fHstartdz;     // Start dz.
      TH1F* fHenddx;       // End dx.
      TH1F* fHenddy;       // End dy.
      TH1F* fHenddz;       // End dz.
      TH2F* fHlvsl;        // MC vs. reco length.
      TH1F* fHdl;          // Delta(length).
      TH2F* fHpvsp;        // MC vs. reco momentum.
      TH2F* fHpvspc;       // MC vs. reco momentum (contained tracks).
      TH1F* fHdp;          // Momentum difference.
      TH1F* fHdpc;         // Momentum difference (contained tracks).
      TH1F* fHppull;       // Momentum pull.
      TH1F* fHppullc;      // Momentum pull (contained tracks).
      // Pure MC particle histograms (efficiency denominator).
      TH1F* fHmcstartx;    // Starting x position.
      TH1F* fHmcstarty;    // Starting y position.
      TH1F* fHmcstartz;    // Starting z position.
      TH1F* fHmcendx;      // Ending x position.
      TH1F* fHmcendy;      // Ending y position.
      TH1F* fHmcendz;      // Ending z position.
      TH1F* fHmctheta;     // Theta.
      TH1F* fHmcphi;       // Phi.
      TH1F* fHmctheta_xz;  // Theta_xz.
      TH1F* fHmctheta_yz;  // Theta_yz.
      TH1F* fHmcmom;       // Momentum.
      TH1F* fHmclen;       // Length.
      // Histograms for well-reconstructed matched tracks (efficiency numerator).
      TH1F* fHgstartx;     // Starting x position.
      TH1F* fHgstarty;     // Starting y position.
      TH1F* fHgstartz;     // Starting z position.
      TH1F* fHgendx;       // Ending x position.
      TH1F* fHgendy;       // Ending y position.
      TH1F* fHgendz;       // Ending z position.
      TH1F* fHgtheta;      // Theta.
      TH1F* fHgphi;        // Phi.
      TH1F* fHgtheta_xz;   // Theta_xz.
      TH1F* fHgtheta_yz;   // Theta_yz.
      TH1F* fHgmom;        // Momentum.
      TH1F* fHglen;        // Length.
      // Efficiency histograms.
      TH1F* fHestartx;     // Starting x position.
      TH1F* fHestarty;     // Starting y position.
      TH1F* fHestartz;     // Starting z position.
      TH1F* fHeendx;       // Ending x position.
      TH1F* fHeendy;       // Ending y position.
      TH1F* fHeendz;       // Ending z position.
      TH1F* fHetheta;      // Theta.
      TH1F* fHephi;        // Phi.
      TH1F* fHetheta_xz;   // Theta_xz.
      TH1F* fHetheta_yz;   // Theta_yz.
      TH1F* fHemom;        // Momentum.
      TH1F* fHelen;        // Length.
    };
    std::map<int, MCHists> fMCHistMap;       // Indexed by pdg id.
    std::map<int, RecoHists> fRecoHistMap;   // Indexed by pdg id.

};


AnaTree::AnaTree::AnaTree(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fTrigModuleLabel       (pset.get< std::string >("TrigModuleLabel"))
  , fHitsModuleLabel       (pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel       (pset.get< std::string >("TrackModuleLabel"))
  , fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel"))
  , fTrkSpptAssocModuleLabel    (pset.get< std::string >("TrkSpptAssocModuleLabel"))
  , fHitSpptAssocModuleLabel    (pset.get< std::string >("HitSpptAssocModuleLabel"))
  ,  fSimulationProducerLabel ( pset.get< std::string >("SimulationLabel"))
  ,  fCalorimetryModuleLabel ( pset.get< std::string >("CalorimetryModuleLabel"))
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

  art::ServiceHandle<cheat::BackTracker> bt;
  //  int trkid=1;
  const sim::ParticleList& plist = bt->ParticleList();

  simb::MCParticle *particle=0;
  //  const simb::MCParticle *particle = bt->TrackIDToParticle(trkid);
  for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    particle = ipar->second;
    
    if(!(particle->Process()=="primary" && abs(particle->PdgCode())== abs(fPdg))) continue;

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
    evt.getByLabel(fSimulationProducerLabel, particleHandle2);

    // put it in a more easily usable form
    std::vector< art::Ptr<simb::MCParticle> > particles2;
    art::fill_ptr_vector(particles2, particleHandle2);

    //    int fMC_Ntrack = particles2.size();
    float fMC_startXYZT[1000][4];
    float fMC_endXYZT[1000][4];

    double trkde=0;

    int i=0; // track index
    for (auto const& particle: particles2 ) {
      //        int Ndaughters = particle->NumberDaughters();
      //        vector<int> daughters;
      //        for (int i=0; i<Ndaughters; i++) {
      //            daughters.push_back(particle->Daughter(i));
      //        }
      //        fMC_daughters.push_back(daughters);
      size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
      trkpdg_MC[i]=particle->PdgCode();
      if(!(particle->Process()=="primary" && abs(particle->PdgCode())==abs(fPdg)))
	continue;
      int trackID=particle->TrackId();
      for ( auto const& channel : (*simChannelHandle) )
	{
	  auto channelNumber = channel.Channel();
	  if ( fGeometry->SignalType( channelNumber ) == geo::kCollection )
	    {
	      auto const& timeSlices = channel.TDCIDEMap();
	      for ( auto const& timeSlice : timeSlices )
		{
		  auto const& energyDeposits = timeSlice.second;
		  for ( auto const& energyDeposit : energyDeposits )
		    {
		      if ( energyDeposit.trackID == trackID )
			{
			  
			  trkde+=energyDeposit.numElectrons*1./( 1e6/23.6);//MeV
			}//TrackID
		    }//energyDeposit
		}//timeSlice
	    }//CollectionPlane
	}//simChannel

      trkdq_MC[i]=trkde*1e6/23.6;//back to q
      
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

      int zz =0;
      int zz2 =0;
      bool insideActiveVolume=false;
      for(size_t ii=0;ii<numberTrajectoryPoints;ii++)
	  {
	    const TLorentzVector& tmpPosition=particle->Position(ii);
	    //tmpPosition.GetXYZT(xyztArray);   
	    if((tmpPosition[0]<cryoBound_pos[0]) && (tmpPosition[0]>cryoBound_neg[0])) {
	      if((tmpPosition[1]<cryoBound_pos[1]) && (tmpPosition[1]>cryoBound_neg[1])) {
		if((tmpPosition[2]<cryoBound_pos[2]) && (tmpPosition[2]>cryoBound_neg[2])) {
		  if (!insideActiveVolume) {
		    zz = ii;
		    //std::cout << "Now particle is in cryostat " << std::endl;
		    insideActiveVolume=true;
		  }		  
		  //std::cout << "Temp Pos " << tmpPosition[0] << ", " <<  tmpPosition[1] << ", " <<  tmpPosition[2] << std::endl;
		  tmpPosition.GetXYZT(xyztArray);
		  zz2 = ii;
		}
	      }
	    }
	    
	    if ( (insideActiveVolume) && (zz2 != (int)ii) ) break;
	  }

      const TLorentzVector& positionStart = particle->Position(zz);
      TLorentzVector& positionEnd  =( TLorentzVector&)particle->Position(zz2);     
	//        const TLorentzVector& momentumStart = particle->Momentum(0);
	//        const TLorentzVector& momentumEnd   = particle->Momentum(last);
      TLorentzVector& momentumStart  =( TLorentzVector&)particle->Momentum(0);
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
      double mctime = particle->T();                                 // nsec
      double mcdx = mctime * 1.e-3 * larprop->DriftVelocity();   // cm
      // Calculate the points where this mc particle enters and leaves the
      // fiducial volume, and the length in the fiducial volume.
      TVector3 mcstart;
      TVector3 mcend;
      TVector3 mcstartmom;
      TVector3 mcendmom;
      double plen = length(*particle, mcdx, mcstart, mcend, mcstartmom, mcendmom);
      trklen_cut_MC[i]=plen;
      trkdedx_MC[i]=trkde/trklen_cut_MC[i];
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
    trkdedx[i]=trkde/trklen[i];
    //    trkang[i]=TMath::Cos((v2-v1).Angle(tmpVec.Vect()));
    trkang[i]=TMath::Cos((v2-v1).Angle(tmpVec.Vect()));
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
      
      art::FindMany<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
      if (fmcal.isValid()){
        std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
        //std::cout<<"calo size "<<calos.size()<<std::endl;
        for (size_t jj = 0; jj<calos.size(); ++jj){
	  //          trkke[it1][i][j]    = calos[j]->KineticEnergy();
	  //	  trkrange[i][jj] = calos[jj]->Range();
	  //          trkpitchc[it1][i][j]= calos[j] -> TrkPitchC();
	  //   ntrkhitscal[i][jj] = calos[jj] -> dEdx().size();
	  trkkinE[i][jj]  = (calos[jj] -> KineticEnergy());
	  //std::cout << trkkinE[i][jj] << std::endl;
	  int tt= calos[jj] -> dEdx().size();
	  for(int k = 0; k < tt; ++k) {
            trkdedx2[i][jj][k]  = (calos[jj] -> dEdx())[k];
	    trkdqdx[i][jj][k]   = (calos[jj] -> dQdx())[k];
	    trkpitchHit[i][jj][k]  = (calos[jj] -> TrkPitchVec())[k];   	    
            trkresrg[i][jj][k]  = (calos[jj] -> ResidualRange())[k];
	    trkplaneid[i][jj][k]=(calos[jj]->PlaneID()).Plane;	    
	    //            trkxyz[it1][i][j][k][0] = (calos[jj] -> XYZ())[k].X();
	    //            trkxyz[it1][i][j][k][1] = (calos[jj] -> XYZ())[k].Y();
	    //            trkxyz[it1][i][j][k][2] = (calos[jj] -> XYZ())[k].Z();
          }
        }
      }


      // *********************
      //  End Calorimetric stuff
      //  
      // *********************
      



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
  double mctime = particle->T();                                 // nsec
  double mcdx = mctime * 1.e-3 * larprop->DriftVelocity();   // cm
  // Calculate the points where this mc particle enters and leaves the
  // fiducial volume, and the length in the fiducial volume.
  TVector3 mcstart;
  TVector3 mcend;
  TVector3 mcstartmom;
  TVector3 mcendmom;
  double plen = length(*particle, mcdx, mcstart, mcend, mcstartmom, mcendmom);
  // Get the displacement of this mc particle in the global coordinate system.
  //  TVector3 mcpos = mcstart - pos;
  TVector3 mcpos = pos -mcstart ;
  // Rotate the momentum and position to the
  // track-local coordinate system.
  TVector3 mcmoml = rot * mcstartmom;
  TVector3 mcposl = rot * mcpos;
  trktheta_xz_MC[i] = std::atan2(mcstartmom.X(), mcstartmom.Z());
  trktheta_yz_MC[i] = std::atan2(mcstartmom.Y(), mcstartmom.Z());
  trketa_xy_MC[i] = std::atan2(mcstartmom.X(), mcstartmom.Y());
  trketa_zy_MC[i] = std::atan2(mcstartmom.Z(), mcstartmom.Y());



  trktheta_MC[i]=mcstartmom.Theta();
  trkphi_MC[i]=mcstartmom.Phi();

  trkcolinearity[i] = mcmoml.Z() / mcmoml.Mag();
  double u = mcposl.X();
  double v = mcposl.Y();
  double w = mcposl.Z();
  trkwmatchdisp[i]=w;
  /*  std::cout << "++++++" << std::endl;
  std::cout << "w " << w << std::endl;
 std::cout << "trkcolinearity  " << trkcolinearity[i] << std::endl;
 std::cout << "plen  " << plen << std::endl;
 std::cout << "mcstartmom mag  " << mcstartmom.Mag() << std::endl;
 
 std::cout << "++++++" << std::endl;*/
  double pu = mcmoml.X();
  double pv = mcmoml.Y();
  double pw = mcmoml.Z();
  double dudw = pu / pw;
  double dvdw = pv / pw;
  /*  std::cout << "pu  "<<pu << "pv  " << pv << std::endl;
  std::cout << "pw  "<<pw  << std::endl;

  std::cout << "u  "<<u << "v  " << v << std::endl;*/


  double u0 = u - w * dudw;
  double v0 = v - w * dvdw;
  trkmatchdisp[i]=abs( std::sqrt(u0*u0 + v0*v0));
  art::Ptr<recob::Track> ptrack(trackh, i);
  const recob::Track& track = *ptrack;
  trklenratio[i] = length(track)/plen;
  trklen_L[i]=length(track);
  //
  //**********************
  //
  // End Cut specific quantities
  //  
  //***********************



 
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

  simb::MCParticle* particle=0;
  for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    particle = ipar->second;
  }
  int pdg = particle->PdgCode();
  if (abs(pdg)!=abs(fPdg)) continue;
  TVector3 startmom;
  startmom=particle->Momentum(0).Vect();
  TVector3 mcmomltemp=rot * startmom;
  trkcolin[i]=mcmomltemp.Z()/mcmomltemp.Mag();
  }


  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  nhits = hitlist.size();
  nclust=clusterlist.size();
  for (size_t i = 0; i<hitlist.size(); ++i){
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

  // **********************
  // **********************
  //   Fill Tree:
  // **********************
  // **********************
  
  fTree->Fill();
  // *********************
  // *********************



  // **********************
  // **********************
  //
  //  Histograms:
  //
  //
  //
  // *********************
  // *********************
 
  std::unique_ptr<mf::LogInfo> pdump;
  if(fDump > 0) {
    --fDump;
    pdump = std::unique_ptr<mf::LogInfo>(new mf::LogInfo("TrackAnaCT"));
  }
  // Histograms.

  art::ServiceHandle<cheat::BackTracker> bt2;
  const   sim::ParticleList& pplist=bt2->ParticleList();
  std::vector<const simb::MCParticle*> plist2;
  plist2.reserve(pplist.size());
  if(particle->E() >= 0.001*particle->Mass() + fMinMCKE) {
    // Calculate the x offset due to nonzero mc particle time.
    double mctime = particle->T();                                // nsec
    double mcdx = mctime * 1.e-3 * larprop->DriftVelocity();  // cm
    // Calculate the length of this mc particle inside the fiducial volume.
    TVector3 mcstart;
    TVector3 mcend;
    TVector3 mcstartmom;
    TVector3 mcendmom;
    double plen = length(*particle, mcdx, mcstart, mcend, mcstartmom, mcendmom);
    // Apply minimum fiducial length cut.  Always reject particles that have
    // zero length in the tpc regardless of the configured cut.
    if(plen > 0. && plen > fMinMCLen) {
      // This is a good mc particle (capable of making a track).
      plist2.push_back(particle);
      // Dump MC particle information here.
      if(pdump) {
	//do nothing
      }
      // Fill histograms.
      int pdg=fPdg;
      if(fMCHistMap.count(pdg) == 0) {
	std::ostringstream ostr;
	ostr << "MC" << (fIgnoreSign ? "All" : (pdg > 0 ? "Pos" : "Neg")) << std::abs(pdg);
	fMCHistMap[pdg] = MCHists(ostr.str());
      }

      const MCHists& mchists = fMCHistMap[pdg];
      double mctheta_xz = std::atan2(mcstartmom.X(), mcstartmom.Z());
      double mctheta_yz = std::atan2(mcstartmom.Y(), mcstartmom.Z());

      mchists.fHmcstartx->Fill(mcstart.X());
      mchists.fHmcstarty->Fill(mcstart.Y());
      mchists.fHmcstartz->Fill(mcstart.Z());
      mchists.fHmcendx->Fill(mcend.X());
      mchists.fHmcendy->Fill(mcend.Y());
      mchists.fHmcendz->Fill(mcend.Z());
      mchists.fHmctheta->Fill(mcstartmom.Theta());
      mchists.fHmcphi->Fill(mcstartmom.Phi());
      mchists.fHmctheta_xz->Fill(mctheta_xz);
      mchists.fHmctheta_yz->Fill(mctheta_yz);
      mchists.fHmcmom->Fill(mcstartmom.Mag());
      mchists.fHmclen->Fill(plen);
    }
  } // ...+minMCKE

  //    art::Handle< std::vector<recob::Track> > trackh;
  //    evt.getByLabel(fTrackModuleLabel, trackh);
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
	int hit_tpc=-1;
	if(hits.size()!=0)
	  {
	    geo::WireID tmpWireid=hits.at(0)->WireID();
	    hit_tpc=tmpWireid.TPC;
	    }
	  else hit_tpc=1; */
	int hit_tpc=1;


         ///
        ///
        //
        //
        //
        //
        //
        //
        /*        art::Handle< std::vector<recob::Hit> > hitListHandle;
        std::vector<art::Ptr<recob::Hit> > hitlist;
        if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
        art::fill_ptr_vector(hitlist, hitListHandle);*/
        // Calculate the x offset due to nonzero reconstructed time.
        //double recotime = track.Time() * detprop->SamplingRate();       // nsec
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
          double dpos = bdist(pos,hit_tpc);
          double dend = bdist(end,hit_tpc);
          double tlen = length(track);
          double theta_xz = std::atan2(dir.X(), dir.Z());
          double theta_yz = std::atan2(dir.Y(), dir.Z());
	  trktheta_xz[i] = std::atan2(dir.X(), dir.Z());
	  trktheta_yz[i] = std::atan2(dir.Y(), dir.Z());
 	  trketa_xy[i] = std::atan2(dir.X(), dir.Y());
 	  trketa_zy[i] = std::atan2(dir.Z(), dir.Y());
 	  trktheta[i]=dir.Theta();
 	  trkphi[i]=dir.Phi();
          if(fRecoHistMap.count(0) == 0)
            fRecoHistMap[0] = RecoHists("Reco");
          const RecoHists& rhists = fRecoHistMap[0];
          rhists.fHstartx->Fill(pos.X());
          rhists.fHstarty->Fill(pos.Y());
          rhists.fHstartz->Fill(pos.Z());
          rhists.fHstartd->Fill(dpos);
          rhists.fHendx->Fill(end.X());
          rhists.fHendy->Fill(end.Y());
          rhists.fHendz->Fill(end.Z());
          rhists.fHendd->Fill(dend);
          rhists.fHtheta->Fill(dir.Theta());
          rhists.fHphi->Fill(dir.Phi());
          rhists.fHtheta_xz->Fill(theta_xz);
          rhists.fHtheta_yz->Fill(theta_yz);
          double mom = 0.;
          if(track.NumberFitMomentum() > 0)
            mom = track.VertexMomentum();
          rhists.fHmom->Fill(mom);
          rhists.fHlen->Fill(tlen);
          // Id of matching mc particle.
	  int mcid = -1;
          // Loop over direction.  
          for(int swap=0; swap<2; ++swap) {
            // Analyze reversed tracks only if start momentum = end momentum.
            if(swap != 0 && track.NumberFitMomentum() > 0 &&
               std::abs(track.VertexMomentum() - track.EndMomentum()) > 1.e-3)
              continue;
            // Calculate the global-to-local rotation matrix.
            TMatrixD rot(3,3);
            int start_point = (swap == 0 ? 0 : ntraj-1);
            track.GlobalToLocalRotationAtPoint(start_point, rot);
            // Update track data for reversed case.
            if(swap != 0) {
              rot(1, 0) = -rot(1, 0);
              rot(2, 0) = -rot(2, 0);
              rot(1, 1) = -rot(1, 1);
              rot(2, 1) = -rot(2, 1);
              rot(1, 2) = -rot(1, 2);
              rot(2, 2) = -rot(2, 2);
              pos = track.End();
              dir = -track.EndDirection();
              end = track.Vertex();
              pos[0] += trackdx;
              end[0] += trackdx;
              dpos = bdist(pos,hit_tpc);
              dend = bdist(end,hit_tpc);
              theta_xz = std::atan2(dir.X(), dir.Z());
              theta_yz = std::atan2(dir.Y(), dir.Z());
              if(track.NumberFitMomentum() > 0)
                mom = track.EndMomentum();
            }
	    
            // Get covariance matrix.
            //            const TMatrixD& cov = (swap == 0 ? track.VertexCovariance() : track.EndCovariance());
            // Loop over track-like mc particles.
            for(auto ipart = plist2.begin(); ipart != plist2.end(); ++ipart) {
              const simb::MCParticle* part = *ipart;
              if (!part)
                throw cet::exception("SeedAna") << "no particle! [II]\n";
              int pdg = part->PdgCode();
              if(fIgnoreSign) pdg = std::abs(pdg);
              auto iMCHistMap = fMCHistMap.find(pdg);
              if (iMCHistMap == fMCHistMap.end())
                throw cet::exception("SeedAna") << "no particle with ID=" << pdg << "\n";
              const MCHists& mchists = iMCHistMap->second;
              // Calculate the x offset due to nonzero mc particle time.
              double mctime = part->T();                                 // nsec
              double mcdx = mctime * 1.e-3 * larprop->DriftVelocity();   // cm
              // Calculate the points where this mc particle enters and leaves the
              // fiducial volume, and the length in the fiducial volume.
              TVector3 mcstart;
              TVector3 mcend;
              TVector3 mcstartmom;
              TVector3 mcendmom;
              double plen = length(*part, mcdx, mcstart, mcend, mcstartmom, mcendmom);
              // Get the displacement of this mc particle in the global coordinate system.
              TVector3 mcpos = mcstart - pos;
              // Rotate the momentum and position to the
              // track-local coordinate system.
              TVector3 mcmoml = rot * mcstartmom;
              TVector3 mcposl = rot * mcpos;
              double colinearity = mcmoml.Z() / mcmoml.Mag();
              double u = mcposl.X();
              double v = mcposl.Y();
              double w = mcposl.Z();
              double pu = mcmoml.X();
              double pv = mcmoml.Y();
              double pw = mcmoml.Z();
              double dudw = pu / pw;
              double dvdw = pv / pw;
              double u0 = u - w * dudw;
              double v0 = v - w * dvdw;
              double uv0 = std::sqrt(u0*u0 + v0*v0);
              mchists.fHduvcosth->Fill(colinearity, uv0);
              if(std::abs(uv0) < fMatchDisp) {
                // Fill slope matching histograms.
                mchists.fHmcdudw->Fill(dudw);
                mchists.fHmcdvdw->Fill(dvdw);
                //                mchists.fHdudwpull->Fill(dudw / std::sqrt(cov(2,2)));
                //                mchists.fHdvdwpull->Fill(dvdw / std::sqrt(cov(3,3)));
              }
              mchists.fHcosth->Fill(colinearity);
              if(colinearity > fMatchColinearity) {
                // Fill displacement matching histograms.
                mchists.fHmcu->Fill(u0);
                mchists.fHmcv->Fill(v0);
                mchists.fHmcw->Fill(w);
                //                mchists.fHupull->Fill(u0 / std::sqrt(cov(0,0)));
                //                mchists.fHvpull->Fill(v0 / std::sqrt(cov(1,1)));

                if(std::abs(uv0) < fMatchDisp) {
                  // Fill matching histograms.
                  double mctheta_xz = std::atan2(mcstartmom.X(), mcstartmom.Z());
                  double mctheta_yz = std::atan2(mcstartmom.Y(), mcstartmom.Z());
                  mchists.fHstartdx->Fill(pos.X() - mcstart.X());
                  mchists.fHstartdy->Fill(pos.Y() - mcstart.Y());
                  mchists.fHstartdz->Fill(pos.Z() - mcstart.Z());
                  mchists.fHenddx->Fill(end.X() - mcend.X());
                  mchists.fHenddy->Fill(end.Y() - mcend.Y());
                  mchists.fHenddz->Fill(end.Z() - mcend.Z());
                  mchists.fHlvsl->Fill(plen, tlen);
                  mchists.fHdl->Fill(tlen - plen);
                  mchists.fHpvsp->Fill(mcstartmom.Mag(), mom);
                  double dp = mom - mcstartmom.Mag();
                  mchists.fHdp->Fill(dp);
                  //                  mchists.fHppull->Fill(dp / std::sqrt(cov(4,4)));
                  if(std::abs(dpos) >= 5. && std::abs(dend) >= 5.) {
                    mchists.fHpvspc->Fill(mcstartmom.Mag(), mom);
                    mchists.fHdpc->Fill(dp);
                    //                    mchists.fHppullc->Fill(dp / std::sqrt(cov(4,4)));
                  }
                  // Count this track as well-reconstructed if it is matched to an
                  // mc particle (which it is if get here), and if in addition the
                  // starting position (w) matches and the reconstructed track length
                  // is more than 0.5 of the mc particle trajectory length.
                  bool good = std::abs(w) <= fWMatchDisp &&
                    tlen > 0.5 * plen;
                  if(good) {
		    mcid = part->TrackId();
		    // Calculate and fill hit efficiency and purity.
                    std::set<int> tkidset;
                    tkidset.insert(mcid);
		    /*		    double hiteff = 
		      bt->HitCollectionEfficiency(tkidset, trackhits, hitlist, geo::k3D);
		                        double hitpurity = bt->HitCollectionPurity(tkidset, trackhits);
		                        mchists.fHHitEff->Fill(hiteff);
		                        mchists.fHHitPurity->Fill(hitpurity);
		    */

		    // Fill efficiency numerator histograms.

		    mchists.fHgstartx->Fill(mcstart.X());
                    mchists.fHgstarty->Fill(mcstart.Y());
                    mchists.fHgstartz->Fill(mcstart.Z());
                    mchists.fHgendx->Fill(mcend.X());
                    mchists.fHgendy->Fill(mcend.Y());
                    mchists.fHgendz->Fill(mcend.Z());
                    mchists.fHgtheta->Fill(mcstartmom.Theta());
                    mchists.fHgphi->Fill(mcstartmom.Phi());
                    mchists.fHgtheta_xz->Fill(mctheta_xz);
                    mchists.fHgtheta_yz->Fill(mctheta_yz);
                    mchists.fHgmom->Fill(mcstartmom.Mag());
                    mchists.fHglen->Fill(plen);
                  }
                }
              }
            }

          }
          // Dump track information here.
	  /*        if(pdump) {
            TVector3 pos = track.Vertex();
            TVector3 dir = track.VertexDirection();
            TVector3 end = track.End();
            pos[0] += trackdx;
            end[0] += trackdx;
            TVector3 enddir = track.EndDirection();
            double pstart = track.VertexMomentum();
            double pend = track.EndMomentum();
            *pdump << "\nOffset"
                   << std::setw(3) << track.ID()
                   << std::setw(6) << mcid
                   << "  "
                   << std::fixed << std::setprecision(2) 
                   << std::setw(10) << trackdx
                   << "\nStart " 
                   << std::setw(3) << track.ID()
                   << std::setw(6) << mcid
                   << "  "
                   << std::fixed << std::setprecision(2) 
                   << std::setw(10) << pos[0]
                   << std::setw(10) << pos[1]
                   << std::setw(10) << pos[2];
            if(pstart > 0.) {
              *pdump << "  "
                     << std::fixed << std::setprecision(3) 
                     << std::setw(10) << dir[0]
                     << std::setw(10) << dir[1]
                     << std::setw(10) << dir[2];
            }
            else
              *pdump << std::setw(32) << " ";
            *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pstart;
            *pdump << "\nEnd   " 
                   << std::setw(3) << track.ID()
                   << std::setw(6) << mcid
                   << "  "
                   << std::fixed << std::setprecision(2)
                   << std::setw(10) << end[0]
                   << std::setw(10) << end[1]
                   << std::setw(10) << end[2];
            if(pend > 0.01) {
              *pdump << "  " 
                     << std::fixed << std::setprecision(3) 
                     << std::setw(10) << enddir[0]
                     << std::setw(10) << enddir[1]
                     << std::setw(10) << enddir[2];
            }
            else 
              *pdump << std::setw(32)<< " ";
            *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pend << "\n";
          }*/
        }
      }
  }// ipar
  
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
  fTree->Branch("trklen_cut_MC",trklen_cut_MC,"trklen_cut_MC[ntracks_reco]/D");
  fTree->Branch("trkmom_MC",trkmom_MC,"trkmom_MC[ntracks_reco]/D");
  fTree->Branch("trkmom_XMC",trkmom_XMC,"trkmom_XMC[ntracks_reco]/D");
  fTree->Branch("trkmom_YMC",trkmom_YMC,"trkmom_YMC[ntracks_reco]/D");
  fTree->Branch("trkmom_ZMC",trkmom_ZMC,"trkmom_ZMC[ntracks_reco]/D");
  fTree->Branch("trkstartdoc_XMC",trkstartdoc_XMC,"trkstartdoc_XMC[ntracks_reco]/D");
  fTree->Branch("trkstartdoc_YMC",trkstartdoc_YMC,"trkstartdoc_YMC[ntracks_reco]/D");
  fTree->Branch("trkstartdoc_ZMC",trkstartdoc_ZMC,"trkstartdoc_ZMC[ntracks_reco]/D");
  fTree->Branch("trkpdg_MC",trkpdg_MC,"trkpdg_MC[ntracks_reco]/D");
  fTree->Branch("trkd2",trkd2,"trkd2[ntracks_reco]/D");
  fTree->Branch("trkcolin",trkcolin,"trkcolin[ntracks_reco]/D");
 fTree->Branch("trktheta_xz_MC",trktheta_xz_MC,"trktheta_xz_MC[ntracks_reco]/D");
  fTree->Branch("trktheta_yz_MC",trktheta_yz_MC,"trktheta_yz_MC[ntracks_reco]/D");
  fTree->Branch("trktheta_MC",trktheta_MC,"trktheta_MC[ntracks_reco]/D");
  fTree->Branch("trkphi_MC",trkphi_MC,"trkphi_MC[ntracks_reco]/D");
  fTree->Branch("trketa_xy_MC",trketa_xy_MC,"trketa_xy_MC[ntracks_reco]/D");
  fTree->Branch("trketa_zy_MC",trketa_zy_MC,"trketa_zy_MC[ntracks_reco]/D");
  fTree->Branch("trktheta_xz",trktheta_xz,"trktheta_xz[ntracks_reco]/D");
  fTree->Branch("trktheta_yz",trktheta_yz,"trktheta_yz[ntracks_reco]/D");
  fTree->Branch("trketa_xy",trketa_xy,"trketa_xy[ntracks_reco]/D");
  fTree->Branch("trketa_zy",trketa_zy,"trketa_zy[ntracks_reco]/D");
  fTree->Branch("trktheta",trktheta,"trktheta[ntracks_reco]/D");
  fTree->Branch("trkphi",trkphi,"trkphi[ntracks_reco]/D");

  fTree->Branch("trkdedx",trkdedx,"trkdedx[ntracks_reco]/D");
  fTree->Branch("trkdedx2",trkdedx2,"trkdedx2[ntracks_reco][3][1000]/D");
  fTree->Branch("trkdqdx",trkdqdx,"trkdqdx[ntracks_reco][3][1000]/D");
  fTree->Branch("trkpitchHit",trkpitchHit,"trkpitchHit[ntracks_reco][3][1000]/D"); 
  fTree->Branch("trkkinE",trkkinE,"trkkinE[ntracks_reco][3]/D"); 
  fTree->Branch("trkplaneid",trkplaneid,"trkplaneid[ntracks_reco][3][1000]/D");
  fTree->Branch("trkresrg",trkresrg,"trkresrg[ntracks_reco][3][1000]/D");
  fTree->Branch("trkdedx_MC",trkdedx_MC,"trkdedx_MC[ntracks_reco]/D");
  fTree->Branch("trkdq_MC",trkdq_MC,"trkdq_MC[ntracks_reco]/D");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/D");
  fTree->Branch("trklen",trklen,"trklen[ntracks_reco]/D");
  fTree->Branch("trklen_L",trklen_L,"trklen_L[ntracks_reco]/D");
  fTree->Branch("trkid",trkid,"trkid[ntracks_reco]/D");
  fTree->Branch("mcang_x",&mcang_x,"mcang_x/D");
  fTree->Branch("mcang_y",&mcang_y,"mcang_y/D");
  fTree->Branch("mcang_z",&mcang_z,"mcang_z/D");

  fTree->Branch("mcpos_x",&mcpos_x,"mcpos_x/D");
  fTree->Branch("mcpos_y",&mcpos_y,"mcpos_y/D");
  fTree->Branch("mcpos_z",&mcpos_z,"mcpos_z/D");

  fTree->Branch("trkang",trkang,"trkang[ntracks_reco]/D");
  fTree->Branch("trkcolinearity",trkcolinearity,"trkcolinearity[ntracks_reco]/D");
  fTree->Branch("trkmatchdisp",trkmatchdisp,"trkmatchdisp[ntracks_reco]/D");
  fTree->Branch("trkwmatchdisp",trkwmatchdisp,"trkwmatchdisp[ntracks_reco]/D");
  fTree->Branch("trklenratio",trklenratio,"trklenratio[ntracks_reco]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/D");
  fTree->Branch("ntrkhits",ntrkhits,"ntrkhits[ntracks_reco]/I");
  fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/D");
  fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/D");
  fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/D");
  fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][3]/D");
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

AnaTree::AnaTree::RecoHists::RecoHists():
    //
    // Purpose: Default constructor.
    //
    fHstartx(0),
    fHstarty(0),
    fHstartz(0),
    fHstartd(0),
    fHendx(0),
    fHendy(0),
    fHendz(0),
    fHendd(0),
    fHtheta(0),
    fHphi(0),
    fHtheta_xz(0),
    fHtheta_yz(0),
    fHmom(0),
    fHlen(0)
    ,fHHitChg(0)
    ,fHHitWidth(0)
    ,fHHitPdg(0)
    ,fHHitTrkId(0)
    ,fModeFrac(0)
    ,fNTrkIdTrks(0)
    ,fNTrkIdTrks2(0)
    ,fNTrkIdTrks3(0)
  {}
AnaTree::AnaTree::RecoHists::RecoHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Make sure all histogram pointers are initially zero.
    *this = RecoHists();
    // Get services.
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;
    // Make histogram directory.
    art::TFileDirectory topdir = tfs->mkdir("trkana", "TrackAnaCT histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);
    // Book histograms.
    fHstartx = dir.make<TH1F>("xstart", "X Start Position",
                              100, -2.*geom->Cryostat(0).HalfWidth(), 4.*geom->Cryostat(0).HalfWidth());
    fHstarty = dir.make<TH1F>("ystart", "Y Start Position",
                              100, -geom->Cryostat(0).HalfHeight(), geom->Cryostat(0).HalfHeight());
    fHstartz = dir.make<TH1F>("zstart", "Z Start Position",
                              100, 0., geom->Cryostat(0).Length());
    fHstartd = dir.make<TH1F>("dstart", "Start Position Distance to Boundary",
                              100, -10., geom->Cryostat(0).HalfWidth());
    fHendx = dir.make<TH1F>("xend", "X End Position",
                            100, -2.*geom->Cryostat(0).HalfWidth(), 4.*geom->Cryostat(0).HalfWidth());
    fHendy = dir.make<TH1F>("yend", "Y End Position",
                            100, -geom->Cryostat(0).HalfHeight(), geom->Cryostat(0).HalfHeight());
    fHendz = dir.make<TH1F>("zend", "Z End Position",
                            100, 0., geom->Cryostat(0).Length());
    fHendd = dir.make<TH1F>("dend", "End Position Distance to Boundary",
                            100, -10., geom->Cryostat(0).HalfWidth());
    fHtheta = dir.make<TH1F>("theta", "Theta", 100, 0., 3.142);
    fHphi = dir.make<TH1F>("phi", "Phi", 100, -3.142, 3.142);
    fHtheta_xz = dir.make<TH1F>("theta_xz", "Theta_xz", 100, -3.142, 3.142);
    fHtheta_yz = dir.make<TH1F>("theta_yz", "Theta_yz", 100, -3.142, 3.142);
    fHmom = dir.make<TH1F>("mom", "Momentum", 100, 0., 10.);
    fHlen = dir.make<TH1F>("len", "Track Length", 100, 0., 3.0 * geom->Cryostat(0).Length());
    fHHitChg = dir.make<TH1F>("hchg", "Hit Charge (ADC counts)", 100, 0., 4000.);
    fHHitWidth = dir.make<TH1F>("hwid", "Hit Width (ticks)", 40, 0., 20.);
    fHHitPdg = dir.make<TH1F>("hpdg", "Hit Pdg code",5001, -2500.5, +2500.5);
    fHHitTrkId = dir.make<TH1F>("htrkid", "Hit Track ID", 401, -200.5, +200.5);
    fModeFrac = dir.make<TH1F>("hmodefrac", "quasi-Purity: Fraction of component tracks with the Track mode value", 20, 0.01, 1.01);
    fNTrkIdTrks = dir.make<TH1F>("hntrkids", "quasi-Efficiency: Number of stitched tracks in which TrkId appears", 20, 0., +10.0);
    fNTrkIdTrks2 = dir.make<TH2F>("hntrkids2", "Number of stitched tracks in which TrkId appears vs KE [GeV]", 20, 0., +10.0, 20, 0.0, 1.5);
    fNTrkIdTrks3 = dir.make<TH2F>("hntrkids3", "MC Track vs Reco Track, wtd by nhits on Collection Plane", 10, -0.5, 9.5, 10, -0.5, 9.5);
    fNTrkIdTrks3->GetXaxis()->SetTitle("Sorted-by-Descending-CPlane-Hits outer Track Number");
    fNTrkIdTrks3->GetYaxis()->SetTitle("Sorted-by-Descending-True-Length G4Track");
    
  }
  // MCHists methods.
  AnaTree::AnaTree::MCHists::MCHists() :
    //
    // Purpose: Default constructor.
    //
    fHduvcosth(0),
    fHcosth(0),
    fHmcu(0),
    fHmcv(0),
    fHmcw(0),
    fHupull(0),
    fHvpull(0),
    fHmcdudw(0),
    fHmcdvdw(0),
    fHdudwpull(0),
    fHdvdwpull(0),
    fHHitEff(0),
    fHHitPurity(0),
    fHstartdx(0),
    fHstartdy(0),
    fHstartdz(0),
    fHenddx(0),
    fHenddy(0),
    fHenddz(0),
    fHlvsl(0),
    fHdl(0),
    fHpvsp(0),
    fHpvspc(0),
    fHdp(0),
    fHdpc(0),
    fHppull(0),
    fHppullc(0),
    fHmcstartx(0),
    fHmcstarty(0),
    fHmcstartz(0),
    fHmcendx(0),
    fHmcendy(0),
    fHmcendz(0),
    fHmctheta(0),
    fHmcphi(0),
    fHmctheta_xz(0),
    fHmctheta_yz(0),
    fHmcmom(0),
    fHmclen(0),
    fHgstartx(0),
    fHgstarty(0),
    fHgstartz(0),
    fHgendx(0),
    fHgendy(0),
    fHgendz(0),
    fHgtheta(0),
    fHgphi(0),
    fHgtheta_xz(0),
    fHgtheta_yz(0),
    fHgmom(0),
    fHglen(0),
    fHestartx(0),
    fHestarty(0),
    fHestartz(0),
    fHeendx(0),
    fHeendy(0),
    fHeendz(0),
    fHetheta(0),
    fHephi(0),
    fHetheta_xz(0),
    fHetheta_yz(0),
    fHemom(0),
    fHelen(0)
  {}
AnaTree::AnaTree::MCHists::MCHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Make sure all histogram pointers are initially zero.
    *this = MCHists();
    // Get services.
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;
    // Make histogram directory.
    art::TFileDirectory topdir = tfs->mkdir("trkana", "TrackAnaCT histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);
    // Book histograms.
    fHduvcosth = dir.make<TH2F>("duvcosth", "Delta(uv) vs. Colinearity", 
                                100, 0.95, 1., 100, 0., 1.);
    fHcosth = dir.make<TH1F>("colin", "Colinearity", 100, 0.95, 1.);
    fHmcu = dir.make<TH1F>("mcu", "MC Truth U", 100, -5., 5.);
    fHmcv = dir.make<TH1F>("mcv", "MC Truth V", 100, -5., 5.);
    fHmcw = dir.make<TH1F>("mcw", "MC Truth W", 100, -20., 20.);
    fHupull = dir.make<TH1F>("dupull", "U Pull", 100, -20., 20.);
    fHvpull = dir.make<TH1F>("dvpull", "V Pull", 100, -20., 20.);
    fHmcdudw = dir.make<TH1F>("mcdudw", "MC Truth U Slope", 100, -0.2, 0.2);
    fHmcdvdw = dir.make<TH1F>("mcdvdw", "MV Truth V Slope", 100, -0.2, 0.2);
    fHdudwpull = dir.make<TH1F>("dudwpull", "U Slope Pull", 100, -10., 10.);
    fHdvdwpull = dir.make<TH1F>("dvdwpull", "V Slope Pull", 100, -10., 10.);
    fHHitEff = dir.make<TH1F>("hiteff", "MC Hit Efficiency", 100, 0., 1.0001);
    fHHitPurity = dir.make<TH1F>("hitpurity", "MC Hit Purity", 100, 0., 1.0001);
    fHstartdx = dir.make<TH1F>("startdx", "Start Delta x", 100, -10., 10.);
    fHstartdy = dir.make<TH1F>("startdy", "Start Delta y", 100, -10., 10.);
    fHstartdz = dir.make<TH1F>("startdz", "Start Delta z", 100, -10., 10.);
    fHenddx = dir.make<TH1F>("enddx", "End Delta x", 100, -10., 10.);
    fHenddy = dir.make<TH1F>("enddy", "End Delta y", 100, -10., 10.);
    fHenddz = dir.make<TH1F>("enddz", "End Delta z", 100, -10., 10.);
    fHlvsl = dir.make<TH2F>("lvsl", "Reco Length vs. MC Truth Length",
                            100, 0., 1.1 * geom->Cryostat(0).Length(), 100, 0., 1.1 * geom->Cryostat(0).Length());
    fHdl = dir.make<TH1F>("dl", "Track Length Minus MC Particle Length", 100, -50., 50.);
    fHpvsp = dir.make<TH2F>("pvsp", "Reco Momentum vs. MC Truth Momentum",
                            100, 0., 5., 100, 0., 5.);
    fHpvspc = dir.make<TH2F>("pvspc", "Reco Momentum vs. MC Truth Momentum (Contained Tracks)",
                             100, 0., 5., 100, 0., 5.);
    fHdp = dir.make<TH1F>("dp", "Reco-MC Momentum Difference", 100, -5., 5.);
    fHdpc = dir.make<TH1F>("dpc", "Reco-MC Momentum Difference (Contained Tracks)",
                           100, -5., 5.);
    fHppull = dir.make<TH1F>("ppull", "Momentum Pull", 100, -10., 10.);
    fHppullc = dir.make<TH1F>("ppullc", "Momentum Pull (Contained Tracks)", 100, -10., 10.);
    fHmcstartx = dir.make<TH1F>("mcxstart", "MC X Start Position",
                                10, -2.*geom->Cryostat(0).HalfWidth(), 4.*geom->Cryostat(0).HalfWidth());
    fHmcstarty = dir.make<TH1F>("mcystart", "MC Y Start Position",
                                10, -geom->Cryostat(0).HalfHeight(), geom->Cryostat(0).HalfHeight());
    fHmcstartz = dir.make<TH1F>("mczstart", "MC Z Start Position",
                                10, 0., geom->Cryostat(0).Length());
    fHmcendx = dir.make<TH1F>("mcxend", "MC X End Position",
                              10, -2.*geom->Cryostat(0).HalfWidth(), 4.*geom->Cryostat(0).HalfWidth());
    fHmcendy = dir.make<TH1F>("mcyend", "MC Y End Position",
                              10, -geom->Cryostat(0).HalfHeight(), geom->Cryostat(0).HalfHeight());
    fHmcendz = dir.make<TH1F>("mczend", "MC Z End Position",
                              10, 0., geom->Cryostat(0).Length());
    fHmctheta = dir.make<TH1F>("mctheta", "MC Theta", 40, 0., 3.142);
    fHmcphi = dir.make<TH1F>("mcphi", "MC Phi", 40, -3.142, 3.142);
    fHmctheta_xz = dir.make<TH1F>("mctheta_xz", "MC Theta_xz", 40, -3.142, 3.142);
    fHmctheta_yz = dir.make<TH1F>("mctheta_yz", "MC Theta_yz", 40, -3.142, 3.142);
    fHmcmom = dir.make<TH1F>("mcmom", "MC Momentum", 10, 0., 10.);
    fHmclen = dir.make<TH1F>("mclen", "MC Particle Length", 10, 0., 1.1 * geom->Cryostat(0).Length());
    fHgstartx = dir.make<TH1F>("gxstart", "Good X Start Position",
                               10, -2.*geom->Cryostat(0).HalfWidth(), 4.*geom->Cryostat(0).HalfWidth());
    fHgstarty = dir.make<TH1F>("gystart", "Good Y Start Position",
                               10, -geom->Cryostat(0).HalfHeight(), geom->Cryostat(0).HalfHeight());
    fHgstartz = dir.make<TH1F>("gzstart", "Good Z Start Position",
                               10, 0., geom->Cryostat(0).Length());
    fHgendx = dir.make<TH1F>("gxend", "Good X End Position",
                             10, -2.*geom->Cryostat(0).HalfWidth(), 4.*geom->Cryostat(0).HalfWidth());
    fHgendy = dir.make<TH1F>("gyend", "Good Y End Position",
                             10, -geom->Cryostat(0).HalfHeight(), geom->Cryostat(0).HalfHeight());
    fHgendz = dir.make<TH1F>("gzend", "Good Z End Position",
                             10, 0., geom->Cryostat(0).Length());
    fHgtheta = dir.make<TH1F>("gtheta", "Good Theta", 40, 0., 3.142);
    fHgphi = dir.make<TH1F>("gphi", "Good Phi", 40, -3.142, 3.142);
    fHgtheta_xz = dir.make<TH1F>("gtheta_xz", "Good Theta_xz", 40, -3.142, 3.142);
    fHgtheta_yz = dir.make<TH1F>("gtheta_yz", "Good Theta_yz", 40, -3.142, 3.142);
    fHgmom = dir.make<TH1F>("gmom", "Good Momentum", 10, 0., 10.);
    fHglen = dir.make<TH1F>("glen", "Good Particle Length", 10, 0., 1.1 * geom->Cryostat(0).Length());
    fHestartx = dir.make<TH1F>("exstart", "Efficiency vs. X Start Position",
                               10, -2.*geom->Cryostat(0).HalfWidth(), 4.*geom->Cryostat(0).HalfWidth());
    fHestarty = dir.make<TH1F>("eystart", "Efficiency vs. Y Start Position",
                               10, -geom->Cryostat(0).HalfHeight(), geom->Cryostat(0).HalfHeight());
    fHestartz = dir.make<TH1F>("ezstart", "Efficiency vs. Z Start Position",
                               10, 0., geom->Cryostat(0).Length());
    fHeendx = dir.make<TH1F>("exend", "Efficiency vs. X End Position",
                             10, -2.*geom->Cryostat(0).HalfWidth(), 4.*geom->Cryostat(0).HalfWidth());
    fHeendy = dir.make<TH1F>("eyend", "Efficiency vs. Y End Position",
                             10, -geom->Cryostat(0).HalfHeight(), geom->Cryostat(0).HalfHeight());
    fHeendz = dir.make<TH1F>("ezend", "Efficiency vs. Z End Position",
                             10, 0., geom->Cryostat(0).Length());
    fHetheta = dir.make<TH1F>("etheta", "Efficiency vs. Theta", 40, 0., 3.142);
    fHephi = dir.make<TH1F>("ephi", "Efficiency vs. Phi", 40, -3.142, 3.142);
    fHetheta_xz = dir.make<TH1F>("etheta_xz", "Efficiency vs. Theta_xz", 40, -3.142, 3.142);
    fHetheta_yz = dir.make<TH1F>("etheta_yz", "Efficiency vs. Theta_yz", 40, -3.142, 3.142);
    fHemom = dir.make<TH1F>("emom", "Efficiency vs. Momentum", 10, 0., 10.);
    fHelen = dir.make<TH1F>("elen", "Efficiency vs. Particle Length",
    10, 0., 1.1 * geom->Cryostat(0).Length());
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
    trklen_cut_MC[i] = -99999;
    trkmom_MC[i] = -99999;
    trkmom_XMC[i] = -99999;
    trkmom_YMC[i] = -99999;
    trkmom_ZMC[i] = -99999;
    trkstartdoc_XMC[i] = -99999;
    trkstartdoc_YMC[i] = -99999;
    trkstartdoc_ZMC[i] = -99999;
    trkpdg_MC[i] = -99999;
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
    
    for(int ii=0;ii<3;ii++)
      {
    trkkinE[i][ii] = -99999;
	for(int k=0;k<1000;k++)
	  {
	    trkdedx2[i][ii][k] = -99999;
	    trkdqdx[i][ii][k] = -99999;
	    trkpitchHit[i][ii][k] = -99999;
	    trkplaneid[i][ii][k] = -99999;
	    trkresrg[i][ii][k] = -99999;
	  }
      }
    trkdedx[i] = -99999;
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

}

void AnaTree::AnaTree::endJob()
  //
  // Purpose: End of job.
  //
  {

    // Fill efficiency histograms.
    for(std::map<int, MCHists>::const_iterator i = fMCHistMap.begin();
        i != fMCHistMap.end(); ++i) {
      const MCHists& mchists = i->second;
      effcalc(mchists.fHgstartx, mchists.fHmcstartx, mchists.fHestartx);
      effcalc(mchists.fHgstarty, mchists.fHmcstarty, mchists.fHestarty);
      effcalc(mchists.fHgstartz, mchists.fHmcstartz, mchists.fHestartz);
      effcalc(mchists.fHgendx, mchists.fHmcendx, mchists.fHeendx);
      effcalc(mchists.fHgendy, mchists.fHmcendy, mchists.fHeendy);
      effcalc(mchists.fHgendz, mchists.fHmcendz, mchists.fHeendz);
      effcalc(mchists.fHgtheta, mchists.fHmctheta, mchists.fHetheta);
      effcalc(mchists.fHgphi, mchists.fHmcphi, mchists.fHephi);
      effcalc(mchists.fHgtheta_xz, mchists.fHmctheta_xz, mchists.fHetheta_xz);
      effcalc(mchists.fHgtheta_yz, mchists.fHmctheta_yz, mchists.fHetheta_yz);
      effcalc(mchists.fHgmom, mchists.fHmcmom, mchists.fHemom);
      effcalc(mchists.fHglen, mchists.fHmclen, mchists.fHelen);
    }
  }


DEFINE_ART_MODULE(AnaTree::AnaTree)
