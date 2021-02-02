////////////////////////////////////////////////////////////////////////
// Class:       pandoraAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        pandoraAnalysis_module.cc
//
// Generated at Fri Aug  7 15:01:35 2020 by Maria Brigida Brunetti using cetskelgen
// from cetlib version v3_10_00.
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

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <TTree.h>
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>


namespace test {
  class pandoraAnalysis;
}


class test::pandoraAnalysis : public art::EDAnalyzer {
public:
  explicit pandoraAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pandoraAnalysis(pandoraAnalysis const&) = delete;
  pandoraAnalysis(pandoraAnalysis&&) = delete;
  pandoraAnalysis& operator=(pandoraAnalysis const&) = delete;
  pandoraAnalysis& operator=(pandoraAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree *fTree;

  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;

  unsigned int fNMCParticles;
  unsigned int fNPFParticles;

/*  static const int kNMaxMCParticles = 200;
  static const int kNMaxPFParticles = 200;
  static const int kNMaxPFPClusters = 3;*/

  static const int kNMaxMCParticles = 20000;
  static const int kNMaxPFParticles = 2000;
  static const int kNMaxPFPClusters = 100;
  static const int kNViews = 3;


  bool fMCIsPrimary[kNMaxMCParticles];
  int fMCParticlePdgCode[kNMaxMCParticles];
  double fMCParticleTrueEnergy [kNMaxMCParticles]; 
  int fMCParticleTrackID[kNMaxMCParticles]; 
  int fMCParticleParentTrackID[kNMaxMCParticles]; 
  std::string fMCParticleStartProcess[kNMaxMCParticles]; 
  std::string fMCParticleEndProcess[kNMaxMCParticles]; 
  int fMCParticleNTrajectoryPoint[kNMaxMCParticles]; 
  double fMCParticleStartPositionX[kNMaxMCParticles];
  double fMCParticleStartPositionY[kNMaxMCParticles];
  double fMCParticleStartPositionZ[kNMaxMCParticles];
  double fMCParticleStartPositionT[kNMaxMCParticles];
  double fMCParticleStartMomentumX[kNMaxMCParticles];
  double fMCParticleStartMomentumY[kNMaxMCParticles];
  double fMCParticleStartMomentumZ[kNMaxMCParticles];
  double fMCParticleStartMomentumE[kNMaxMCParticles];
  double fMCParticleEndPositionX[kNMaxMCParticles];
  double fMCParticleEndPositionY[kNMaxMCParticles];
  double fMCParticleEndPositionZ[kNMaxMCParticles];
  double fMCParticleEndPositionT[kNMaxMCParticles];
  double fMCParticleEndMomentumX[kNMaxMCParticles];
  double fMCParticleEndMomentumY[kNMaxMCParticles];
  double fMCParticleEndMomentumZ[kNMaxMCParticles];
  double fMCParticleEndMomentumE[kNMaxMCParticles];
  double fMCParticleVertexTime[kNMaxMCParticles];
  double fMCParticleEndTime[kNMaxMCParticles];
  int fMCParticleNHits[kNMaxMCParticles];
  //int fMCPfoMatchedPosition[kNMaxMCParticles];

  int fPFPID[kNMaxPFParticles];
  bool fPFPIsPrimary[kNMaxPFParticles];
  int fPFPTrueParticleMatchedID[kNMaxPFParticles];
  int fPFPTrueParticleMatchedPosition[kNMaxPFParticles];
  int fPFPParentID[kNMaxPFParticles];
  int fPFPPdgCode[kNMaxPFParticles];
  int fPFPNChildren[kNMaxPFParticles];
  int fPFPNClusters[kNMaxPFParticles];
  int fPFPNHits[kNMaxPFParticles];
  int fPFPNHitsView[kNMaxPFParticles][kNViews];
  int fPFPTrueParticleMatchedIDView[kNMaxPFParticles][kNViews];
  int fPFPTrueParticleMatchedPositionView[kNMaxPFParticles][kNViews];
  int fPFPNTracks[kNMaxPFParticles];
  bool fPFPIsTrack[kNMaxPFParticles];
  bool fPFPIsShower[kNMaxPFParticles];
  int fPFPCluPlane[kNMaxPFParticles][kNMaxPFPClusters];
  int fPFPCluView[kNMaxPFParticles][kNMaxPFPClusters];
  int fPFPCluNHits[kNMaxPFParticles][kNMaxPFPClusters];
  double fPFPCluIntegral[kNMaxPFParticles][kNMaxPFPClusters];


  int fPFPTrackID[kNMaxPFParticles];
  double fPFPTrackLength[kNMaxPFParticles];
  double fPFPTrackStartX[kNMaxPFParticles];
  double fPFPTrackStartY[kNMaxPFParticles];
  double fPFPTrackStartZ[kNMaxPFParticles];
  double fPFPTrackVertexX[kNMaxPFParticles];
  double fPFPTrackVertexY[kNMaxPFParticles];
  double fPFPTrackVertexZ[kNMaxPFParticles];
  double fPFPTrackEndX[kNMaxPFParticles];
  double fPFPTrackEndY[kNMaxPFParticles];
  double fPFPTrackEndZ[kNMaxPFParticles];
  double fPFPTrackTheta[kNMaxPFParticles];
  double fPFPTrackPhi[kNMaxPFParticles];
  double fPFPTrackZenithAngle[kNMaxPFParticles];
  double fPFPTrackAzimuthAngle[kNMaxPFParticles];
  double fPFPTrackStartDirectionX[kNMaxPFParticles];
  double fPFPTrackStartDirectionY[kNMaxPFParticles];
  double fPFPTrackStartDirectionZ[kNMaxPFParticles];
  double fPFPTrackVertexDirectionX[kNMaxPFParticles];
  double fPFPTrackVertexDirectionY[kNMaxPFParticles];
  double fPFPTrackVertexDirectionZ[kNMaxPFParticles];
  double fPFPTrackEndDirectionX[kNMaxPFParticles];
  double fPFPTrackEndDirectionY[kNMaxPFParticles];
  double fPFPTrackEndDirectionZ[kNMaxPFParticles];
  float fPFPTrackChi2[kNMaxPFParticles];
  int fPFPTrackNdof[kNMaxPFParticles];
 
  int    fPFPShowerID[kNMaxPFParticles];
  int    fPFPShowerBestPlane[kNMaxPFParticles];
  double fPFPShowerDirectionX[kNMaxPFParticles];
  double fPFPShowerDirectionY[kNMaxPFParticles];
  double fPFPShowerDirectionZ[kNMaxPFParticles];
  double fPFPShowerDirectionErrX[kNMaxPFParticles];
  double fPFPShowerDirectionErrY[kNMaxPFParticles];
  double fPFPShowerDirectionErrZ[kNMaxPFParticles];
  double fPFPShowerStartX[kNMaxPFParticles];
  double fPFPShowerStartY[kNMaxPFParticles];
  double fPFPShowerStartZ[kNMaxPFParticles];
  double fPFPShowerStartErrX[kNMaxPFParticles];
  double fPFPShowerStartErrY[kNMaxPFParticles];
  double fPFPShowerStartErrZ[kNMaxPFParticles];
  double fPFPShowerLength[kNMaxPFParticles];
  double fPFPShowerOpenAngle[kNMaxPFParticles];

  double fPFPCompleteness[kNMaxMCParticles];
  double fPFPPurity[kNMaxMCParticles];

  std::string fTruthLabel;
  std::string fHitLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPFParticleLabel;
  bool fRollUpUnsavedIDs;

  const geo::Geometry* fGeom;
};


test::pandoraAnalysis::pandoraAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fTruthLabel = p.get<std::string>("TruthLabel");
  fHitLabel = p.get<std::string>("HitLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fShowerLabel = p.get<std::string>("ShowerLabel");
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
  fRollUpUnsavedIDs = p.get<bool>("RollUpUnsavedIDs"); 
} 

void test::pandoraAnalysis::analyze(art::Event const& e)
{
  const art::ServiceHandle<cheat::BackTrackerService> btServ;
  fEventID = e.id().event();
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  fNMCParticles = 0;
  fNPFParticles = 0;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  std::cout << "=============== EVENT ID " << fEventID << " == RUN ID " << fRunID << " == SUBRUN ID " << fSubRunID << " ================" << std::endl;

  for(int iMc=0; iMc<kNMaxMCParticles; iMc++){
    fMCIsPrimary[iMc]=0;
    fMCParticlePdgCode[iMc]=0;
    fMCParticleTrueEnergy[iMc]=0;
    fMCParticleTrackID[iMc]=999999;
    fMCParticleParentTrackID[iMc]=0;
    fMCParticleStartProcess[iMc]="";
    fMCParticleEndProcess[iMc]="";
    fMCParticleNTrajectoryPoint[iMc]=0;
    fMCParticleStartPositionX[iMc]=0;
    fMCParticleStartPositionY[iMc]=0;
    fMCParticleStartPositionZ[iMc]=0;
    fMCParticleStartPositionT[iMc]=0;
    fMCParticleStartMomentumX[iMc]=0;
    fMCParticleStartMomentumY[iMc]=0;
    fMCParticleStartMomentumZ[iMc]=0;
    fMCParticleStartMomentumE[iMc]=0;
    fMCParticleEndPositionX[iMc]=0;
    fMCParticleEndPositionY[iMc]=0;
    fMCParticleEndPositionZ[iMc]=0;
    fMCParticleEndPositionT[iMc]=0;
    fMCParticleEndMomentumX[iMc]=0;
    fMCParticleEndMomentumY[iMc]=0;
    fMCParticleEndMomentumZ[iMc]=0;
    fMCParticleEndMomentumE[iMc]=0;
    fMCParticleVertexTime[iMc]=0;
    fMCParticleEndTime[iMc]=0;
    fMCParticleNHits[iMc]=0;
    //fMCPfoMatchedPosition[iMc]=999999;
  }

  for(int iPfp=0; iPfp<kNMaxPFParticles; iPfp++){

   fPFPID[iPfp]=0;
   fPFPTrueParticleMatchedID[iPfp]=999999;
   fPFPTrueParticleMatchedPosition[iPfp]=999999;

   fPFPNHits[iPfp]=0;
   for (unsigned int iView = 0; iView < kNViews; iView++)
   {
        fPFPNHitsView[iPfp][iView] = 999999;
        fPFPTrueParticleMatchedIDView[iPfp][iView] = 999999;
        fPFPTrueParticleMatchedPositionView[iPfp][iView] = 999999;
   }

   fPFPNClusters[iPfp]=0;
   fPFPIsTrack[iPfp]=0;
   fPFPIsShower[iPfp]=0;

   fPFPTrackID[iPfp]=0;
   fPFPTrackLength[iPfp]=0;
   fPFPTrackStartX[iPfp]=0;
   fPFPTrackStartY[iPfp]=0;
   fPFPTrackStartZ[iPfp]=0;
   fPFPTrackVertexX[iPfp]=0;
   fPFPTrackVertexY[iPfp]=0;
   fPFPTrackVertexZ[iPfp]=0;
   fPFPTrackEndX[iPfp]=0;
   fPFPTrackEndY[iPfp]=0;
   fPFPTrackEndZ[iPfp]=0;
   fPFPTrackTheta[iPfp]=0;
   fPFPTrackPhi[iPfp]=0;
   fPFPTrackZenithAngle[iPfp]=0;
   fPFPTrackAzimuthAngle[iPfp]=0;
   fPFPTrackStartDirectionX[iPfp]=0;
   fPFPTrackStartDirectionY[iPfp]=0;
   fPFPTrackStartDirectionZ[iPfp]=0;
   fPFPTrackVertexDirectionX[iPfp]=0;
   fPFPTrackVertexDirectionY[iPfp]=0;
   fPFPTrackVertexDirectionZ[iPfp]=0;
   fPFPTrackEndDirectionX[iPfp]=0;
   fPFPTrackEndDirectionY[iPfp]=0;
   fPFPTrackEndDirectionZ[iPfp]=0;
   fPFPTrackChi2[iPfp]=0;
   fPFPTrackNdof[iPfp]=0;
   
   fPFPShowerID[iPfp]=0;
   fPFPShowerBestPlane[iPfp]=0;
   fPFPShowerDirectionX[iPfp]=0;
   fPFPShowerDirectionY[iPfp]=0;
   fPFPShowerDirectionZ[iPfp]=0;
   fPFPShowerDirectionErrX[iPfp]=0;
   fPFPShowerDirectionErrY[iPfp]=0;
   fPFPShowerDirectionErrZ[iPfp]=0;
   fPFPShowerStartX[iPfp]=0;
   fPFPShowerStartY[iPfp]=0;
   fPFPShowerStartZ[iPfp]=0;
   fPFPShowerStartErrX[iPfp]=0;
   fPFPShowerStartErrY[iPfp]=0;
   fPFPShowerStartErrZ[iPfp]=0;
   fPFPShowerLength[iPfp]=0;
   fPFPShowerOpenAngle[iPfp]=0;

   for(int iClu=0; iClu<kNMaxPFPClusters; iClu++){
     fPFPCluPlane[iPfp][iClu]=0; 
     fPFPCluView[iPfp][iClu]=0; 
     fPFPCluNHits[iPfp][iClu]=0; 
     fPFPCluIntegral[iPfp][iClu]=0; 
   }

   fPFPCompleteness[iPfp]=0;
   fPFPPurity[iPfp]=0;

  }

  //Get all hits
  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit> > allHits;
  if(e.getByLabel(fHitLabel,hitHandle))
  {art::fill_ptr_vector(allHits, hitHandle);}

  //Fill MC particle to hits map
  std::map<int,int> trueParticleHits;
  for (const auto& hit: allHits){
      TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
      if (TruthMatchUtils::Valid(g4ID))
          trueParticleHits[g4ID]++;
  }

  //Access the truth information
  if(!e.isRealData()){
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    if(mcParticles.isValid()){
 
     fNMCParticles=mcParticles->size();
      bool isMCPrimary(false);
      for(unsigned int iMc=0; iMc< mcParticles->size(); iMc++){
        const simb::MCParticle trueParticle = mcParticles->at(iMc);
        fMCParticleTrueEnergy[iMc]=trueParticle.E();
        fMCParticlePdgCode[iMc]=trueParticle.PdgCode();
        fMCParticleTrackID[iMc]=trueParticle.TrackId();
        fMCParticleVertexTime[iMc]=trueParticle.T();
        fMCParticleEndTime[iMc]=trueParticle.EndT();
        fMCParticleParentTrackID[iMc]=trueParticle.Mother();
        fMCParticleStartProcess[iMc]=trueParticle.Process();
        fMCParticleEndProcess[iMc]=trueParticle.EndProcess();
        fMCParticleNTrajectoryPoint[iMc]=trueParticle.NumberTrajectoryPoints();
        trueParticle.Process()=="primary"? isMCPrimary=true:isMCPrimary=false;
        fMCIsPrimary[iMc] = isMCPrimary;
        fMCParticleNHits[iMc]=trueParticleHits[trueParticle.TrackId()];
        fMCParticleStartPositionX[iMc] = trueParticle.Position().X();
        fMCParticleStartPositionY[iMc] = trueParticle.Position().Y();
        fMCParticleStartPositionZ[iMc] = trueParticle.Position().Z();
        fMCParticleStartPositionT[iMc] = trueParticle.Position().T();
        fMCParticleEndPositionX[iMc] = trueParticle.EndPosition().X();
        fMCParticleEndPositionY[iMc] = trueParticle.EndPosition().Y();
        fMCParticleEndPositionZ[iMc] = trueParticle.EndPosition().Z();
        fMCParticleEndPositionT[iMc] = trueParticle.EndPosition().T();

        fMCParticleStartMomentumX[iMc] = trueParticle.Momentum().X();
        fMCParticleStartMomentumY[iMc] = trueParticle.Momentum().Y();
        fMCParticleStartMomentumZ[iMc] = trueParticle.Momentum().Z();
        fMCParticleStartMomentumE[iMc] = trueParticle.Momentum().E();
        fMCParticleEndMomentumX[iMc] = trueParticle.EndMomentum().X();
        fMCParticleEndMomentumY[iMc] = trueParticle.EndMomentum().Y();
        fMCParticleEndMomentumZ[iMc] = trueParticle.EndMomentum().Z();
        fMCParticleEndMomentumE[iMc] = trueParticle.EndMomentum().E();

      }
    }
  }

  //Access the PFParticles from Pandora
  const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);
  fNPFParticles = pfparticleVect.size();
  if(!fNPFParticles) {
    std::cout << "No PFParticles found!" << std::endl;
    return;
  }

  //Access the Clusters
  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  std::vector<art::Ptr<recob::Cluster>> clusterVect;
  if (e.getByLabel(fPFParticleLabel,clusterHandle))
    art::fill_ptr_vector(clusterVect,clusterHandle);

  art::FindManyP<recob::Cluster> clusterParticleAssoc(pfparticleVect, e, fPFParticleLabel);

  art::Handle<std::vector<recob::Track>> trackHandle;
  if (!(e.getByLabel(fTrackLabel, trackHandle))){
    std::cout<<"Unable to find std::vector<recob::Track> with module label: " << fTrackLabel << std::endl;
    return;
  }

  std::vector<art::Ptr<recob::Track> > trackList;
  if (e.getByLabel(fTrackLabel,trackHandle)){
    art::fill_ptr_vector(trackList, trackHandle);
  }

  //std::map<int,int> pfpToMcMap;
  int iPfp(0);
  for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){

    fPFPID[iPfp]=iPfp;
    fPFPIsPrimary[iPfp]=pfp->IsPrimary();
    fPFPPdgCode[iPfp]=pfp->PdgCode();
    fPFPNChildren[iPfp]=pfp->NumDaughters();
    (pfp->IsPrimary())? fPFPParentID[iPfp]=-1 : fPFPParentID[iPfp]=pfp->Parent();

    std::vector<art::Ptr<recob::Cluster>> pfpClusters = clusterParticleAssoc.at(pfp.key());
    fPFPNClusters[iPfp]=pfpClusters.size();
    if(!pfpClusters.empty()){
      int iClu(0);
      for(const art::Ptr<recob::Cluster> &clu:pfpClusters){
	fPFPCluPlane[iPfp][iClu]=clu->Plane().asPlaneID().Plane;
	fPFPCluView[iPfp][iClu]=clu->View();
        fPFPCluNHits[iPfp][iClu]=clu->NHits();
	fPFPCluIntegral[iPfp][iClu]=clu->Integral();
        iClu++;
        if (iClu == kNMaxPFPClusters)
            break;
      }
    }

    std::vector<art::Ptr<recob::Hit>> pfpHits;

    if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp,e,fPFParticleLabel,fTrackLabel)){
      fPFPIsTrack[iPfp]=true;
      art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp,e,fPFParticleLabel, fTrackLabel); 

      fPFPTrackID[iPfp]=track->ID();
      fPFPTrackLength[iPfp]=track->Length();
      fPFPTrackStartX[iPfp]=track->Start().X();
      fPFPTrackStartY[iPfp]=track->Start().Y();
      fPFPTrackStartZ[iPfp]=track->Start().Z();
      fPFPTrackVertexX[iPfp]=track->Vertex().X();
      fPFPTrackVertexY[iPfp]=track->Vertex().Y();
      fPFPTrackVertexZ[iPfp]=track->Vertex().Z();
      fPFPTrackEndX[iPfp]=track->End().X();
      fPFPTrackEndY[iPfp]=track->End().Y();
      fPFPTrackEndZ[iPfp]=track->End().Z();
      fPFPTrackTheta[iPfp]=track->Theta();
      fPFPTrackPhi[iPfp]=track->Phi();
      fPFPTrackZenithAngle[iPfp]=track->ZenithAngle();
      fPFPTrackAzimuthAngle[iPfp]=track->AzimuthAngle();
      fPFPTrackStartDirectionX[iPfp]=track->StartDirection().X();
      fPFPTrackStartDirectionY[iPfp]=track->StartDirection().Y();
      fPFPTrackStartDirectionZ[iPfp]=track->StartDirection().Z();
      fPFPTrackVertexDirectionX[iPfp]=track->VertexDirection().X();
      fPFPTrackVertexDirectionY[iPfp]=track->VertexDirection().Y();
      fPFPTrackVertexDirectionZ[iPfp]=track->VertexDirection().Z();
      fPFPTrackEndDirectionX[iPfp]=track->EndDirection().X();
      fPFPTrackEndDirectionY[iPfp]=track->EndDirection().Y();
      fPFPTrackEndDirectionZ[iPfp]=track->EndDirection().Z();
      fPFPTrackChi2[iPfp]=track->Chi2();
      fPFPTrackNdof[iPfp]=track->Ndof();

      pfpHits = dune_ana::DUNEAnaTrackUtils::GetHits(track,e,fTrackLabel);
      
      fPFPNHits[iPfp]=pfpHits.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 0);
      fPFPNHitsView[iPfp][0] = pfpHitsView0.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 1);
      fPFPNHitsView[iPfp][1] = pfpHitsView1.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 2);
      fPFPNHitsView[iPfp][2] = pfpHitsView2.size();

      if(!e.isRealData()) {
        TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHits,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID)){
          fPFPTrueParticleMatchedID[iPfp] = g4ID;

	  int pos(999999); for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) {
	    if(fMCParticleTrackID[ipos]==g4ID)pos=ipos;
 	  }
          fPFPTrueParticleMatchedPosition[iPfp] = pos;

        }
        TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView0,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][0] = g4IDView0;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView0)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][0] = ipos;
                break;
            }
        }
        TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView1,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][1] = g4IDView1;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView1)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][1] = ipos;
                break;
            }
        }
        TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView2,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][2] = g4IDView2;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView2)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][2] = ipos;
                break;
            }
        }

      }
    }

    if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp,e,fPFParticleLabel,fShowerLabel)){
      fPFPIsShower[iPfp]=true;
      art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp,e,fPFParticleLabel, fShowerLabel); 
      fPFPShowerID[iPfp]=shower->ID();
      fPFPShowerBestPlane[iPfp]=shower->best_plane();
      fPFPShowerDirectionX[iPfp]=shower->Direction().X();
      fPFPShowerDirectionY[iPfp]=shower->Direction().Y();
      fPFPShowerDirectionZ[iPfp]=shower->Direction().Z();
      fPFPShowerDirectionErrX[iPfp]=shower->DirectionErr().X();
      fPFPShowerDirectionErrY[iPfp]=shower->DirectionErr().Y();
      fPFPShowerDirectionErrZ[iPfp]=shower->DirectionErr().Z();
      fPFPShowerStartX[iPfp]=shower->ShowerStart().X();
      fPFPShowerStartY[iPfp]=shower->ShowerStart().Y();
      fPFPShowerStartZ[iPfp]=shower->ShowerStart().Z();
      fPFPShowerStartErrX[iPfp]=shower->ShowerStartErr().X();
      fPFPShowerStartErrY[iPfp]=shower->ShowerStartErr().Y();
      fPFPShowerStartErrZ[iPfp]=shower->ShowerStartErr().Z();
      fPFPShowerLength[iPfp]=shower->Length();
      fPFPShowerOpenAngle[iPfp]=shower->OpenAngle();

      pfpHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower,e,fShowerLabel);
      fPFPNHits[iPfp]=pfpHits.size();

      if(!e.isRealData()) {
        TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHits,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID)){
          fPFPTrueParticleMatchedID[iPfp] = g4ID;
	  int pos(999999); for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) {
	    if(fMCParticleTrackID[ipos]==g4ID)pos=ipos;
  	  }
          fPFPTrueParticleMatchedPosition[iPfp] = pos;
        }
      }
    }

    if(!e.isRealData()) {
      int pfpNHits = fPFPNHits[iPfp];
      int nHitsBestMatchTrueParticle = fMCParticleNHits[fPFPTrueParticleMatchedPosition[iPfp]]; 
      //int nSharedHits = 
      //std::map<int, int> pfpTrueHitsMap = GetTruePrimaryHits(clockData, trueParticles, truePrimaries, pfpHits);
      //int   pfpHitsTrueHits    = pfpTrueHitsMap.at(trueId.first);

      //Fill shared MC particle to hits map for this specific Pfp
      int nSharedTrueParticlePfpHits(0);
      for (const auto& hit: pfpHits){
        int trackID     = 0;
        float hitEnergy = 0;
        std::vector<sim::TrackIDE> trackIDEs = btServ->HitToTrackIDEs(clockData, hit);
        for (const auto& ide: trackIDEs) {
          if (ide.energy > hitEnergy){
            hitEnergy = ide.energy;
            trackID   = TMath::Abs(ide.trackID);
          }
        }
        if(trackID==fPFPTrueParticleMatchedID[iPfp])++nSharedTrueParticlePfpHits;
      }
      fPFPCompleteness[iPfp]     = (float)nSharedTrueParticlePfpHits / pfpNHits;
      fPFPPurity[iPfp]           = (float)nSharedTrueParticlePfpHits / nHitsBestMatchTrueParticle;
    }
    iPfp++;
  }
  fTree->Fill();
}

void test::pandoraAnalysis::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("pandoraOutput","Pandora Output Tree");

  //Event branches
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("runID",&fRunID,"runID/i");
  fTree->Branch("subrunID",&fSubRunID,"subrunID/i");
  fTree->Branch("nMCParticles",&fNMCParticles,"nMCParticles/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");

  //MC truth branches
  //fTree->Branch("mcIsMCPrimary",&fMCIsPrimary);
  fTree->Branch("mcIsMCPrimary",&fMCIsPrimary,"MCIsPrimary[nMCParticles]/O");
  fTree->Branch("mcParticlePdgCode",&fMCParticlePdgCode,"MCParticlePdgCode[nMCParticles]/I");
  fTree->Branch("mcParticleTrueEnergy",&fMCParticleTrueEnergy,"MCParticleTrueEnergy[nMCParticles]/D");
  fTree->Branch("mcParticleTrackID",&fMCParticleTrackID,"MCParticleTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleParentTrackID",&fMCParticleParentTrackID,"MCParticleParentTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleStartProcess",&fMCParticleStartProcess,"MCParticleStartProcess[nMCParticles]/C");
  fTree->Branch("mcParticleEndProcess",&fMCParticleEndProcess, "MCParticleEndProcess[nMCParticles]/C");
  fTree->Branch("mcParticleNTrajectoryPoints",&fMCParticleNTrajectoryPoint, "MCParticleNTrajectoryPoint[nMCParticles]/I");
  fTree->Branch("mcParticleStartPositionX",&fMCParticleStartPositionX,"MCParticleStartPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionY",&fMCParticleStartPositionY,"MCParticleStartPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionZ",&fMCParticleStartPositionZ,"MCParticleStartPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionT",&fMCParticleStartPositionT,"MCParticleStartPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumX",&fMCParticleStartMomentumX, "MCParticleStartMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumY",&fMCParticleStartMomentumY, "MCParticleStartMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumZ",&fMCParticleStartMomentumZ, "MCParticleStartMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumE",&fMCParticleStartMomentumE, "MCParticleStartMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionX",&fMCParticleEndPositionX, "MCParticleEndPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionY",&fMCParticleEndPositionY, "MCParticleEndPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionZ",&fMCParticleEndPositionZ, "MCParticleEndPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionT",&fMCParticleEndPositionT, "MCParticleEndPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumX",&fMCParticleEndMomentumX, "MCParticleEndMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumY",&fMCParticleEndMomentumY, "MCParticleEndMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumZ",&fMCParticleEndMomentumZ, "MCParticleEndMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumE",&fMCParticleEndMomentumE, "MCParticleEndMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleVertexTime",&fMCParticleVertexTime, "MCParticleVertexTime[nMCParticles]/D");
  fTree->Branch("mcParticleEndTime",&fMCParticleEndTime, "MCParticleEndTime[nMCParticles]/D");
  fTree->Branch("mcParticleNHits", &fMCParticleNHits, "MCParticleNHits[nMCParticles]/I");
  //fTree->Branch("mcPfoMatchedPosition", &fMCPfoMatchedPosition, "fMCPfoMatchedPosition[nMCParticles]/I");

  //PFP branches
  fTree->Branch("pfpTrueParticleMatchedID",&fPFPTrueParticleMatchedID,"PFPTrueParticleMatchedID[nPFParticles]/I");
  fTree->Branch("pfpTrueParticleMatchedPosition",&fPFPTrueParticleMatchedPosition,"PFPTrueParticleMatchedPosition[nPFParticles]/I");
  fTree->Branch("pfpIsPrimary",&fPFPIsPrimary,"PFPIsPrimary[nPFParticles]/O");
  fTree->Branch("pfpID",&fPFPID, "PFPID[nPFParticles]/I");
  fTree->Branch("pfpParentID",&fPFPParentID, "PFPParentID[nPFParticles]/I");
  fTree->Branch("pfpPdgCode",&fPFPPdgCode, "PFPPdgCode[nPFParticles]/I");
  fTree->Branch("pfpNChildren",&fPFPNChildren,"PFPNChildren[nPFParticles]/I");
  fTree->Branch("pfpNClusters",&fPFPNClusters,"PFPNClusters[nPFParticles]/I");
  fTree->Branch("pfpNHits",&fPFPNHits,"PFPNHits[nPFParticles]/I");
  fTree->Branch("pfpNHitsView",&fPFPNHitsView,"PFPNHitsView[nPFParticles][3]/I");
  fTree->Branch("pfpTrueParticleMatchedIDView",&fPFPTrueParticleMatchedIDView,"PFPTrueParticleMatchedIDView[nPFParticles][3]/I");
  fTree->Branch("pfpTrueParticleMatchedPositionView",&fPFPTrueParticleMatchedPositionView,"PFPTrueParticleMatchedPositionView[nPFParticles][3]/I");
  fTree->Branch("pfpNTracks",&fPFPNTracks,"PFPNTracks[nPFParticles]/I");
  fTree->Branch("pfpIsTrack",&fPFPIsTrack,"PFPIsTrack[nPFParticles]/O");
  fTree->Branch("pfpIsShower",&fPFPIsShower,"PFPIsShower[nPFParticles]/O");
  fTree->Branch("pfpTrackID", &fPFPTrackID,"PFPNClusters[nPFParticles]/I");
  fTree->Branch("pfpTrackLength",&fPFPTrackLength,"PFPTrackLength[nPFParticles]/D");
  fTree->Branch("pfpTrackStartX",&fPFPTrackStartX,"PFPTrackStartX[nPFParticles]/D");
  fTree->Branch("pfpTrackStartY",&fPFPTrackStartY,"PFPTrackStartY[nPFParticles]/D");
  fTree->Branch("pfpTrackStartZ",&fPFPTrackStartZ,"PFPTrackStartZ[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexX",&fPFPTrackVertexX,"PFPTrackVertexX[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexY",&fPFPTrackVertexY,"PFPTrackVertexY[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexZ",&fPFPTrackVertexZ,"PFPTrackVertexZ[nPFParticles]/D");
  fTree->Branch("pfpTrackEndX",&fPFPTrackEndX,"PFPTrackEndX[nPFParticles]/D");
  fTree->Branch("pfpTrackEndY",&fPFPTrackEndY,"PFPTrackEndY[nPFParticles]/D");
  fTree->Branch("pfpTrackEndZ",&fPFPTrackEndZ,"PFPTrackEndZ[nPFParticles]/D");
  fTree->Branch("pfpTrackTheta",&fPFPTrackTheta,"PFPTrackTheta[nPFParticles]/D");
  fTree->Branch("pfpTrackPhi",&fPFPTrackPhi,"PFPTrackPhi[nPFParticles]/D");
  fTree->Branch("pfpTrackZenithAngle",&fPFPTrackZenithAngle,"PFPTrackZenithAngle[nPFParticles]/D");
  fTree->Branch("pfpTrackAzimuthAngle",&fPFPTrackAzimuthAngle,"PFPTrackAzimuthAngle[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionX",&fPFPTrackStartDirectionX,"PFPTrackStartDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionY",&fPFPTrackStartDirectionY,"PFPTrackStartDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionZ",&fPFPTrackStartDirectionZ,"PFPTrackStartDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionX",&fPFPTrackVertexDirectionX,"PFPTrackVertexDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionY",&fPFPTrackVertexDirectionY,"PFPTrackVertexDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionZ",&fPFPTrackVertexDirectionZ,"PFPTrackVertexDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionX",&fPFPTrackEndDirectionX,"PFPTrackEndDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionY",&fPFPTrackEndDirectionY,"PFPTrackEndDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionZ",&fPFPTrackEndDirectionZ,"PFPTrackEndDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackChi2",&fPFPTrackChi2,"PFPTrackChi2[nPFParticles]/F");
  fTree->Branch("pfpTrackStartNdof",&fPFPTrackNdof,"PFPTrackNdof[nPFParticles]/I");

  fTree->Branch("pfpCluPlane",fPFPCluPlane,"PFPCluPlane[nPFParticles][100]/I");
  fTree->Branch("pfpCluView",fPFPCluView,"PFPCluView[nPFParticles][100]/I");
  fTree->Branch("pfpCluNHits",fPFPCluNHits,"PFPCluNHits[nPFParticles][100]/I");
  fTree->Branch("pfpCluIntegral",fPFPCluIntegral,"PFPCluIntegral[nPFParticles][100]/D");

  fTree->Branch("pfpShowerID",&fPFPShowerID,"PFPShowerID[nPFParticles]/I");
  fTree->Branch("pfpShowerBestPlane",&fPFPShowerBestPlane,"PFPShowerBestPlane[nPFParticles]/I");
  fTree->Branch("pfpShowerDirectionX",&fPFPShowerDirectionX,"PFPShowerDirectionX[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionY",&fPFPShowerDirectionY,"PFPShowerDirectionY[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionZ",&fPFPShowerDirectionZ,"PFPShowerDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrX",&fPFPShowerDirectionErrX,"PFPShowerDirectionErrX[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrY",&fPFPShowerDirectionErrY,"PFPShowerDirectionErrY[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrZ",&fPFPShowerDirectionErrZ,"PFPShowerDirectionErrZ[nPFParticles]/D");
  fTree->Branch("pfpShowerStartX",&fPFPShowerStartX,"PFPShowerStartX[nPFParticles]/D");
  fTree->Branch("pfpShowerStartY",&fPFPShowerStartY,"PFPShowerStartY[nPFParticles]/D");
  fTree->Branch("pfpShowerStartZ",&fPFPShowerStartZ,"PFPShowerStartZ[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrX",&fPFPShowerStartErrX,"PFPShowerStartErrX[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrY",&fPFPShowerStartErrY,"PFPShowerStartErrY[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrZ",&fPFPShowerStartErrZ,"PFPShowerStartErrZ[nPFParticles]/D");
  fTree->Branch("pfpShowerLength",&fPFPShowerLength,"PFPShowerLength[nPFParticles]/D");
  fTree->Branch("pfpShowerOpeningAngle",&fPFPShowerOpenAngle,"PFPShowerOpenAngle[nPFParticles]/D");

  fTree->Branch("pfpCompleteness", &fPFPCompleteness, "PFPCompleteness[nPFParticles]/D");
  fTree->Branch("pfpPurity", &fPFPPurity, "PFPPurity[nPFParticles]/D");
}
void test::pandoraAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::pandoraAnalysis)
