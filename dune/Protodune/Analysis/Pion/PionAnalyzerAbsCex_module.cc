////////////////////////////////////////////////////////////////////////
// Class:       PionAnalyzerAbsCex
// Plugin Type: analyzer (art v3_00_00)
// File:        PionAnalyzerAbsCex_module.cc
//
// Generated at Tue Jan  8 09:12:19 2019 by Jacob Calcutt using cetskelgen
// from cetlib version v3_04_00.
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

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEShowerUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "art/Framework/Services/Optional/TFileService.h"

// ROOT includes
#include "TTree.h"

namespace pionana {
  class PionAnalyzerAbsCex;
}


class pionana::PionAnalyzerAbsCex : public art::EDAnalyzer {
public:
  explicit PionAnalyzerAbsCex(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PionAnalyzerAbsCex(PionAnalyzerAbsCex const&) = delete;
  PionAnalyzerAbsCex(PionAnalyzerAbsCex&&) = delete;
  PionAnalyzerAbsCex& operator=(PionAnalyzerAbsCex const&) = delete;
  PionAnalyzerAbsCex& operator=(PionAnalyzerAbsCex&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();

private:

  // Declare member data here.
  //const art::InputTag fSpacePointModuleLabel;
  const art::InputTag fBeamModuleLabel;
  const art::InputTag fTrackModuleLabel;

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;

  // beam information
  std::vector<double> beamPosx;
  std::vector<double> beamPosy;
  std::vector<double> beamPosz;
  
  std::vector<double> beamDirx;
  std::vector<double> beamDiry;
  std::vector<double> beamDirz;

  std::vector<double> beamMomentum;

  double tof;
  short ckov0status;
  short ckov1status;

  // fcl parameters for PFP particles
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  bool fVerbose;
  protoana::ProtoDUNEDataUtils dataUtil;

  // define parameters for primary tracks
  std::vector<double> primtrk_startx;
  std::vector<double> primtrk_starty;
  std::vector<double> primtrk_startz;
  std::vector<double> primtrk_Dirx;
  std::vector<double> primtrk_Diry;
  std::vector<double> primtrk_Dirz;
  std::vector<double> primtrklen;
  std::vector<double> primtrkID;
  std::vector<int> primtrk_trktag;

  double cosine_beam_primtrk;
 
  std::vector<int> pdg_code;
  std::vector<int> n_daughter;
  std::vector<int> n_beamparticle;
  std::vector<int> isPrimary;
  std::vector<int> pfp_self;
  //std::vector<int> pfp_parent;
  std::vector<int> pfp_daughter;
  std::vector<int> pfp_daughter_PDGs;
  std::vector<int> pfp_daughter_TrackIDs;
  std::vector<int> pfp_daughter_ShowerIDs;

};


pionana::PionAnalyzerAbsCex::PionAnalyzerAbsCex(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),

  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fVerbose(p.get<bool>("Verbose")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pionana::PionAnalyzerAbsCex::analyze(art::Event const& evt)
{
  //reset containers
  reset();  


//  art::Handle< std::vector<recob::Track> > trackListHandle;
//  std::vector<art::Ptr<recob::Track> > tracklist;
//  if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
//  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryTag);
//  art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle, evt, "pandoraTrack");


  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

//  art::Handle< std::vector<recob::Track> > trkHandle;
//  std::vector< art::Ptr<recob::Track> > trks;
//  if (evt.getByLabel(fTrackModuleLabel, trkHandle))
//    art::fill_ptr_vector(trks, trkHandle);


  art::Handle< std::vector<beam::ProtoDUNEBeamEvent> > pdbeamHandle;
  std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
  if (evt.getByLabel(fBeamModuleLabel, pdbeamHandle))
    art::fill_ptr_vector(beaminfo, pdbeamHandle);
  else{
    std::cout<<"No beam information from "<<fBeamModuleLabel<<std::endl;
  }

  tof = -1;
  ckov0status = -1;
  ckov1status = -1;

  if (beaminfo.size()){
     
    auto const & beamevent = beaminfo[0];

    if (beamevent->GetTimingTrigger() == 12 && beamevent->CheckIsMatched()){

      //Get TOF info if valid
      if (beamevent->GetTOFChan() != -1){
        tof = beamevent->GetTOF();
      }
      
      //JC: IMPORTANT! Try filtering out for events with Momentum/Track vector size == 1

      //Get beam particle trajectory info
      auto & tracks = beamevent->GetBeamTracks();
      std::cout<<"###############################################################"<<std::endl;
      std::cout<<"ToF:"<<tof<<" [ns]"<<std::endl;
      std::cout<<"beam trk size:"<<tracks.size()<<std::endl;
      for (size_t i = 0; i<tracks.size(); ++i){
        beamPosx.push_back(tracks[i].End().X());
        beamPosy.push_back(tracks[i].End().Y());
        beamPosz.push_back(tracks[i].End().Z());
        beamDirx.push_back(tracks[i].StartDirection().X());
        beamDiry.push_back(tracks[i].StartDirection().Y());
        beamDirz.push_back(tracks[i].StartDirection().Z());
      
        std::cout<<"run/subrun/evt:"<<run<<"/"<<subrun<<"/"<<event<<std::endl;	
        std::cout<<"beamPosx/beamPosy/beamPosz:"<<tracks[i].End().X()<<"/"<<tracks[i].End().Y()<<"/"<<tracks[i].End().Z()<<std::endl;
        std::cout<<"beamDirx/beamDiry/beamDirz:"<<tracks[i].StartDirection().X()<<"/"<<tracks[i].StartDirection().Y()<<"/"<<tracks[i].StartDirection().Z()<<std::endl;
        std::cout<<"###############################################################"<<std::endl;
      }
      //Get reconstructed beam momentum info
      auto & beammom = beamevent->GetRecoBeamMomenta();
      std::cout<<"==============================================================="<<std::endl;
      std::cout<<"beam mom size:"<<beammom.size()<<std::endl;
      for (size_t i = 0; i<beammom.size(); ++i){
        beamMomentum.push_back(beammom[i]);
        std::cout<<"beam mom["<<i<<"]:"<<beammom[i]<<" [GeV]"<<std::endl;
      }
      std::cout<<"==============================================================="<<std::endl;
      
      std::cout<<"\n*******************************************************"<<std::endl;
      std::cout<<"Moving on to the PFParticle section..."<<std::endl;	
      
      protoana::ProtoDUNEPFParticleUtils pfpUtil;
      auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
      
      
      // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
      // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These	
      // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those 
      // PFParticles considered as primary
      
      bool beamTriggerEvent = dataUtil.IsBeamTrigger(evt);
      if(beamTriggerEvent){
        std::cout << "This data event has a beam trigger" << std::endl;
      }
      
      
      const std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
      if(beamParticles.size() == 0){
        std::cerr << "We found no beam particles for this event... moving on" << std::endl;
        return;
      }
      
      
      n_beamparticle.push_back(beamParticles.size());
      std::cout<<"we have "<<beamParticles.size()<<" beam particle(s)"<<std::endl;

      // Trying to get the associated tracks from the PFParticle
      art::FindManyP<recob::Track> tracksFromPFParticles(recoParticles, evt, "pandoraTrack");
        
   
      for(const recob::PFParticle* particle : beamParticles){

      
        const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
        const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);

        if(thisTrack) {
          std::cout << "Beam particle is track-like" << std::endl;
          primtrk_trktag.push_back(1);
        }

        if(thisShower) { 
          std::cout << "Beam particle is shower-like" << std::endl;
          primtrk_trktag.push_back(-1);
        }
       
        pdg_code.push_back(particle->PdgCode());
        n_daughter.push_back(particle->NumDaughters());
        isPrimary.push_back(particle->IsPrimary());
        pfp_self.push_back(particle->Self());
          
        std::cout << "pdg code:"     << particle->PdgCode()      << std::endl;
        std::cout << "IsPrimary:"    << particle->IsPrimary()    << std::endl;
        std::cout << "NumDaughters:" << particle->NumDaughters() << std::endl;
        std::cout << "Self:"         << particle->Self()         << std::endl;	
        std::cout << "Parent:"       << particle->Parent()       << std::endl;


        if ( particle->NumDaughters() > 0 ) {
          for ( int ii = 0; ii < particle->NumDaughters(); ++ii ) {
            std::cout << "Daughter[" << ii << "]:" << particle->Daughter(ii) << std::endl;
            pfp_daughter.push_back(particle->Daughter(ii));
          }
        }
        else {
          pfp_daughter.push_back(-99);
        }


        // Find the particle vertex. We need the tracker tag here because we need to do a bit of
        // additional work if the PFParticle is track-like to find the vertex. 
        const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);

        std::cout << "Primary Beam Track Vertex: " << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << std::endl;
        
        primtrk_startx.push_back(vtx.X());
        primtrk_starty.push_back(vtx.Y());
        primtrk_startz.push_back(vtx.Z());
        
        // Get track direction
        if (thisTrack) {
          auto trackdir = thisTrack->StartDirection();
          std::cout<<"run/subrun/event:"<<run<<"/"<<subrun<<"/"<<event<<std::endl;	
          std::cout<<"trkDirx/trkDiry/trkDirz:"<<trackdir.X()<<"/"<<trackdir.Y()<<"/"<<trackdir.Z()<<std::endl;
          primtrk_Dirx.push_back(trackdir.X());
          primtrk_Diry.push_back(trackdir.Y());
          primtrk_Dirz.push_back(trackdir.Z());
        
          primtrklen.push_back(thisTrack->Length()); 
          std::cout<<"trk length: "<<thisTrack->Length()<<" [cm]"<<std::endl;
          primtrkID.push_back(thisTrack->ID());
          std::cout<<"trk ID: "<<thisTrack->ID()<<""<<std::endl; //HY::Fix me::trk ID seems wrong 
        
          if (tracks.size()){
              cosine_beam_primtrk = tracks[0].StartDirection().X()*trackdir.X()
              + tracks[0].StartDirection().Y()*trackdir.Y()
              + tracks[0].StartDirection().Z()*trackdir.Z();
          }
        }
        
        // Now we can look for the interaction point of the particle if one exists, i.e where the particle
        // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
        // check this by making sure we get a sensible position
        const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
        
        
        // Let's get the daughter PFParticles... we can do this simply without the utility
        for(const int daughterID : particle->Daughters()){
          // Daughter ID is the element of the original recoParticle vector
          const recob::PFParticle *daughterParticle = &(recoParticles->at(daughterID));
          const recob::Track* daughterTrack = pfpUtil.GetPFParticleTrack(*daughterParticle,evt,fPFParticleTag,fTrackerTag);

          std::cout << "Daughter "    << daughterID << " has "       << daughterParticle->NumDaughters() << " daughters" << std::endl;
          std::cout << "\tTrack ID: " << daughterTrack->ID()         << std::endl;
          std::cout << "\tPID: "      << daughterParticle->PdgCode() << std::endl;
          
          const TVector3 daughterVertex          = pfpUtil.GetPFParticleVertex( *daughterParticle, evt, fPFParticleTag, fTrackerTag );
          const TVector3 daughterSecondaryVertex = pfpUtil.GetPFParticleSecondaryVertex( *daughterParticle, evt, fPFParticleTag, fTrackerTag );

          double X, Y, Z;
          X = daughterVertex.X();
          Y = daughterVertex.Y();
          Z = daughterVertex.Z();

          std::cout << "\tVertex: " << X << ", " << Y << ", " << Z << std::endl;

          X = daughterSecondaryVertex.X();
          Y = daughterSecondaryVertex.Y();
          Z = daughterSecondaryVertex.Z();

          std::cout << "\tSecondaryVertex: " << X << ", " << Y << ", " << Z << std::endl;

          pfp_daughter_PDGs.push_back( daughterParticle->PdgCode() );
        }
        
        // For actually studying the objects, it is easier to have the daughters in their track and shower forms.
        // We can use the utility to get a vector of track-like and a vector of shower-like daughters
        const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
        const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  

        std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;

        for( size_t iD = 0; iD < trackDaughters.size(); ++iD ){
          std::cout << "Daughter track: " << iD << " " << trackDaughters[iD]->ID() << std::endl;
          pfp_daughter_TrackIDs.push_back( trackDaughters[iD]->ID() );
        }
        for( size_t iD = 0; iD < showerDaughters.size(); ++iD ){
          std::cout << "Daughter shower: " << iD << " " << showerDaughters[iD]->ID() << std::endl;
          pfp_daughter_ShowerIDs.push_back( showerDaughters[iD]->ID() );
        }

       
      }
      std::cout<<"*******************************************************"<<std::endl;
    } 

    if (beamevent->GetBITrigger() == 1){ 
      //Get CKov status
      ckov0status = beamevent->GetCKov0Status();
      ckov1status = beamevent->GetCKov1Status();
    } 
  }

  fTree->Fill();
}

void pionana::PionAnalyzerAbsCex::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("beamPosx",&beamPosx);
  fTree->Branch("beamPosy",&beamPosy);
  fTree->Branch("beamPosz",&beamPosz);
  fTree->Branch("beamDirx",&beamDirx);
  fTree->Branch("beamDiry",&beamDiry);
  fTree->Branch("beamDirz",&beamDirz);
  fTree->Branch("beamMomentum",&beamMomentum);
  fTree->Branch("tof", &tof, "tof/D");
  fTree->Branch("ckov0status", &ckov0status, "ckov0status/S");
  fTree->Branch("ckov1status", &ckov1status, "ckov1status/S");

  fTree->Branch("cosine_beam_primtrk", &cosine_beam_primtrk, "cosine_beam_primtrk/D");
  fTree->Branch("primtrk_startx",&primtrk_startx);
  fTree->Branch("primtrk_starty",&primtrk_starty);
  fTree->Branch("primtrk_startz",&primtrk_startz);
  fTree->Branch("primtrk_Dirx",&primtrk_Dirx);
  fTree->Branch("primtrk_Diry",&primtrk_Diry);
  fTree->Branch("primtrk_Dirz",&primtrk_Dirz);
  fTree->Branch("primtrklen",&primtrklen);
  fTree->Branch("primtrkID",&primtrkID);
  fTree->Branch("pdg_code", &pdg_code);
  fTree->Branch("n_beamparticle", &n_beamparticle);
  fTree->Branch("n_daughter", &n_daughter);
  fTree->Branch("isPrimary", &isPrimary);
  fTree->Branch("pfp_self", &pfp_self);
  //fTree->Branch("pfp_parent", &pfp_parent);
  fTree->Branch("pfp_daughter", &pfp_daughter);
  fTree->Branch("pfp_daughter_PDGs", &pfp_daughter_PDGs);
  fTree->Branch("pfp_daughter_TrackIDs", &pfp_daughter_TrackIDs);
  fTree->Branch("pfp_daughter_ShowerIDs", &pfp_daughter_ShowerIDs);
  fTree->Branch("primtrk_trktag", &primtrk_trktag);

}

void pionana::PionAnalyzerAbsCex::endJob()
{

}

void pionana::PionAnalyzerAbsCex::reset()
{
  beamPosx.clear();
  beamPosy.clear();
  beamPosz.clear();
  beamDirx.clear();
  beamDiry.clear();
  beamDirz.clear();
  beamMomentum.clear();

  primtrk_startx.clear();
  primtrk_starty.clear();
  primtrk_startz.clear();  
  primtrk_Dirx.clear();
  primtrk_Diry.clear();
  primtrk_Dirz.clear();
  primtrklen.clear();
  primtrkID.clear();
  primtrk_trktag.clear();

  pdg_code.clear();
  n_beamparticle.clear();
  n_daughter.clear();

  isPrimary.clear();
  pfp_self.clear();
  //pfp_parent.clear();
  pfp_daughter.clear();
  pfp_daughter_PDGs.clear();
  pfp_daughter_TrackIDs.clear();
  pfp_daughter_ShowerIDs.clear();

  cosine_beam_primtrk=-99;

}

DEFINE_ART_MODULE(pionana::PionAnalyzerAbsCex)
