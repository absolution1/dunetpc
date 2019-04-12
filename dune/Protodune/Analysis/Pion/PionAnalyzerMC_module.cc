////////////////////////////////////////////////////////////////////////
// Class:       PionAnalyzerMC
// Plugin Type: analyzer (art v3_00_00)
// File:        PionAnalyzerMC_module.cc
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
#include "dune/Protodune/Analysis/ProtoDUNEBeamlineUtils.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "art/Framework/Services/Optional/TFileService.h"

// ROOT includes
#include "TTree.h"

namespace pionana {
  class PionAnalyzerMC;
}


class pionana::PionAnalyzerMC : public art::EDAnalyzer {
public:
  explicit PionAnalyzerMC(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PionAnalyzerMC(PionAnalyzerMC const&) = delete;
  PionAnalyzerMC(PionAnalyzerMC&&) = delete;
  PionAnalyzerMC& operator=(PionAnalyzerMC const&) = delete;
  PionAnalyzerMC& operator=(PionAnalyzerMC&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();

private:

  // Declare member data here.
  const art::InputTag fTrackModuleLabel;

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;

  int MC;

  int nPiPlus_truth, nPiMinus_truth, nPi0_truth;
  int nProton_truth, nNeutron_truth;
  int MC_true_PDG;
  int geantGood_PDG;
  std::string geantGood_EndProcess;
  std::string MC_true_EndProcess;
  std::string MC_true_Process;

  int MC_origin;
  std::vector< bool > MC_daughter_good_reco;
  std::vector< int > MC_daughter_true_PDGs;
  std::vector< int > MC_daughter_true_IDs;

  std::vector< bool > MC_shower_daughter_good_reco;
  std::vector< int > MC_shower_daughter_true_PDGs;
  std::vector< int > MC_shower_daughter_true_IDs;

  std::vector< int > geantGood_daughter_PDGs;
  std::vector< int > geantGood_daughter_IDs;
  std::vector< double > geantGood_daughter_lens;
  std::vector< double > MC_daughter_true_lens;
  std::vector< double > MC_shower_daughter_true_lens;
  std::vector< double > MC_daughter_reco_lens;
  std::vector< bool > MC_daughter_shower_good_reco;
  bool MC_good_reco;


  double startX, startY, startZ;
  double endX, endY, endZ;
  double len;
  int broken_candidate;
  double trackDirX, trackDirY, trackDirZ;
  int beamTrackID;
  std::vector< int > daughter_trackID;
  std::vector< int > daughter_showerID;

  std::vector< double > dEdX, dQdX, resRange;
  std::vector< std::vector< double > > daughter_dEdX, daughter_dQdX, daughter_resRange;
  double vtxX, vtxY, vtxZ; 
  std::vector< double > daughter_startX, daughter_endX;
  std::vector< double > daughter_startY, daughter_endY;
  std::vector< double > daughter_startZ, daughter_endZ;
  std::vector< double > daughter_shower_startX;
  std::vector< double > daughter_shower_startY;
  std::vector< double > daughter_shower_startZ;
  std::vector< double > daughter_len;

  int nTrackDaughters, nShowerDaughters;

  int type;
  int nBeamParticles;

  std::string fCalorimetryTag;
  
  std::string fTrackerTag;    
  std::string fShowerTag;     
  std::string fPFParticleTag; 
  std::string fGeneratorTag;
  bool fVerbose;             
  fhicl::ParameterSet dataUtil;
};


pionana::PionAnalyzerMC::PionAnalyzerMC(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),

  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils"))


{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pionana::PionAnalyzerMC::analyze(art::Event const& evt)
{
  //reset containers
  reset();  


  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

   
  // Get the PFParticle utility
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  nBeamParticles = beamParticles.size();

  if(beamParticles.size() == 0){
    std::cerr << "We found no beam particles for this event... moving on" << std::endl;
    return;
  }

  std::cout << "Found " << nBeamParticles << " beamParticles" << std::endl;

  // We can now look at these particles
  const recob::PFParticle* particle = beamParticles.at(0);
    const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
    if(thisTrack != 0x0){
      std::cout << "Beam particle is track-like " << thisTrack->ID() << std::endl;
      type = 13;
    }
    if(thisShower != 0x0){
      std::cout << "Beam particle is shower-like" << std::endl;
      type = 11;
    }



    // Find the particle vertex. We need the tracker tag here because we need to do a bit of
    // additional work if the PFParticle is track-like to find the vertex. 
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    std::cout << "Vertex: " << vtx[0] << " " << vtx[1] << " " << vtx[2] << std::endl;

    //Get the TPC ID for the beginning of the track
    if( thisTrack ){


      startX = thisTrack->Trajectory().Start().X();
      startY = thisTrack->Trajectory().Start().Y();
      startZ = thisTrack->Trajectory().Start().Z();
      endX = thisTrack->Trajectory().End().X();
      endY = thisTrack->Trajectory().End().Y();
      endZ = thisTrack->Trajectory().End().Z();

      auto startDir = thisTrack->StartDirection();
      auto endDir   = thisTrack->EndDirection();

      //try flipping
      if( startZ > endZ ){
        std::cout << "startZ > endZ: " << startZ << " " << endZ << std::endl;
        endX = thisTrack->Trajectory().Start().X();
        endY = thisTrack->Trajectory().Start().Y();
        endZ = thisTrack->Trajectory().Start().Z();
        startX = thisTrack->Trajectory().End().X();
        startY = thisTrack->Trajectory().End().Y();
        startZ = thisTrack->Trajectory().End().Z();
        
        trackDirX =  -1. * endDir.X(); 
        trackDirY =  -1. * endDir.Y(); 
        trackDirZ =  -1. * endDir.Z(); 

      }
      else{
        std::cout << "endZ > startZ: " << startZ << " " << endZ << std::endl;
        trackDirX =  startDir.X(); 
        trackDirY =  startDir.Y(); 
        trackDirZ =  startDir.Z(); 
      }
      std::cout << "StartDir: " << startDir.Z() << std::endl;
      std::cout << "EndDir: "   << endDir.Z() << std::endl;

      std::cout << "trackDirX: " << trackDirX << std::endl;
      std::cout << "trackDirY: " << trackDirY << std::endl;
      std::cout << "trackDirZ: " << trackDirZ << std::endl;

      len  = thisTrack->Length();    
      beamTrackID = thisTrack->ID();
      
      std::cout << "Start: " << startX << " " << startY << " " << startZ << std::endl;
      std::cout << "End: " << endX << " " << endY << " " << endZ << std::endl;
      std::cout << "len: " << len << std::endl;

      // Now we can look for the interaction point of the particle if one exists, i.e where the particle
      // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
      // check this by making sure we get a sensible position
      const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
      vtxX = interactionVtx.X();
      vtxY = interactionVtx.Y();
      vtxZ = interactionVtx.Z();
      std::cout << "Interaction Vertex: " << vtxX << " " <<  vtxY << " " <<  vtxZ << std::endl;
      

      // Let's get the daughter PFParticles... we can do this simply without the utility
      for(const int daughterID : particle->Daughters()){
        // Daughter ID is the element of the original recoParticle vector
        const recob::PFParticle *daughterParticle = &(recoParticles->at(daughterID));
        std::cout << "Daughter " << daughterID << " PDG: " << daughterParticle->PdgCode() << std::endl; 
      }
 
      // For actually studying the objects, it is easier to have the daughters in their track and shower forms.
      // We can use the utility to get a vector of track-like and a vector of shower-like daughters
      const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
      const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  
      std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;
      std::cout << std::endl;

      for( size_t i = 0; i < trackDaughters.size(); ++i ){
        std::cout << "Track daughter " << i << " has len " << trackDaughters[i]->Length() << std::endl; 
        daughter_len.push_back( trackDaughters[i]->Length() );
      }

      for( size_t i = 0; i < showerDaughters.size(); ++i ){
        std::cout << "Shower daughter " << i << " Starts at " << showerDaughters[i]->ShowerStart().X() << " " << showerDaughters[i]->ShowerStart().Y() << " " << showerDaughters[i]->ShowerStart().Z() << std::endl;
        daughter_showerID.push_back( showerDaughters[i]->ID() );
        daughter_shower_startX.push_back( showerDaughters[i]->ShowerStart().X() );
        daughter_shower_startY.push_back( showerDaughters[i]->ShowerStart().Y() );
        daughter_shower_startZ.push_back( showerDaughters[i]->ShowerStart().Z() );
      }

      nTrackDaughters = trackDaughters.size();
      nShowerDaughters = showerDaughters.size();



      art::ServiceHandle<cheat::BackTrackerService> bt_serv;
      art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;

      MC = 1;
      // Get the truth utility to help us out
      protoana::ProtoDUNETruthUtils truthUtil;
      // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
      // simulation has fGeneratorTag = "generator"
      auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
      std::cout << "MCTruth origin: ";
      switch( (*mcTruths)[0].Origin() ){
        case simb::kCosmicRay: 
          std::cout << "Cosmic" << std::endl;
          break;
        case simb::kSingleParticle:
          std::cout << "Beam" << std::endl;
          break;
        default:
          std::cout << "Other" << std::endl;
      }

      // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
      // of these, so we pass the first element into the function to get the good particle
      const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
      if(geantGoodParticle){
        geantGood_PDG = geantGoodParticle->PdgCode();
        std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() << std::endl;
        std::cout << "ID: " << geantGoodParticle->TrackId() << std::endl;
        std::cout << "Mother: " << geantGoodParticle->Mother() << std::endl;
      }
      std::cout << "Has " << geantGoodParticle->NumberDaughters() << " daughters" << std::endl;


      //This gets the MCParticle contributing the most to the track
      const simb::MCParticle* trueParticle = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackerTag);
      if( trueParticle ){ 
        std::cout << "True particle from beam track: " << trueParticle->PdgCode() << std::endl;
        MC_true_PDG = trueParticle->PdgCode();
        MC_true_EndProcess = trueParticle->EndProcess();
        MC_true_Process = trueParticle->Process();

        MC_origin = pi_serv->TrackIdToMCTruth_P(trueParticle->TrackId())->Origin();
        std::cout << "Track-to-True ID: " << trueParticle->TrackId() << std::endl;

        //What created the MCPraticle that created the track?
        std::cout << "Track-to-True origin: ";
        switch( MC_origin ){
          case simb::kCosmicRay: 
            std::cout << "Cosmic" << std::endl;
            break;
          case simb::kSingleParticle:
            std::cout << "Beam" << std::endl;
            break;
          default:
            std::cout << "Other" << std::endl;
        }

        //Reconstructed beam track corresponds to the true beam particle
        if( trueParticle->TrackId() == geantGoodParticle->TrackId() ){
          MC_good_reco = true;
        }

        std::cout << "MCParticle trajectory:" << std::endl;
        //Get the beginning of the MCParticle and project this to z=0 as if it was the beamline
        std::cout << trueParticle->Trajectory().X(0) << " " << trueParticle->Trajectory().Y(0) << " " << trueParticle->Trajectory().Z(0) << std::endl;
        size_t nT = trueParticle->Trajectory().size();
        if( nT > 0 ) std::cout << trueParticle->Trajectory().X(nT-1) << " " << trueParticle->Trajectory().Y(nT-1) << " " << trueParticle->Trajectory().Z(nT-1) << std::endl;
        
      }
      else std::cout << "Couldn't get trueParticle" << std::endl;
       

      const sim::ParticleList & plist = pi_serv->ParticleList(); 


      for( int i = 0; i < geantGoodParticle->NumberDaughters(); ++i ){
        int daughterID = geantGoodParticle->Daughter(i);

        //Skip photons, neutrons, the nucleus
      //  if( plist[ daughterID ]->PdgCode() == 22 || plist[ daughterID ]->PdgCode() == 2112 || plist[ daughterID ]->PdgCode() > 1000000000  ) continue;

        std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
        auto part = plist[ daughterID ];
        int pid = part->PdgCode();
        geantGood_daughter_PDGs.push_back(pid);
        geantGood_daughter_IDs.push_back( part->TrackId() );      
        geantGood_daughter_lens.push_back( part->Trajectory().TotalLength() );
        std::cout << "PID: " << pid << std::endl;
        std::cout << "Start: " << part->Position(0).X() << " " << part->Position(0).Y() << " " << part->Position(0).Z() << std::endl;
        std::cout << "End: " << part->EndPosition().X() << " " << part->EndPosition().Y() << " " << part->EndPosition().Z() << std::endl;
        std::cout << "Len: " << part->Trajectory().TotalLength() << std::endl;

        if( pid == 211 )  nPiPlus_truth++;
        if( pid == -211 ) nPiMinus_truth++;
        if( pid == 111 )  nPi0_truth++;
        if( pid == 2212 ) nProton_truth++;

      }

      //Add in check between true good particles and the reconstructed tracks
      for( size_t i = 0; i < trackDaughters.size(); ++i ){

        auto daughterTrackFromRecoTrack = trackDaughters[i];

        bool found_daughter = false;
        const simb::MCParticle* daughterParticleFromRecoTrack = truthUtil.GetMCParticleFromRecoTrack(*daughterTrackFromRecoTrack, evt, fTrackerTag); 

        if( daughterParticleFromRecoTrack ){
          int loc = 0;
          for( size_t j = 0; j < geantGood_daughter_IDs.size(); ++j ){
            if ( geantGood_daughter_IDs[j] == daughterParticleFromRecoTrack->TrackId() ){
              found_daughter = true;
              loc = j;
              break;
            }
          }

          MC_daughter_good_reco.push_back( found_daughter );
          if( found_daughter ){
            MC_daughter_true_PDGs.push_back( geantGood_daughter_PDGs[loc] ); 
            MC_daughter_true_IDs.push_back( geantGood_daughter_IDs[loc] );
            MC_daughter_true_lens.push_back( geantGood_daughter_lens[loc] );
          }
          else{

            MC_daughter_true_PDGs.push_back( daughterParticleFromRecoTrack->PdgCode() );
            MC_daughter_true_IDs.push_back( daughterParticleFromRecoTrack->TrackId() );
            MC_daughter_true_lens.push_back( daughterParticleFromRecoTrack->Trajectory().TotalLength() );
          }
        }

      }

      for( size_t i = 0; i < showerDaughters.size(); ++i ){
        auto daughterShowerFromRecoTrack = showerDaughters[i];

        bool found_daughter = false;
        const simb::MCParticle* daughterParticleFromRecoShower = truthUtil.GetMCParticleFromRecoShower(*daughterShowerFromRecoTrack, evt, fShowerTag); 
        std::cout << "Shower daughter " << i << " " << daughterParticleFromRecoShower << std::endl;

        if( daughterParticleFromRecoShower ){
          std::cout << "Shower daughter ID: " << daughterParticleFromRecoShower->TrackId() << std::endl;
          int loc = 0;
          for( size_t j = 0; j < geantGood_daughter_IDs.size(); ++j ){
            if ( geantGood_daughter_IDs[j] == daughterParticleFromRecoShower->TrackId() ){
              found_daughter = true;
              loc = j;
              break;
            }
          }

          MC_shower_daughter_good_reco.push_back( found_daughter );
          if( found_daughter ){
            MC_shower_daughter_true_PDGs.push_back( geantGood_daughter_PDGs[loc] ); 
            MC_shower_daughter_true_IDs.push_back( geantGood_daughter_IDs[loc] );
            MC_shower_daughter_true_lens.push_back( geantGood_daughter_lens[loc] );
          }
          else{
            MC_shower_daughter_true_PDGs.push_back( daughterParticleFromRecoShower->PdgCode() );
            MC_shower_daughter_true_IDs.push_back( daughterParticleFromRecoShower->TrackId() );
            MC_shower_daughter_true_lens.push_back( daughterParticleFromRecoShower->Trajectory().TotalLength() );
          }
        }

      }


      geantGood_EndProcess = geantGoodParticle->EndProcess();



      //Want to see all of the truth particles that contributed to this track
      std::vector< std::pair< const simb::MCParticle*, double > > contribParts = truthUtil.GetAllMCParticlesFromRecoTrack(*thisTrack, evt, fTrackerTag);
      std::cout << contribParts.size() << " Truth Particles Contributed to this track" << std::endl;
      for( size_t ip = 0; ip < contribParts.size(); ++ip ){
        auto part = contribParts.at( ip ).first;
        double energy = contribParts.at( ip ).second;
        std::cout << ip << " " << part->TrackId() << " " << part->PdgCode() << " " << energy << std::endl;
      }


      //Calorimetry 
      //
      protoana::ProtoDUNETrackUtils trackUtil;
      std::vector< anab::Calorimetry> calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
      std::cout << "Planes: " << calo[0].PlaneID().toString() << " " << calo[1].PlaneID().toString()  << " " << calo[2].PlaneID().toString() << std::endl;
      auto calo_dQdX = calo[0].dQdx();
      auto calo_dEdX = calo[0].dEdx();
      auto calo_range = calo[0].ResidualRange();
      for( size_t i = 0; i < calo_dQdX.size(); ++i ){

        std::cout << calo_dQdX[i] << " " << calo_dEdX[i] << " " << calo_range[i] << std::endl;
        dQdX.push_back( calo_dQdX[i] );
        dEdX.push_back( calo_dEdX[i] );
        resRange.push_back( calo_range[i] );
      }


      //Go through the track-like daughters and save their calorimetry
      for( size_t i = 0; i < trackDaughters.size(); ++i ){
        auto daughterTrack = trackDaughters.at(i);
        
        daughter_startX.push_back( daughterTrack->Trajectory().Start().X() );
        daughter_startY.push_back( daughterTrack->Trajectory().Start().Y() );
        daughter_startZ.push_back( daughterTrack->Trajectory().Start().Z() );
        daughter_endX.push_back( daughterTrack->Trajectory().End().X() );
        daughter_endY.push_back( daughterTrack->Trajectory().End().Y() );
        daughter_endZ.push_back( daughterTrack->Trajectory().End().Z() );

        daughter_trackID.push_back( daughterTrack->ID() );


        std::vector< anab::Calorimetry > dummy_calo = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
        auto dummy_dQdx = dummy_calo[0].dQdx();
        auto dummy_dEdx = dummy_calo[0].dEdx();
        auto dummy_Range = dummy_calo[0].ResidualRange();
 
        daughter_dQdX.push_back( std::vector<double>() );   
        daughter_resRange.push_back( std::vector<double>() );
        daughter_dEdX.push_back( std::vector<double>() );

        for( size_t i = 0; i < dummy_dQdx.size(); ++i ){
          daughter_dQdX.back().push_back( dummy_dQdx[i] );
          daughter_resRange.back().push_back( dummy_Range[i] );
          daughter_dEdX.back().push_back( dummy_dEdx[i] );
        }
      }

    } 


  fTree->Fill();
}

void pionana::PionAnalyzerMC::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");
  fTree->Branch("startX", &startX);
  fTree->Branch("startY", &startY);
  fTree->Branch("startZ", &startZ);
  fTree->Branch("endX", &endX);
  fTree->Branch("endY", &endY);
  fTree->Branch("endZ", &endZ);
  fTree->Branch("len", &len);
  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("type", &type);
  fTree->Branch("nBeamParticles", &nBeamParticles);

  fTree->Branch("trackDirX", &trackDirX);
  fTree->Branch("trackDirY", &trackDirY);
  fTree->Branch("trackDirZ", &trackDirZ);
  fTree->Branch("beamTrackID", &beamTrackID);
  fTree->Branch("daughter_trackID", &daughter_trackID);
  fTree->Branch("daughter_showerID", &daughter_showerID);

  fTree->Branch("MC", &MC);
  fTree->Branch("dQdX", &dQdX);
  fTree->Branch("dEdX", &dEdX);
  fTree->Branch("resRange", &resRange);
  fTree->Branch("daughter_dQdX", &daughter_dQdX);
  fTree->Branch("daughter_dEdX", &daughter_dEdX);
  fTree->Branch("daughter_resRange", &daughter_resRange);
  fTree->Branch("daughter_len", &daughter_len);
  fTree->Branch("daughter_startX", &daughter_startX);
  fTree->Branch("daughter_startY", &daughter_startY);
  fTree->Branch("daughter_startZ", &daughter_startZ);
  fTree->Branch("vtxX", &vtxX);
  fTree->Branch("vtxY", &vtxY);
  fTree->Branch("vtxZ", &vtxZ);
  fTree->Branch("daughter_endX", &daughter_endX);
  fTree->Branch("daughter_endY", &daughter_endY);
  fTree->Branch("daughter_endZ", &daughter_endZ);
  fTree->Branch("daughter_shower_startX", &daughter_shower_startX);
  fTree->Branch("daughter_shower_startY", &daughter_shower_startY);
  fTree->Branch("daughter_shower_startZ", &daughter_shower_startZ);
  fTree->Branch("nTrackDaughters", &nTrackDaughters);
  fTree->Branch("nShowerDaughters", &nShowerDaughters);
  fTree->Branch("nProton_truth", &nProton_truth);
  fTree->Branch("nPi0_truth", &nPi0_truth);
  fTree->Branch("nPiPlus_truth", &nPiPlus_truth);
  fTree->Branch("nPiMinus_truth", &nPiMinus_truth);
  fTree->Branch("MC_true_PDG", &MC_true_PDG);
  fTree->Branch("geantGood_PDG", &geantGood_PDG);
  fTree->Branch("geantGood_EndProcess", &geantGood_EndProcess);
  fTree->Branch("MC_true_EndProcess", &MC_true_EndProcess);
  fTree->Branch("MC_true_Process", &MC_true_Process);
  fTree->Branch("MC_origin", &MC_origin);
  fTree->Branch("MC_good_reco", &MC_good_reco);
  fTree->Branch("MC_daughter_good_reco", &MC_daughter_good_reco);
  fTree->Branch("MC_daughter_true_PDGs", &MC_daughter_true_PDGs);
  fTree->Branch("MC_daughter_true_IDs", &MC_daughter_true_IDs);

  fTree->Branch("MC_shower_daughter_good_reco", &MC_shower_daughter_good_reco);
  fTree->Branch("MC_shower_daughter_true_PDGs", &MC_shower_daughter_true_PDGs);
  fTree->Branch("MC_shower_daughter_true_IDs", &MC_shower_daughter_true_IDs);

  fTree->Branch("geantGood_daughter_PDGs", &geantGood_daughter_PDGs);
  fTree->Branch("geantGood_daughter_IDs", &geantGood_daughter_IDs);
  fTree->Branch("geantGood_daughter_lens", &geantGood_daughter_lens);
  fTree->Branch("MC_daughter_true_lens", &MC_daughter_true_lens);
  fTree->Branch("MC_shower_daughter_true_lens", &MC_shower_daughter_true_lens);
  fTree->Branch("MC_daughter_reco_lens", &MC_daughter_reco_lens);
  fTree->Branch("MC_daughter_shower_good_reco", &MC_daughter_shower_good_reco);


}

void pionana::PionAnalyzerMC::endJob()
{

}

void pionana::PionAnalyzerMC::reset()
{
  startX = -1;
  startY = -1;
  startZ = -1;
  endX = -1;
  endY = -1;
  endZ = -1;

  len = -1;
  broken_candidate = 0;
  type = -1;
  nBeamParticles = 0;

  MC = 0;
  nProton_truth = 0;
  nPi0_truth = 0;
  nPiPlus_truth = 0;
  nPiMinus_truth = 0;
  MC_true_PDG = 0;
  geantGood_PDG = 0;
  geantGood_EndProcess ="";
  MC_true_EndProcess ="";
  MC_true_Process ="";
  MC_origin = -1;

  MC_good_reco = false;
  MC_daughter_good_reco.clear();
  MC_daughter_true_PDGs.clear();
  MC_daughter_true_IDs.clear();

  MC_shower_daughter_good_reco.clear();
  MC_shower_daughter_true_PDGs.clear();
  MC_shower_daughter_true_IDs.clear();

  geantGood_daughter_PDGs.clear();
  geantGood_daughter_lens.clear();
  MC_daughter_true_lens.clear();
  MC_shower_daughter_true_lens.clear();
  MC_daughter_reco_lens.clear();
  geantGood_daughter_IDs.clear();
  MC_daughter_shower_good_reco.clear();

  nTrackDaughters = -1;
  nShowerDaughters = -1;

  dQdX.clear();
  dEdX.clear();
  vtxX = -1.;
  vtxY = -1.;
  vtxZ = -1.;
  daughter_startX.clear();
  daughter_startY.clear();
  daughter_startZ.clear();
  daughter_endX.clear();
  daughter_endY.clear();
  daughter_endZ.clear();

  daughter_shower_startX.clear();
  daughter_shower_startY.clear();
  daughter_shower_startZ.clear();

  resRange.clear();

  daughter_dQdX.clear();
  daughter_dEdX.clear();
  daughter_resRange.clear();
  daughter_len.clear();

  beamTrackID = -1;
  daughter_trackID.clear();
  daughter_showerID.clear();

}

DEFINE_ART_MODULE(pionana::PionAnalyzerMC)
