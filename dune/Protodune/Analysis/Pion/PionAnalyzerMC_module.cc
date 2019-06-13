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

#include "art_root_io/TFileService.h"

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




  /******************************/
  //Truth level info of the primary beam particle
  //that generated the event
  int true_beam_PDG;
  std::string true_beam_EndProcess;
  double true_beam_EndVertex_X;
  double true_beam_EndVertex_Y;
  double true_beam_EndVertex_Z;
  double true_beam_Start_X;
  double true_beam_Start_Y;
  double true_beam_Start_Z;

    //Truth level info of the daughter MCParticles coming out of the 
    //true primary particle
  std::vector< int > true_beam_daughter_PDGs;
  std::vector< int > true_beam_daughter_IDs;
  std::vector< double > true_beam_daughter_lens;

  //How many of each true particle came out of the true primary beam particle?
  int nPiPlus_truth, nPiMinus_truth, nPi0_truth;
  int nProton_truth, nNeutron_truth;
  /*****************************/
  



  //Reconstructed track info
  double startX, startY, startZ;
  double endX, endY, endZ;
  double vtxX, vtxY, vtxZ; 
  double len;
  double trackDirX, trackDirY, trackDirZ;
  std::vector< double > dEdX, dQdX, resRange;
  int beamTrackID;

  std::string reco_beam_truth_EndProcess; //What process ended the reco beam particle
  std::string reco_beam_truth_Process;    //What process created the reco beam particle
  
  int reco_beam_truth_PDG; //What is the PDG of the true MC particle contributing the most to 
                           //the reconstructed beam track

  
  bool reco_beam_good; //Does the true particle contributing most to the 
                       //reconstructed beam track coincide with the actual
                       //beam particle that generated the event
                       

  int reco_beam_truth_origin; //What is the origin of the reconstructed beam track?



  //Truth-level info of the reconstructed particles coming out of the 
  //reconstructed beam track
  std::vector< bool >   reco_beam_truth_daughter_good_reco;
  std::vector< int >    reco_beam_truth_daughter_true_PDGs;
  std::vector< int >    reco_beam_truth_daughter_true_IDs;
  std::vector< double > reco_beam_truth_daughter_true_lens;

  std::vector< double > reco_beam_truth_daughter_shower_true_lens;
  std::vector< bool >   reco_beam_truth_daughter_shower_good_reco;
  std::vector< int >    reco_beam_truth_daughter_shower_true_PDGs;
  std::vector< int >    reco_beam_truth_daughter_shower_true_IDs;


  //Reco-level info of the reconstructed daughters coming out of the
  //reconstructed beam tracl
  std::vector< int > reco_daughter_trackID;
  std::vector< int > reco_daughter_showerID;
  std::vector< std::vector< double > > reco_daughter_dEdX, reco_daughter_dQdX, reco_daughter_resRange;
  std::vector< double > reco_daughter_startX, reco_daughter_endX;
  std::vector< double > reco_daughter_startY, reco_daughter_endY;
  std::vector< double > reco_daughter_startZ, reco_daughter_endZ;
  std::vector< double > reco_daughter_shower_startX;
  std::vector< double > reco_daughter_shower_startY;
  std::vector< double > reco_daughter_shower_startZ;
  std::vector< double > reco_daughter_len;

  int nTrackDaughters, nShowerDaughters;

  int type;
  int nBeamParticles;



  //FCL pars
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

   
  // Get various utilities 
  protoana::ProtoDUNEPFParticleUtils                    pfpUtil;
  protoana::ProtoDUNETrackUtils                         trackUtil;
  art::ServiceHandle<cheat::BackTrackerService>         bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  protoana::ProtoDUNETruthUtils                         truthUtil;
  ////////////////////////////////////////
  


  // This gets the true beam particle that generated the event
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  const simb::MCParticle* true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
  ////////////////////////////


  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  nBeamParticles = beamParticles.size();

  if(beamParticles.size() == 0){
    std::cerr << "We found no beam particles for this event... moving on" << std::endl;
    return;
  }
  //////////////////////////////////////////////////////////////////


  // Get the reconstructed PFParticle tagged as beam by Pandora
  const recob::PFParticle* particle = beamParticles.at(0);


  //Reco track vertex and direction/length
  const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
  std::cout << "PFParticle Vertex: " << vtx[0] << " " << vtx[1] << " " << vtx[2] << std::endl;



  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
  const simb::MCParticle* trueParticle = 0x0;
  if(thisShower != 0x0){
    std::cout << "Beam particle is shower-like" << std::endl;
    type = 11;
    //This gets the MCParticle contributing the most to the shower
    trueParticle = truthUtil.GetMCParticleFromRecoShower(*thisShower, evt, fShowerTag);
  }
  if(thisTrack != 0x0){

    //This gets the MCParticle contributing the most to the track
    trueParticle = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackerTag);


    std::cout << "Beam particle is track-like " << thisTrack->ID() << std::endl;
    type = 13;

    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    vtxX = interactionVtx.X();
    vtxY = interactionVtx.Y();
    vtxZ = interactionVtx.Z();
    std::cout << "Secondary Vertex: " << vtxX << " " <<  vtxY << " " <<  vtxZ << std::endl;
    ////////////////////////////////////////////
    
  }
  //////////////////////////////////////////////////

  if( trueParticle ){

    //Check that this is the correct true particle
    if( trueParticle->TrackId() == true_beam_particle->TrackId() ){
      reco_beam_good = true;
    }

  }
  //////////////////////////////////////////////////////////////

  reco_beam_truth_PDG = trueParticle->PdgCode();

  reco_beam_truth_Process = trueParticle->Process();
  reco_beam_truth_origin = pi_serv->TrackIdToMCTruth_P(trueParticle->TrackId())->Origin();
  //What created the MCPraticle that created the track?
  std::cout << "Track-to-True origin: ";
  switch( reco_beam_truth_origin ){
    case simb::kCosmicRay: 
      std::cout << "Cosmic" << std::endl;
      break;
    case simb::kSingleParticle:
      std::cout << "Beam" << std::endl;
      break;
    default:
      std::cout << "Other" << std::endl;
  }




  //Some truth information
  //
  //
  //What created the true particle
  true_beam_EndProcess = true_beam_particle->EndProcess();

  ///////////////////////////////////////
    
  
  true_beam_PDG         = true_beam_particle->PdgCode();
  true_beam_EndVertex_X = true_beam_particle->EndX();
  true_beam_EndVertex_Y = true_beam_particle->EndY();
  true_beam_EndVertex_Z = true_beam_particle->EndZ();
  true_beam_Start_X     = true_beam_particle->Position(0).X();
  true_beam_Start_Y     = true_beam_particle->Position(0).Y();
  true_beam_Start_Z     = true_beam_particle->Position(0).Z();


  //Look at true daughters coming out of the true beam particle
  std::cout << "Has " << true_beam_particle->NumberDaughters() << " daughters" << std::endl;

  const sim::ParticleList & plist = pi_serv->ParticleList(); 
  for( int i = 0; i < true_beam_particle->NumberDaughters(); ++i ){
    int daughterID = true_beam_particle->Daughter(i);

    std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
    auto part = plist[ daughterID ];
    int pid = part->PdgCode();
    true_beam_daughter_PDGs.push_back(pid);
    true_beam_daughter_IDs.push_back( part->TrackId() );      
    true_beam_daughter_lens.push_back( part->Trajectory().TotalLength() );
    std::cout << "PID: " << pid << std::endl;
    std::cout << "Start: " << part->Position(0).X() << " " << part->Position(0).Y() << " " << part->Position(0).Z() << std::endl;
    std::cout << "End: " << part->EndPosition().X() << " " << part->EndPosition().Y() << " " << part->EndPosition().Z() << std::endl;
    std::cout << "Len: " << part->Trajectory().TotalLength() << std::endl;

    if( pid == 211 )  ++nPiPlus_truth;
    if( pid == -211 ) ++nPiMinus_truth;
    if( pid == 111 )  ++nPi0_truth;
    if( pid == 2212 ) ++nProton_truth;
    if( pid == 2112 ) ++nNeutron_truth;

  }


  if( thisTrack ){

/*
    if( trueParticle ){

      //Check that this is the correct true particle
      if( trueParticle->TrackId() == true_beam_particle->TrackId() ){
        reco_beam_good = true;
      }

    }
    //////////////////////////////////////////////////////////////
    
    
    beamTrackID = thisTrack->ID();
    

    reco_beam_truth_Process = trueParticle->Process();
    reco_beam_truth_origin = pi_serv->TrackIdToMCTruth_P(trueParticle->TrackId())->Origin();
    //What created the MCPraticle that created the track?
    std::cout << "Track-to-True origin: ";
    switch( reco_beam_truth_origin ){
      case simb::kCosmicRay: 
        std::cout << "Cosmic" << std::endl;
        break;
      case simb::kSingleParticle:
        std::cout << "Beam" << std::endl;
        break;
      default:
        std::cout << "Other" << std::endl;
    }

    reco_beam_truth_PDG = trueParticle->PdgCode();
    */

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
    std::cout << "Start: " << startX << " " << startY << " " << startZ << std::endl;
    std::cout << "End: " << endX << " " << endY << " " << endZ << std::endl;
    std::cout << "len: " << len << std::endl;
    ////////////////////////////////////////////////////////////////




    //Primary Track Calorimetry 
    std::vector< anab::Calorimetry> calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
    std::cout << "Planes: " << calo[0].PlaneID().toString() << " " << calo[1].PlaneID().toString()  << " " << calo[2].PlaneID().toString() << std::endl;
    auto calo_dQdX = calo[0].dQdx();
    auto calo_dEdX = calo[0].dEdx();
    auto calo_range = calo[0].ResidualRange();
    for( size_t i = 0; i < calo_dQdX.size(); ++i ){
      dQdX.push_back( calo_dQdX[i] );
      dEdX.push_back( calo_dEdX[i] );
      resRange.push_back( calo_range[i] );
    }
    ////////////////////////////////////////////



    //Looking at reco daughters from the reco beam track
    const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
    const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  
    std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;
    std::cout << std::endl;


    nTrackDaughters = trackDaughters.size();

    for( size_t i = 0; i < trackDaughters.size(); ++i ){
      std::cout << "Track daughter " << i << " has len " << trackDaughters[i]->Length() << std::endl; 
      auto daughterTrack = trackDaughters.at(i);
      
      reco_daughter_len.push_back( daughterTrack->Length() );
      reco_daughter_startX.push_back( daughterTrack->Trajectory().Start().X() );
      reco_daughter_startY.push_back( daughterTrack->Trajectory().Start().Y() );
      reco_daughter_startZ.push_back( daughterTrack->Trajectory().Start().Z() );
      reco_daughter_endX.push_back( daughterTrack->Trajectory().End().X() );
      reco_daughter_endY.push_back( daughterTrack->Trajectory().End().Y() );
      reco_daughter_endZ.push_back( daughterTrack->Trajectory().End().Z() );
      reco_daughter_trackID.push_back( daughterTrack->ID() );


      std::vector< anab::Calorimetry > dummy_calo = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
      auto dummy_dQdx = dummy_calo[0].dQdx();
      auto dummy_dEdx = dummy_calo[0].dEdx();
      auto dummy_Range = dummy_calo[0].ResidualRange();
 
      reco_daughter_dQdX.push_back( std::vector<double>() );   
      reco_daughter_resRange.push_back( std::vector<double>() );
      reco_daughter_dEdX.push_back( std::vector<double>() );

      for( size_t j = 0; j < dummy_dQdx.size(); ++j ){
        reco_daughter_dQdX.back().push_back( dummy_dQdx[j] );
        reco_daughter_resRange.back().push_back( dummy_Range[j] );
        reco_daughter_dEdX.back().push_back( dummy_dEdx[j] );
      }


      ///Try to match the reconstructed daughter tracks from the reco'd/tagged beam track
      //   to the true daughter particles coming out of the true beam track
      bool found_daughter = false;
      const simb::MCParticle* daughterParticleFromRecoTrack = truthUtil.GetMCParticleFromRecoTrack(*daughterTrack, evt, fTrackerTag); 

      if( daughterParticleFromRecoTrack ){
        int loc = 0;
        for( size_t j = 0; j < true_beam_daughter_IDs.size(); ++j ){
          
          //do the checking
          if ( true_beam_daughter_IDs[j] == daughterParticleFromRecoTrack->TrackId() ){
            found_daughter = true;
            loc = j;
            break;
          }
        }

        reco_beam_truth_daughter_good_reco.push_back( found_daughter );
        if( found_daughter ){
          reco_beam_truth_daughter_true_PDGs.push_back( true_beam_daughter_PDGs[loc] ); 
          reco_beam_truth_daughter_true_IDs.push_back( true_beam_daughter_IDs[loc] );
          reco_beam_truth_daughter_true_lens.push_back( true_beam_daughter_lens[loc] );
        }
        else{

          reco_beam_truth_daughter_true_PDGs.push_back( daughterParticleFromRecoTrack->PdgCode() );
          reco_beam_truth_daughter_true_IDs.push_back( daughterParticleFromRecoTrack->TrackId() );
          reco_beam_truth_daughter_true_lens.push_back( daughterParticleFromRecoTrack->Trajectory().TotalLength() );
        }
      }
      else{
        reco_beam_truth_daughter_good_reco.push_back( false );
        reco_beam_truth_daughter_true_PDGs.push_back( -1 );
        reco_beam_truth_daughter_true_IDs.push_back( -1 );
        reco_beam_truth_daughter_true_lens.push_back( -1 );
      }
        
    }

    nShowerDaughters = showerDaughters.size();

    for( size_t i = 0; i < showerDaughters.size(); ++i ){
      std::cout << "Shower daughter " << i << " Starts at " << showerDaughters[i]->ShowerStart().X() << " " << showerDaughters[i]->ShowerStart().Y() << " " << showerDaughters[i]->ShowerStart().Z() << std::endl;
      reco_daughter_showerID.push_back( showerDaughters[i]->ID() );
      reco_daughter_shower_startX.push_back( showerDaughters[i]->ShowerStart().X() );
      reco_daughter_shower_startY.push_back( showerDaughters[i]->ShowerStart().Y() );
      reco_daughter_shower_startZ.push_back( showerDaughters[i]->ShowerStart().Z() );

      auto daughterShowerFromRecoTrack = showerDaughters[i];

      bool found_daughter = false;
      const simb::MCParticle* daughterParticleFromRecoShower = truthUtil.GetMCParticleFromRecoShower(*daughterShowerFromRecoTrack, evt, fShowerTag); 
      std::cout << "Shower daughter " << i << " " << daughterParticleFromRecoShower << std::endl;

      if( daughterParticleFromRecoShower ){
        std::cout << "Shower daughter ID: " << daughterParticleFromRecoShower->TrackId() << std::endl;
        int loc = 0;
        for( size_t j = 0; j < true_beam_daughter_IDs.size(); ++j ){
          if ( true_beam_daughter_IDs[j] == daughterParticleFromRecoShower->TrackId() ){
            found_daughter = true;
            loc = j;
            break;
          }
        }

        reco_beam_truth_daughter_shower_good_reco.push_back( found_daughter );
        if( found_daughter ){
          reco_beam_truth_daughter_shower_true_PDGs.push_back( true_beam_daughter_PDGs[loc] ); 
          reco_beam_truth_daughter_shower_true_IDs.push_back( true_beam_daughter_IDs[loc] );
          reco_beam_truth_daughter_shower_true_lens.push_back( true_beam_daughter_lens[loc] );
        }
        else{
          reco_beam_truth_daughter_shower_true_PDGs.push_back( daughterParticleFromRecoShower->PdgCode() );
          reco_beam_truth_daughter_shower_true_IDs.push_back( daughterParticleFromRecoShower->TrackId() );
          reco_beam_truth_daughter_shower_true_lens.push_back( daughterParticleFromRecoShower->Trajectory().TotalLength() );
        }
      }

    }

  }
  else if( thisShower ){
    beamTrackID = thisShower->ID();
  }



  fTree->Fill();
}

void pionana::PionAnalyzerMC::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("type", &type);
  fTree->Branch("MC", &MC);

  fTree->Branch("startX", &startX);
  fTree->Branch("startY", &startY);
  fTree->Branch("startZ", &startZ);
  fTree->Branch("endX", &endX);
  fTree->Branch("endY", &endY);
  fTree->Branch("endZ", &endZ);
  fTree->Branch("len", &len);
  fTree->Branch("trackDirX", &trackDirX);
  fTree->Branch("trackDirY", &trackDirY);
  fTree->Branch("trackDirZ", &trackDirZ);
  fTree->Branch("vtxX", &vtxX);
  fTree->Branch("vtxY", &vtxY);
  fTree->Branch("vtxZ", &vtxZ);
  fTree->Branch("beamTrackID", &beamTrackID);
  fTree->Branch("dQdX", &dQdX);
  fTree->Branch("dEdX", &dEdX);
  fTree->Branch("resRange", &resRange);
  fTree->Branch("nBeamParticles", &nBeamParticles);


  fTree->Branch("reco_daughter_trackID", &reco_daughter_trackID);
  fTree->Branch("reco_daughter_showerID", &reco_daughter_showerID);
  fTree->Branch("reco_daughter_dQdX", &reco_daughter_dQdX);
  fTree->Branch("reco_daughter_dEdX", &reco_daughter_dEdX);
  fTree->Branch("reco_daughter_resRange", &reco_daughter_resRange);
  fTree->Branch("reco_daughter_len", &reco_daughter_len);
  fTree->Branch("reco_daughter_startX", &reco_daughter_startX);
  fTree->Branch("reco_daughter_startY", &reco_daughter_startY);
  fTree->Branch("reco_daughter_startZ", &reco_daughter_startZ);
  fTree->Branch("reco_daughter_endX", &reco_daughter_endX);
  fTree->Branch("reco_daughter_endY", &reco_daughter_endY);
  fTree->Branch("reco_daughter_endZ", &reco_daughter_endZ);
  fTree->Branch("reco_daughter_shower_startX", &reco_daughter_shower_startX);
  fTree->Branch("reco_daughter_shower_startY", &reco_daughter_shower_startY);
  fTree->Branch("reco_daughter_shower_startZ", &reco_daughter_shower_startZ);
  fTree->Branch("nTrackDaughters", &nTrackDaughters);
  fTree->Branch("nShowerDaughters", &nShowerDaughters);

  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_EndProcess", &true_beam_EndProcess);
  fTree->Branch("true_beam_EndVertex_X", &true_beam_EndVertex_X);
  fTree->Branch("true_beam_EndVertex_Y", &true_beam_EndVertex_Y);
  fTree->Branch("true_beam_EndVertex_Z", &true_beam_EndVertex_Z);
  fTree->Branch("true_beam_Start_X", &true_beam_Start_X);
  fTree->Branch("true_beam_Start_Y", &true_beam_Start_Y);
  fTree->Branch("true_beam_Start_Z", &true_beam_Start_Z);
  fTree->Branch("nPi0_truth", &nPi0_truth);
  fTree->Branch("nPiPlus_truth", &nPiPlus_truth);
  fTree->Branch("nProton_truth", &nProton_truth);
  fTree->Branch("nNeutron_truth", &nNeutron_truth);
  fTree->Branch("nPiMinus_truth", &nPiMinus_truth);

  fTree->Branch("true_beam_daughter_PDGs", &true_beam_daughter_PDGs);
  fTree->Branch("true_beam_daughter_IDs", &true_beam_daughter_IDs);
  fTree->Branch("true_beam_daughter_lens", &true_beam_daughter_lens);

  fTree->Branch("reco_beam_truth_EndProcess", &reco_beam_truth_EndProcess);
  fTree->Branch("reco_beam_truth_Process", &reco_beam_truth_Process);
  fTree->Branch("reco_beam_truth_origin", &reco_beam_truth_origin);
  fTree->Branch("reco_beam_truth_PDG", &reco_beam_truth_PDG);
  fTree->Branch("reco_beam_good", &reco_beam_good);


  fTree->Branch("reco_beam_truth_daughter_good_reco", &reco_beam_truth_daughter_good_reco);
  fTree->Branch("reco_beam_truth_daughter_true_PDGs", &reco_beam_truth_daughter_true_PDGs);
  fTree->Branch("reco_beam_truth_daughter_true_IDs", &reco_beam_truth_daughter_true_IDs);

  fTree->Branch("reco_beam_truth_daughter_shower_good_reco", &reco_beam_truth_daughter_shower_good_reco);
  fTree->Branch("reco_beam_truth_daughter_shower_true_PDGs", &reco_beam_truth_daughter_shower_true_PDGs);
  fTree->Branch("reco_beam_truth_daughter_shower_true_IDs", &reco_beam_truth_daughter_shower_true_IDs);

  fTree->Branch("reco_beam_truth_daughter_true_lens", &reco_beam_truth_daughter_true_lens);
  fTree->Branch("reco_beam_truth_daughter_shower_true_lens", &reco_beam_truth_daughter_shower_true_lens);
  fTree->Branch("reco_beam_truth_daughter_shower_good_reco", &reco_beam_truth_daughter_shower_good_reco);


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
  type = -1;
  nBeamParticles = 0;

  MC = 0;
  nProton_truth = 0;
  nNeutron_truth = 0;
  nPi0_truth = 0;
  nPiPlus_truth = 0;
  nPiMinus_truth = 0;
  reco_beam_truth_PDG = 0;
  true_beam_PDG = 0;
  true_beam_EndProcess ="";
  true_beam_EndVertex_X = 0.;
  true_beam_EndVertex_Y = 0.;
  true_beam_EndVertex_Z = 0.;
  true_beam_Start_X = 0.;
  true_beam_Start_Y = 0.;
  true_beam_Start_Z = 0.;
  reco_beam_truth_EndProcess ="";
  reco_beam_truth_Process ="";
  reco_beam_truth_origin = -1;

  reco_beam_good = false;
  reco_beam_truth_daughter_good_reco.clear();
  reco_beam_truth_daughter_true_PDGs.clear();
  reco_beam_truth_daughter_true_IDs.clear();

  reco_beam_truth_daughter_shower_good_reco.clear();
  reco_beam_truth_daughter_shower_true_PDGs.clear();
  reco_beam_truth_daughter_shower_true_IDs.clear();

  true_beam_daughter_PDGs.clear();
  true_beam_daughter_lens.clear();
  reco_beam_truth_daughter_true_lens.clear();
  reco_beam_truth_daughter_shower_true_lens.clear();
  true_beam_daughter_IDs.clear();
  reco_beam_truth_daughter_shower_good_reco.clear();

  nTrackDaughters = -1;
  nShowerDaughters = -1;

  dQdX.clear();
  dEdX.clear();
  vtxX = -1.;
  vtxY = -1.;
  vtxZ = -1.;
  reco_daughter_startX.clear();
  reco_daughter_startY.clear();
  reco_daughter_startZ.clear();
  reco_daughter_endX.clear();
  reco_daughter_endY.clear();
  reco_daughter_endZ.clear();

  reco_daughter_shower_startX.clear();
  reco_daughter_shower_startY.clear();
  reco_daughter_shower_startZ.clear();

  resRange.clear();

  reco_daughter_dQdX.clear();
  reco_daughter_dEdX.clear();
  reco_daughter_resRange.clear();
  reco_daughter_len.clear();

  beamTrackID = -1;
  reco_daughter_trackID.clear();
  reco_daughter_showerID.clear();

}

DEFINE_ART_MODULE(pionana::PionAnalyzerMC)
