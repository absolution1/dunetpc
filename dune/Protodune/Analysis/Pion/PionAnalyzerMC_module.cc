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
#include "lardataobj/RecoBase/TrackHitMeta.h"
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

#include "lardata/ArtDataHelper/MVAReader.h"

#include "art_root_io/TFileService.h"
#include "TProfile.h"
#include "TFile.h"

// ROOT includes
#include "TTree.h"

namespace pionana {
  class PionAnalyzerMC;
  bool sort_IDEs( const sim::IDE * i1, const sim::IDE * i2){
    return( i1->z < i2->z ); 
  }
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
  int true_beam_ID;
  std::string true_beam_EndProcess;
  double true_beam_EndVertex_X;
  double true_beam_EndVertex_Y;
  double true_beam_EndVertex_Z;
  double true_beam_Start_X;
  double true_beam_Start_Y;
  double true_beam_Start_Z;

  double true_beam_Start_DirX;
  double true_beam_Start_DirY;
  double true_beam_Start_DirZ;

  double true_beam_Start_Px;
  double true_beam_Start_Py;
  double true_beam_Start_Pz;
  double true_beam_Start_P;

    //Truth level info of the daughter MCParticles coming out of the 
    //true primary particle
  std::vector< int > true_beam_daughter_PDGs;
  std::vector< int > true_beam_daughter_IDs;
  std::vector< double > true_beam_daughter_lens;


  //Decay products from pi0s
  std::vector< int > true_beam_Pi0_decay_PDGs, true_beam_Pi0_decay_IDs;
  std::vector< int > true_beam_grand_daughter_PDGs, true_beam_grand_daughter_IDs, true_beam_grand_daughter_ParIDs;

  //How many of each true particle came out of the true primary beam particle?
  int nPiPlus_truth, nPiMinus_truth, nPi0_truth;
  int nProton_truth, nNeutron_truth, nNucleus_truth;
  /*****************************/
  



  //Reconstructed track info
  double startX, startY, startZ;
  double endX, endY, endZ;
  double vtxX, vtxY, vtxZ; 
  double len;
  double trackDirX, trackDirY, trackDirZ;
  double trackEndDirX, trackEndDirY, trackEndDirZ;
  std::vector< double > dEdX, dQdX, resRange;
  int beamTrackID;

  std::string reco_beam_truth_EndProcess; //What process ended the reco beam particle
  std::string reco_beam_truth_Process;    //What process created the reco beam particle
  
  int reco_beam_truth_PDG; //What is the PDG of the true MC particle contributing the most to 
                           //the reconstructed beam track
  int reco_beam_truth_ID;
  
  bool reco_beam_good; //Does the true particle contributing most to the 
                       //reconstructed beam track coincide with the actual
                       //beam particle that generated the event
                       
  bool elastic_candidate, found_elastic_scatter;
  int  nElasticScatters;
  std::vector< double > elastic_costheta, elastic_X, elastic_Y, elastic_Z;
  int nReconstructed;
  std::vector< int > reconstructedIDs;
  bool daughter_is_primary;

  int reco_beam_truth_origin; //What is the origin of the reconstructed beam track?

  double reco_beam_Chi2_proton;
  int    reco_beam_Chi2_ndof;

  std::vector< double > reco_daughter_Chi2_proton;
  std::vector< int >    reco_daughter_Chi2_ndof;
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

  double reco_beam_truth_End_Px;
  double reco_beam_truth_End_Py;
  double reco_beam_truth_End_Pz;
  double reco_beam_truth_End_E;
  double reco_beam_truth_End_P;

  double reco_beam_truth_Start_Px;
  double reco_beam_truth_Start_Py;
  double reco_beam_truth_Start_Pz;
  double reco_beam_truth_Start_E;
  double reco_beam_truth_Start_P;


  //Info from matching IDEs
  double IDE_toEnd;
  double reco_toEnd;
  bool   found_ides; 
  double total_traj_length;


  //Reco-level info of the reconstructed daughters coming out of the
  //reconstructed beam tracl
  std::vector< int > reco_daughter_trackID;
  std::vector< double > reco_daughter_completeness;
  std::vector< int > reco_daughter_truth_PDG;
  std::vector< int > reco_daughter_truth_ID;
  std::vector< int > reco_daughter_truth_Origin;
  std::vector< int > reco_daughter_truth_ParID;
  std::vector< std::string > reco_daughter_truth_Process;
  std::vector< int > reco_daughter_showerID;
  std::vector< int > reco_daughter_shower_truth_PDG;
  std::vector< int > reco_daughter_shower_truth_ID;
  std::vector< int > reco_daughter_shower_truth_Origin;
  std::vector< int > reco_daughter_shower_truth_ParID;
  std::vector< std::vector< double > > reco_daughter_dEdX, reco_daughter_dQdX, reco_daughter_resRange;
  std::vector< std::vector< double > > reco_daughter_shower_dEdX, reco_daughter_shower_dQdX, reco_daughter_shower_resRange;
  std::vector< double > reco_daughter_startX, reco_daughter_endX;
  std::vector< double > reco_daughter_startY, reco_daughter_endY;
  std::vector< double > reco_daughter_startZ, reco_daughter_endZ;
  std::vector< double > reco_daughter_deltaR;
  std::vector< double > reco_daughter_shower_startX;
  std::vector< double > reco_daughter_shower_startY;
  std::vector< double > reco_daughter_shower_startZ;
  std::vector< double > reco_daughter_shower_len;
  std::vector< double > reco_daughter_len;
  std::vector< double > reco_daughter_track_score;
  std::vector< double > reco_daughter_em_score;
  std::vector< double > reco_daughter_none_score;
  std::vector< double > reco_daughter_michel_score;

  std::vector< double > reco_daughter_shower_track_score;
  std::vector< double > reco_daughter_shower_em_score;
  std::vector< double > reco_daughter_shower_none_score;
  std::vector< double > reco_daughter_shower_michel_score;

  int nTrackDaughters, nShowerDaughters;

  int type;
  int nBeamParticles;

  std::map< int, TProfile* > templates;

  //FCL pars
  std::string fCalorimetryTag;
  std::string fTrackerTag;    
  std::string fHitTag;    
  std::string fShowerTag;     
  std::string fPFParticleTag; 
  std::string fGeneratorTag;
  std::string dEdX_template_name;
  TFile dEdX_template_file;
  bool fVerbose;             
  fhicl::ParameterSet dataUtil;
};


pionana::PionAnalyzerMC::PionAnalyzerMC(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),

  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  dEdX_template_name(p.get<std::string>("dEdX_template_name")),
  dEdX_template_file( dEdX_template_name.c_str(), "OPEN" ),
  fVerbose(p.get<bool>("Verbose")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils"))
{

  templates[ 211 ]  = (TProfile*)dEdX_template_file.Get( "dedx_range_pi"  );
  templates[ 321 ]  = (TProfile*)dEdX_template_file.Get( "dedx_range_ka"  );
  templates[ 13 ]   = (TProfile*)dEdX_template_file.Get( "dedx_range_mu"  );
  templates[ 2212 ] = (TProfile*)dEdX_template_file.Get( "dedx_range_pro" );

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
  protoana::ProtoDUNEShowerUtils                        showerUtil;
  art::ServiceHandle<cheat::BackTrackerService>         bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  protoana::ProtoDUNETruthUtils                         truthUtil;
  ////////////////////////////////////////
  


  // This gets the true beam particle that generated the event
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  const simb::MCParticle* true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
  if( !true_beam_particle ){
    std::cout << "No true beam particle" << std::endl;
    return;
  }
  ////////////////////////////
  
  

  // Helper to get hits and the 4 associated CNN outputs
  // CNN Outputs: EM, Track, Michel, Empty
  // outputNames: track, em, none, michel
  anab::MVAReader<recob::Hit,4> hitResults(evt, /*fNNetModuleLabel*/ "emtrkmichelid:emtrkmichel" );


  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,fTrackerTag);

  auto recoShowers = evt.getValidHandle< std::vector< recob::Shower > >(fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);

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


  std::cout << "particle " << particle << std::endl;

  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
  const simb::MCParticle* trueParticle = 0x0;

  std::cout << "thisTrack " << thisTrack << std::endl;
  std::cout << "thisShower " << thisShower << std::endl;
  std::cout << "trueParticle " << trueParticle << std::endl;


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

    reco_beam_truth_PDG = trueParticle->PdgCode();
    reco_beam_truth_ID = trueParticle->TrackId();

    reco_beam_truth_Process = trueParticle->Process();
    reco_beam_truth_EndProcess = trueParticle->EndProcess();
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


    reco_beam_truth_Start_Px = trueParticle->Px();
    reco_beam_truth_Start_Py = trueParticle->Py();
    reco_beam_truth_Start_Pz = trueParticle->Pz();
    reco_beam_truth_Start_P  = sqrt( reco_beam_truth_Start_Px*reco_beam_truth_Start_Px 
                                   + reco_beam_truth_Start_Py*reco_beam_truth_Start_Py 
                                   + reco_beam_truth_Start_Pz*reco_beam_truth_Start_Pz );
    reco_beam_truth_Start_E = trueParticle->E();

    size_t np = trueParticle->NumberTrajectoryPoints();
    if( np > 1 ){
      reco_beam_truth_End_Px = trueParticle->Px( np - 2 );
      reco_beam_truth_End_Py = trueParticle->Py( np - 2 );
      reco_beam_truth_End_Pz = trueParticle->Pz( np - 2 );
      reco_beam_truth_End_P  = sqrt( reco_beam_truth_End_Px*reco_beam_truth_End_Px 
                                   + reco_beam_truth_End_Py*reco_beam_truth_End_Py 
                                   + reco_beam_truth_End_Pz*reco_beam_truth_End_Pz );
      reco_beam_truth_End_E  = trueParticle->E( np - 2 );
    }

  }
  //////////////////////////////////////////////////////////////




  //Some truth information
  true_beam_EndProcess = true_beam_particle->EndProcess();
  
  true_beam_PDG         = true_beam_particle->PdgCode();
  true_beam_ID          = true_beam_particle->TrackId();
  true_beam_EndVertex_X = true_beam_particle->EndX();
  true_beam_EndVertex_Y = true_beam_particle->EndY();
  true_beam_EndVertex_Z = true_beam_particle->EndZ();
  true_beam_Start_X     = true_beam_particle->Position(0).X();
  true_beam_Start_Y     = true_beam_particle->Position(0).Y();
  true_beam_Start_Z     = true_beam_particle->Position(0).Z();

  true_beam_Start_Px    = true_beam_particle->Px();
  true_beam_Start_Py    = true_beam_particle->Py();
  true_beam_Start_Pz    = true_beam_particle->Pz();

  true_beam_Start_P     = true_beam_particle->P();

  true_beam_Start_DirX  = true_beam_Start_Px / true_beam_Start_P;
  true_beam_Start_DirY  = true_beam_Start_Py / true_beam_Start_P;
  true_beam_Start_DirZ  = true_beam_Start_Pz / true_beam_Start_P;
  

  //Go through and look at the reconstructed track->true IDs. Look for a match to the true
  //incident beam particle
  for( auto const & tr : *recoTracks ){
   const simb::MCParticle* true_p = truthUtil.GetMCParticleFromRecoTrack(tr, evt, fTrackerTag); 

   if(true_p){
     if( true_p->TrackId() == true_beam_particle->TrackId() ){
       ++nReconstructed;
       reconstructedIDs.push_back( tr.ID() );
     }
   }
  }



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

    if( pid == 211  ) ++nPiPlus_truth;
    if( pid == -211 ) ++nPiMinus_truth;
    if( pid == 111  ) ++nPi0_truth;
    if( pid == 2212 ) ++nProton_truth;
    if( pid == 2112 ) ++nNeutron_truth;
    if( pid > 2212  ) ++nNucleus_truth; 

    //Look for the gammas coming out of the pi0s
    if( pid == 111 ){
      std::cout << "Found pi0. Looking at true daughters" << std::endl;
      for( int j = 0; j < part->NumberDaughters(); ++j ){
        int pi0_decay_daughter_ID = part->Daughter(j);
        auto pi0_decay_part = plist[ pi0_decay_daughter_ID ];
        true_beam_Pi0_decay_PDGs.push_back( pi0_decay_part->PdgCode() );
        true_beam_Pi0_decay_IDs.push_back( pi0_decay_part->TrackId() );
      }
    }

    for( int j = 0; j < part->NumberDaughters(); ++j ){
      int grand_daughter_ID = part->Daughter(j);
      auto grand_daughter_part = plist[ grand_daughter_ID ];
      true_beam_grand_daughter_PDGs.push_back( grand_daughter_part->PdgCode() );
      true_beam_grand_daughter_IDs.push_back(  grand_daughter_part->TrackId() );
      true_beam_grand_daughter_ParIDs.push_back(  part->TrackId() );
    }

  }

  if( thisTrack ){

    beamTrackID = thisTrack->ID();

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

      trackEndDirX =  -1. * startDir.X(); 
      trackEndDirY =  -1. * startDir.Y(); 
      trackEndDirZ =  -1. * startDir.Z(); 
    }
    else{
      std::cout << "endZ > startZ: " << startZ << " " << endZ << std::endl;
      trackDirX    =  startDir.X(); 
      trackDirY    =  startDir.Y(); 
      trackDirZ    =  startDir.Z(); 
      trackEndDirX =  endDir.X(); 
      trackEndDirY =  endDir.Y(); 
      trackEndDirZ =  endDir.Z(); 
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


    //Go through the true processes within the MCTrajectory
    const simb::MCTrajectory & true_beam_trajectory = true_beam_particle->Trajectory();
    auto true_beam_proc_map = true_beam_trajectory.TrajectoryProcesses();
    std::cout << "Processes: " << std::endl;

    std::vector< size_t > elastic_indices;

    for( auto itProc = true_beam_proc_map.begin(); itProc != true_beam_proc_map.end(); ++itProc ){
      int index = itProc->first;
      std::string process = true_beam_trajectory.KeyToProcess(itProc->second);
      std::cout << index << " " << process << std::endl;

      if( process == "hadElastic" ){
        elastic_indices.push_back( index );

        ++nElasticScatters;

        double process_X = true_beam_trajectory.X( index );
        double process_Y = true_beam_trajectory.Y( index );
        double process_Z = true_beam_trajectory.Z( index );

        std::cout << "At "  
                  << process_X << " " 
                  << process_Y << " " 
                  << process_Z << std::endl;

        double delta = sqrt( 
          (process_X - endX)*(process_X - endX) +
          (process_Y - endY)*(process_Y - endY) +
          (process_Z - endZ)*(process_Z - endZ)
        );

        double PX      = true_beam_trajectory.Px( index );
        double next_PX = true_beam_trajectory.Px( index + 1 );
        double PY      = true_beam_trajectory.Py( index );
        double next_PY = true_beam_trajectory.Py( index + 1 );
        double PZ      = true_beam_trajectory.Pz( index );
        double next_PZ = true_beam_trajectory.Pz( index + 1 );

        double total_P = sqrt( PX*PX + PY*PY + PZ*PZ );
        double total_next_P = sqrt( next_PX*next_PX + next_PY*next_PY + next_PZ*next_PZ );

        //Get the angle between the direction of this step and the next
        elastic_costheta.push_back(
          ( ( PX * next_PX ) + ( PY * next_PY ) + ( PZ * next_PZ ) ) / ( total_P * total_next_P )
        );

        elastic_X.push_back( process_X );
        elastic_Y.push_back( process_Y );
        elastic_Z.push_back( process_Z );

        //for now: 10cm
        if( delta < 10. ){
          std::cout << "Found possible elastic scatter with vertex" << std::endl;
          elastic_candidate = true;
        }
      }
    }


    //Testing
    std::cout << "N Traj Pts: " << true_beam_trajectory.size() << std::endl;
    std::cout << true_beam_particle->TrackId() << std::endl;

    std::cout << "N Reco Traj Pts: " << thisTrack->NumberTrajectoryPoints() << std::endl;

    //Get the sim::IDEs for each true trajectory point
    std::map< size_t, std::vector< const sim::IDE * > > trueTrajPtsToSimIDEs = truthUtil.GetSimIDEs( *true_beam_particle );

    //Looking at the hits in the beam track
    //std::map< size_t, std::vector< const recob::Hit * > > trajPtsToHits = trackUtil.GetRecoHitsFromTrajPoints( *thisTrack, evt, fTrackerTag );
    std::map< size_t, const recob::Hit * > trajPtsToHits = trackUtil.GetRecoHitsFromTrajPoints( *thisTrack, evt, fTrackerTag );


    // Get the IDEs from last hit in the track
    auto reco_it = trajPtsToHits.rbegin();
    //auto last_hit = reco_it->second[0];
    std::cout << "Last index " << reco_it->first << std::endl;

    //Finding distance between end of the reconstucted track and the end of the true trajectory
    //by looking at the sim::IDEs

    //Loop backward over the hits in the reconstructed track
    //Note: the loops have the extra condition of occuring while 
    //the IDEs have not been matched
    const sim::IDE * matchedIDE = NULL;  //Not sure if I need these
    const recob::Hit * matchedHit = NULL;

    //For the true sim IDEs: save the index for the traj pt (TruePt_index) and the IDE in the vector (IDE_index)
    //For the reco hits: just save the index of the reco traj point; 
    size_t IDE_index = 0; 
    size_t TruePt_index = 0;
    size_t Hit_index = 0;
    std::map< size_t, std::vector< const sim::IDE * > >::iterator matchedTruePt = trueTrajPtsToSimIDEs.end();
    for( auto itRecoHits = trajPtsToHits.rbegin(); ( itRecoHits != trajPtsToHits.rend() && !found_ides ); ++itRecoHits ){
      //auto theHit = itRecoHits->second[0];      
      auto theHit = itRecoHits->second;
      Hit_index = itRecoHits->first;
        
      std::vector< const sim::IDE * > recoHitIDEs = bt_serv->HitToSimIDEs_Ps( *theHit );     
      
      //Now go backward through the pts in the true trajectory and try to find a match between sim::IDEs
      for( auto itTrueIDEs = trueTrajPtsToSimIDEs.rbegin(); ( itTrueIDEs != trueTrajPtsToSimIDEs.rend() && !found_ides ); ++itTrueIDEs ){
        std::vector< const sim::IDE * > trueIDEs = itTrueIDEs->second;
        TruePt_index = itTrueIDEs->first;

        //Skip empty points
        if( !trueIDEs.size() ) continue;

        //Go backward through the true IDEs, then look at the reco IDEs
        //Note: this for loop uses unsigned ints, so when --i tries to 
        //decrement 0, it will go to max<unsigned int> which is > trueIDEs.size()
        //and will exit the loop
        for( size_t i = trueIDEs.size() - 1; ( i < trueIDEs.size() && !found_ides ); --i ){
          const sim::IDE * theTrueIDE = trueIDEs[i];
          
          for( size_t j = 0; j < recoHitIDEs.size(); ++j ){
            //Check if the pointers are the same value
            if( recoHitIDEs[j] == theTrueIDE ){
              matchedIDE = theTrueIDE; 
              matchedHit = theHit;
              IDE_index  = i;
              found_ides = true; 
              matchedTruePt = trueTrajPtsToSimIDEs.find( TruePt_index );
              break;
            }
          }
        }
      }
    }

    if( found_ides ){
      std::cout << "Found the IDEs " << Hit_index << " " << TruePt_index << " " << IDE_index << std::endl;
      std::cout << matchedIDE << " " << matchedHit << std::endl;

      //Get the distance to the end of the reco track
      reco_toEnd = thisTrack->Length( Hit_index ); 
      std::cout << "To end of reco track: " << reco_toEnd << std::endl;

      //Calculate the integrated distance (along the trajectory) from the matched IDE point to the end of the trajectory
      //using only the trajectory points
      double IDE_x = matchedIDE->x;
      double IDE_y = matchedIDE->y;
      double IDE_z = matchedIDE->z;

      if( true_beam_trajectory.size() > ( TruePt_index + 1 ) ){

        IDE_toEnd += sqrt( std::pow( ( IDE_x - true_beam_trajectory.X( TruePt_index + 1 ) ), 2 ) + 
                           std::pow( ( IDE_y - true_beam_trajectory.Y( TruePt_index + 1 ) ), 2 ) + 
                           std::pow( ( IDE_z - true_beam_trajectory.Z( TruePt_index + 1 ) ), 2 ) );

        for( size_t i = TruePt_index + 1; i < true_beam_trajectory.size() - 1; ++i ){

          double delta = sqrt( std::pow( ( true_beam_trajectory.X( i ) - true_beam_trajectory.X( i + 1 ) ), 2 ) + 
                               std::pow( ( true_beam_trajectory.Y( i ) - true_beam_trajectory.Y( i + 1 ) ), 2 ) + 
                               std::pow( ( true_beam_trajectory.Z( i ) - true_beam_trajectory.Z( i + 1 ) ), 2 ) );

          IDE_toEnd += delta;
        }
      }
      
      std::cout << "To end of true trajectory: " << IDE_toEnd << std::endl;
      
      
      //total length is the length starting from where the sim::IDE began
      for( size_t i = trueTrajPtsToSimIDEs.begin()->first; i < true_beam_trajectory.size() - 1; ++i ){
        total_traj_length += sqrt( std::pow( ( true_beam_trajectory.X( i ) - true_beam_trajectory.X( i + 1 ) ), 2 )  +
                                   std::pow( ( true_beam_trajectory.Y( i ) - true_beam_trajectory.Y( i + 1 ) ), 2 )  +
                                   std::pow( ( true_beam_trajectory.Z( i ) - true_beam_trajectory.Z( i + 1 ) ), 2 )  );
      }
      
      std::cout << "Total length: " << total_traj_length << std::endl;


      

    }
    else{
      std::cout << "Did not find IDEs " << matchedIDE << " " << matchedHit << std::endl; 
    }

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

      double delta_start = sqrt( ( reco_daughter_startX.back() - endX )*( reco_daughter_startX.back() - endX ) + 
                                 ( reco_daughter_startY.back() - endY )*( reco_daughter_startY.back() - endY ) +
                                 ( reco_daughter_startZ.back() - endZ )*( reco_daughter_startZ.back() - endZ ) );

      double delta_end =   sqrt( ( reco_daughter_endX.back() - endX )*( reco_daughter_endX.back() - endX ) + 
                                 ( reco_daughter_endY.back() - endY )*( reco_daughter_endY.back() - endY ) +
                                 ( reco_daughter_endZ.back() - endZ )*( reco_daughter_endZ.back() - endZ ) );

      reco_daughter_deltaR.push_back( ( delta_start < delta_end ) ? delta_start : delta_end );

      reco_daughter_trackID.push_back( daughterTrack->ID() );

      reco_daughter_completeness.push_back( truthUtil.GetCompleteness( *daughterTrack, evt, fTrackerTag, fHitTag ) );

      //For this reco daughter track, get the actual particle contributing to it
      const simb::MCParticle* trueDaughterParticle = truthUtil.GetMCParticleFromRecoTrack(*daughterTrack, evt, fTrackerTag);
      if( trueDaughterParticle ){
        reco_daughter_truth_PDG.push_back( trueDaughterParticle->PdgCode() );
        reco_daughter_truth_ID.push_back( trueDaughterParticle->TrackId() );
        if( trueDaughterParticle->TrackId() == true_beam_particle->TrackId() ){
          daughter_is_primary = true;
        }
          
        reco_daughter_truth_Origin.push_back( 
          pi_serv->TrackIdToMCTruth_P(trueDaughterParticle->TrackId())->Origin()
        );
        reco_daughter_truth_ParID.push_back( trueDaughterParticle->Mother() );
        reco_daughter_truth_Process.push_back( trueDaughterParticle->Process() );
      }
      else{
        reco_daughter_truth_PDG.push_back( -1 );
        reco_daughter_truth_ID.push_back( -1 );
        reco_daughter_truth_Origin.push_back( -1 );
        reco_daughter_truth_ParID.push_back( -1 );
        reco_daughter_truth_Process.push_back( "empty" );
      }

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

      std::pair< double, int > this_chi2_ndof = trackUtil.Chi2PID( reco_daughter_dEdX.back(), reco_daughter_resRange.back(), templates[ 2212 ] );
      reco_daughter_Chi2_proton.push_back( this_chi2_ndof.first );
      reco_daughter_Chi2_ndof.push_back( this_chi2_ndof.second );


      //Go through the hits of each daughter and check the CNN scores 
      //Look at the hits from the track:
      std::cout << "Getting hits from daughter" << std::endl;
      std::vector< art::Ptr< recob::Hit > > daughter_hits = findHits.at( daughterTrack->ID() );

      double track_total = 0.;
      double em_total = 0.;
      double none_total = 0.;
      double michel_total = 0.;   
      for( size_t h = 0; h < daughter_hits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults.getOutput( daughter_hits[h ] );
        track_total  += cnn_out[ hitResults.getIndex("track") ];
        em_total     += cnn_out[ hitResults.getIndex("em") ];
        none_total   += cnn_out[ hitResults.getIndex("none") ];
        michel_total += cnn_out[ hitResults.getIndex("michel") ];
      }

      if( daughter_hits.size() > 0 ){
        reco_daughter_track_score.push_back( track_total / daughter_hits.size() );
        reco_daughter_em_score.push_back( em_total / daughter_hits.size() );
        reco_daughter_none_score.push_back( none_total / daughter_hits.size() );
        reco_daughter_michel_score.push_back( michel_total / daughter_hits.size() );
      }
      else{
        reco_daughter_track_score.push_back( -999. );
        reco_daughter_em_score.push_back( -999. );
        reco_daughter_none_score.push_back( -999. );
        reco_daughter_michel_score.push_back( -999. );
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
      auto daughterShowerFromRecoTrack = showerDaughters[i];

      std::cout << "Shower daughter " << i << " Starts at " << daughterShowerFromRecoTrack->ShowerStart().X() << " " << daughterShowerFromRecoTrack->ShowerStart().Y() << " " << daughterShowerFromRecoTrack->ShowerStart().Z() << std::endl;
      reco_daughter_showerID.push_back( daughterShowerFromRecoTrack->ID() );
      reco_daughter_shower_startX.push_back( daughterShowerFromRecoTrack->ShowerStart().X() );
      reco_daughter_shower_startY.push_back( daughterShowerFromRecoTrack->ShowerStart().Y() );
      reco_daughter_shower_startZ.push_back( daughterShowerFromRecoTrack->ShowerStart().Z() );

      reco_daughter_shower_len.push_back( daughterShowerFromRecoTrack->Length() );

      std::cout << "Looking at shower dedx: ";
      std::cout << daughterShowerFromRecoTrack->dEdx().size() << std::endl;

      //For this reco daughter track, get the actual particle contributing to it
      const simb::MCParticle* trueDaughterParticle = truthUtil.GetMCParticleFromRecoShower(*daughterShowerFromRecoTrack, evt, fShowerTag);
      if( trueDaughterParticle ){
        reco_daughter_shower_truth_PDG.push_back( trueDaughterParticle->PdgCode() );
        reco_daughter_shower_truth_ID.push_back( trueDaughterParticle->TrackId() );
        reco_daughter_shower_truth_Origin.push_back( 
          pi_serv->TrackIdToMCTruth_P(trueDaughterParticle->TrackId())->Origin()
        );
        reco_daughter_shower_truth_ParID.push_back( trueDaughterParticle->Mother() );

      }
      else{
        reco_daughter_shower_truth_PDG.push_back( -1 );
        reco_daughter_shower_truth_ID.push_back( -1 );
        reco_daughter_shower_truth_Origin.push_back( -1 );
        reco_daughter_shower_truth_ParID.push_back( -1 );

      }

      std::vector< anab::Calorimetry > dummy_calo = showerUtil.GetRecoShowerCalorimetry(*daughterShowerFromRecoTrack, evt, fShowerTag, "showercalo");

      reco_daughter_shower_dQdX.push_back( std::vector<double>() );   
      reco_daughter_shower_resRange.push_back( std::vector<double>() );
      reco_daughter_shower_dEdX.push_back( std::vector<double>() );

      if( dummy_calo.size() ){
        std::cout << "Plane: " << dummy_calo[0].PlaneID().toString() << std::endl;      
        auto dummy_dQdx = dummy_calo[0].dQdx();
        auto dummy_dEdx = dummy_calo[0].dEdx();
        auto dummy_Range = dummy_calo[0].ResidualRange();
 
        for( size_t j = 0; j < dummy_dQdx.size(); ++j ){
          reco_daughter_shower_dQdX.back().push_back( dummy_dQdx[j] );
          reco_daughter_shower_resRange.back().push_back( dummy_Range[j] );
          reco_daughter_shower_dEdX.back().push_back( dummy_dEdx[j] );
        }
      }
      //Go through the hits of each daughter and check the CNN scores 
      //Look at the hits from the track:
      int shower_index = showerUtil.GetShowerIndex( *daughterShowerFromRecoTrack, evt, fShowerTag );
      std::vector< art::Ptr< recob::Hit > > daughter_hits = findHitsFromShowers.at( shower_index );
      std::cout << "Got " << daughter_hits.size() << " hits from shower " << daughterShowerFromRecoTrack->ID() << " " << shower_index << std::endl;

      double track_total = 0.;
      double em_total = 0.;
      double none_total = 0.;
      double michel_total = 0.;   
      for( size_t h = 0; h < daughter_hits.size(); ++h ){
        std::array<float,4> cnn_out = hitResults.getOutput( daughter_hits[h ] );
        track_total  += cnn_out[ hitResults.getIndex("track") ];
        em_total     += cnn_out[ hitResults.getIndex("em") ];
        none_total   += cnn_out[ hitResults.getIndex("none") ];
        michel_total += cnn_out[ hitResults.getIndex("michel") ];


        //In addition, try to get the dQ from each hit
        auto theHit = daughter_hits[ h ];
        double dQ = theHit->Integral();
        std::cout << "Hit: " << h << " has dQ " << dQ << std::endl;
        
      }

      if( daughter_hits.size() > 0 ){
        reco_daughter_shower_track_score.push_back( track_total / daughter_hits.size() );
        reco_daughter_shower_em_score.push_back( em_total / daughter_hits.size() );
        reco_daughter_shower_none_score.push_back( none_total / daughter_hits.size() );
        reco_daughter_shower_michel_score.push_back( michel_total / daughter_hits.size() );
      }
      else{
        reco_daughter_shower_track_score.push_back( -999. );
        reco_daughter_shower_em_score.push_back( -999. );
        reco_daughter_shower_none_score.push_back( -999. );
        reco_daughter_shower_michel_score.push_back( -999. );
      }
      



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

    std::pair< double, int > pid_chi2_ndof = trackUtil.Chi2PID( dEdX, resRange, templates[ 2212 ] );
    reco_beam_Chi2_proton = pid_chi2_ndof.first; 
    reco_beam_Chi2_ndof = pid_chi2_ndof.second;
  
    std::cout << "Proton chi2: " << reco_beam_Chi2_proton << std::endl;


    //Looking at shower/track discrimination
    std::cout << "MVA" << std::endl;

    double track_total = 0.;
    double em_total = 0.;
    double none_total = 0.;
    double michel_total = 0.;   
    std::vector< art::Ptr< recob::Hit > > track_hits = findHits.at( thisTrack->ID() );
    for( size_t h = 0; h < track_hits.size(); ++h ){
      std::array<float,4> cnn_out = hitResults.getOutput( track_hits[h ] );
      track_total  += cnn_out[ hitResults.getIndex("track") ];
      em_total     += cnn_out[ hitResults.getIndex("em") ];
      none_total   += cnn_out[ hitResults.getIndex("none") ];
      michel_total += cnn_out[ hitResults.getIndex("michel") ];
    }
    std::cout << "track Total "  << track_total  << " " << track_hits.size() << std::endl;
    std::cout << "em Total "     << em_total     << " " << track_hits.size() << std::endl;
    std::cout << "none Total "   << none_total   << " " << track_hits.size() << std::endl;
    std::cout << "michel Total " << michel_total << " " << track_hits.size() << std::endl;

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

  fTree->Branch("IDE_toEnd", &IDE_toEnd);
  fTree->Branch("found_ides", &found_ides); 
  fTree->Branch("total_traj_length", &total_traj_length);
  fTree->Branch("reco_toEnd", &reco_toEnd);

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
  fTree->Branch("trackEndDirX", &trackEndDirX);
  fTree->Branch("trackEndDirY", &trackEndDirY);
  fTree->Branch("trackEndDirZ", &trackEndDirZ);
  fTree->Branch("vtxX", &vtxX);
  fTree->Branch("vtxY", &vtxY);
  fTree->Branch("vtxZ", &vtxZ);
  fTree->Branch("beamTrackID", &beamTrackID);
  fTree->Branch("dQdX", &dQdX);
  fTree->Branch("dEdX", &dEdX);
  fTree->Branch("resRange", &resRange);
  fTree->Branch("nBeamParticles", &nBeamParticles);


  fTree->Branch("reco_daughter_trackID", &reco_daughter_trackID);
  fTree->Branch("reco_daughter_completeness", &reco_daughter_completeness);
  fTree->Branch("reco_daughter_truth_PDG", &reco_daughter_truth_PDG);
  fTree->Branch("reco_daughter_truth_ID", &reco_daughter_truth_ID);
  fTree->Branch("reco_daughter_truth_Origin", &reco_daughter_truth_Origin);
  fTree->Branch("reco_daughter_truth_ParID", &reco_daughter_truth_ParID);
  fTree->Branch("reco_daughter_truth_Process", &reco_daughter_truth_Process);
  fTree->Branch("reco_daughter_shower_truth_PDG", &reco_daughter_shower_truth_PDG);
  fTree->Branch("reco_daughter_shower_truth_ID", &reco_daughter_shower_truth_ID);
  fTree->Branch("reco_daughter_shower_truth_Origin", &reco_daughter_shower_truth_Origin);
  fTree->Branch("reco_daughter_shower_truth_ParID", &reco_daughter_shower_truth_ParID);

  fTree->Branch("reco_daughter_showerID", &reco_daughter_showerID);
  fTree->Branch("reco_daughter_dQdX", &reco_daughter_dQdX);
  fTree->Branch("reco_daughter_dEdX", &reco_daughter_dEdX);
  fTree->Branch("reco_daughter_resRange", &reco_daughter_resRange);
  fTree->Branch("reco_daughter_shower_dQdX", &reco_daughter_shower_dQdX);
  fTree->Branch("reco_daughter_shower_dEdX", &reco_daughter_shower_dEdX);
  fTree->Branch("reco_daughter_shower_resRange", &reco_daughter_shower_resRange);
  fTree->Branch("reco_daughter_len", &reco_daughter_len);
  fTree->Branch("reco_daughter_startX", &reco_daughter_startX);
  fTree->Branch("reco_daughter_startY", &reco_daughter_startY);
  fTree->Branch("reco_daughter_startZ", &reco_daughter_startZ);
  fTree->Branch("reco_daughter_endX", &reco_daughter_endX);
  fTree->Branch("reco_daughter_endY", &reco_daughter_endY);
  fTree->Branch("reco_daughter_endZ", &reco_daughter_endZ);
  fTree->Branch("reco_daughter_deltaR", &reco_daughter_deltaR);
  fTree->Branch("reco_daughter_shower_startX", &reco_daughter_shower_startX);
  fTree->Branch("reco_daughter_shower_startY", &reco_daughter_shower_startY);
  fTree->Branch("reco_daughter_shower_startZ", &reco_daughter_shower_startZ);
  fTree->Branch("reco_daughter_shower_len", &reco_daughter_shower_len);
  fTree->Branch("nTrackDaughters", &nTrackDaughters);
  fTree->Branch("nShowerDaughters", &nShowerDaughters);

  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_ID", &true_beam_ID);
  fTree->Branch("true_beam_EndProcess", &true_beam_EndProcess);
  fTree->Branch("true_beam_EndVertex_X", &true_beam_EndVertex_X);
  fTree->Branch("true_beam_EndVertex_Y", &true_beam_EndVertex_Y);
  fTree->Branch("true_beam_EndVertex_Z", &true_beam_EndVertex_Z);
  fTree->Branch("true_beam_Start_X", &true_beam_Start_X);
  fTree->Branch("true_beam_Start_Y", &true_beam_Start_Y);
  fTree->Branch("true_beam_Start_Z", &true_beam_Start_Z);

  fTree->Branch("true_beam_Start_Px", &true_beam_Start_Px);
  fTree->Branch("true_beam_Start_Py", &true_beam_Start_Py);
  fTree->Branch("true_beam_Start_Pz", &true_beam_Start_Pz);
  fTree->Branch("true_beam_Start_P", &true_beam_Start_P);

  fTree->Branch("true_beam_Start_DirX", &true_beam_Start_DirX);
  fTree->Branch("true_beam_Start_DirY", &true_beam_Start_DirY);
  fTree->Branch("true_beam_Start_DirZ", &true_beam_Start_DirZ);

  fTree->Branch("nPi0_truth", &nPi0_truth);
  fTree->Branch("nPiPlus_truth", &nPiPlus_truth);
  fTree->Branch("nProton_truth", &nProton_truth);
  fTree->Branch("nNeutron_truth", &nNeutron_truth);
  fTree->Branch("nPiMinus_truth", &nPiMinus_truth);
  fTree->Branch("nNucleus_truth", &nNucleus_truth);

  fTree->Branch("true_beam_daughter_PDGs", &true_beam_daughter_PDGs);
  fTree->Branch("true_beam_daughter_IDs", &true_beam_daughter_IDs);
  fTree->Branch("true_beam_daughter_lens", &true_beam_daughter_lens);
  fTree->Branch("true_beam_Pi0_decay_IDs", &true_beam_Pi0_decay_IDs);
  fTree->Branch("true_beam_Pi0_decay_PDGs", &true_beam_Pi0_decay_PDGs);
  fTree->Branch("true_beam_grand_daughter_IDs", &true_beam_grand_daughter_IDs);
  fTree->Branch("true_beam_grand_daughter_ParIDs", &true_beam_grand_daughter_ParIDs);
  fTree->Branch("true_beam_grand_daughter_PDGs", &true_beam_grand_daughter_PDGs);

  fTree->Branch("reco_beam_truth_EndProcess", &reco_beam_truth_EndProcess);
  fTree->Branch("reco_beam_truth_Process", &reco_beam_truth_Process);
  fTree->Branch("reco_beam_truth_origin", &reco_beam_truth_origin);
  fTree->Branch("reco_beam_truth_PDG", &reco_beam_truth_PDG);
  fTree->Branch("reco_beam_truth_ID", &reco_beam_truth_ID);
  fTree->Branch("reco_beam_good", &reco_beam_good);
  fTree->Branch("elastic_candidate", &elastic_candidate);
  fTree->Branch("found_elastic_scatter", &found_elastic_scatter);
  fTree->Branch("nElasticScatters", &nElasticScatters);
  fTree->Branch("elastic_costheta", &elastic_costheta);
  fTree->Branch("elastic_X", &elastic_X);
  fTree->Branch("elastic_Y", &elastic_Y);
  fTree->Branch("elastic_Z", &elastic_Z);
  fTree->Branch("nReconstructed", &nReconstructed);
  fTree->Branch("reconstructedIDs", &reconstructedIDs);
  fTree->Branch("daughter_is_primary", &daughter_is_primary);

  fTree->Branch("reco_beam_Chi2_proton", &reco_beam_Chi2_proton);
  fTree->Branch("reco_beam_Chi2_ndof", &reco_beam_Chi2_ndof);

  fTree->Branch("reco_daughter_Chi2_proton", &reco_daughter_Chi2_proton);
  fTree->Branch("reco_daughter_Chi2_ndof", &reco_daughter_Chi2_ndof);

  fTree->Branch("reco_daughter_track_score", &reco_daughter_track_score);
  fTree->Branch("reco_daughter_em_score", &reco_daughter_em_score);
  fTree->Branch("reco_daughter_none_score", &reco_daughter_none_score);
  fTree->Branch("reco_daughter_michel_score", &reco_daughter_michel_score);

  fTree->Branch("reco_daughter_shower_track_score", &reco_daughter_shower_track_score);
  fTree->Branch("reco_daughter_shower_em_score", &reco_daughter_shower_em_score);
  fTree->Branch("reco_daughter_shower_none_score", &reco_daughter_shower_none_score);
  fTree->Branch("reco_daughter_shower_michel_score", &reco_daughter_shower_michel_score);

  fTree->Branch("reco_beam_truth_End_Px", &reco_beam_truth_End_Px);
  fTree->Branch("reco_beam_truth_End_Py", &reco_beam_truth_End_Py);
  fTree->Branch("reco_beam_truth_End_Pz", &reco_beam_truth_End_Pz);
  fTree->Branch("reco_beam_truth_End_E", &reco_beam_truth_End_E);
  fTree->Branch("reco_beam_truth_End_P", &reco_beam_truth_End_P);

  fTree->Branch("reco_beam_truth_Start_Px", &reco_beam_truth_Start_Px);
  fTree->Branch("reco_beam_truth_Start_Py", &reco_beam_truth_Start_Py);
  fTree->Branch("reco_beam_truth_Start_Pz", &reco_beam_truth_Start_Pz);
  fTree->Branch("reco_beam_truth_Start_E", &reco_beam_truth_Start_E);
  fTree->Branch("reco_beam_truth_Start_P", &reco_beam_truth_Start_P);

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
  dEdX_template_file.Close();
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

  IDE_toEnd = 0.;
  reco_toEnd = 0.;
  found_ides = false;
  total_traj_length = 0.;

  MC = 0;
  nProton_truth = 0;
  nNeutron_truth = 0;
  nNucleus_truth = 0;
  nPi0_truth = 0;
  nPiPlus_truth = 0;
  nPiMinus_truth = 0;
  reco_beam_truth_PDG = 0;
  reco_beam_truth_ID = 0;
  true_beam_PDG = 0;
  true_beam_ID = 0;
  true_beam_EndProcess ="";
  true_beam_EndVertex_X = 0.;
  true_beam_EndVertex_Y = 0.;
  true_beam_EndVertex_Z = 0.;
  true_beam_Start_X = 0.;
  true_beam_Start_Y = 0.;
  true_beam_Start_Z = 0.;

  true_beam_Start_Px   = 0.; 
  true_beam_Start_Py   = 0.; 
  true_beam_Start_Pz   = 0.; 

  true_beam_Start_P    = 0.; 

  true_beam_Start_DirX = 0.; 
  true_beam_Start_DirY = 0.; 
  true_beam_Start_DirZ = 0.; 

  reco_beam_truth_EndProcess ="";
  reco_beam_truth_Process ="";
  reco_beam_truth_origin = -1;

  reco_beam_truth_End_Px = 0.;
  reco_beam_truth_End_Py = 0.;
  reco_beam_truth_End_Pz = 0.;
  reco_beam_truth_End_E = 0.;
  reco_beam_truth_End_P = 0.;

  reco_beam_truth_Start_Px = 0.;
  reco_beam_truth_Start_Py = 0.;
  reco_beam_truth_Start_Pz = 0.;
  reco_beam_truth_Start_E = 0.;
  reco_beam_truth_Start_P = 0.;

  reco_beam_good = false;
  elastic_candidate = false;
  found_elastic_scatter = false;
  nElasticScatters = 0;
  elastic_costheta.clear();
  elastic_X.clear();
  elastic_Y.clear();
  elastic_Z.clear();
  nReconstructed = 0;
  reconstructedIDs.clear();
  daughter_is_primary = false;
  reco_beam_Chi2_proton = 999.;
  reco_beam_Chi2_ndof = -1;

  reco_daughter_Chi2_proton.clear();
  reco_daughter_Chi2_ndof.clear();

  reco_daughter_track_score.clear();
  reco_daughter_em_score.clear();
  reco_daughter_none_score.clear();
  reco_daughter_michel_score.clear();
  
  reco_daughter_shower_track_score.clear();
  reco_daughter_shower_em_score.clear();
  reco_daughter_shower_none_score.clear();
  reco_daughter_shower_michel_score.clear();

  reco_beam_truth_daughter_good_reco.clear();
  reco_beam_truth_daughter_true_PDGs.clear();
  reco_beam_truth_daughter_true_IDs.clear();

  reco_beam_truth_daughter_shower_good_reco.clear();
  reco_beam_truth_daughter_shower_true_PDGs.clear();
  reco_beam_truth_daughter_shower_true_IDs.clear();

  true_beam_daughter_PDGs.clear();
  true_beam_daughter_lens.clear();
  true_beam_Pi0_decay_IDs.clear();
  true_beam_Pi0_decay_PDGs.clear();
  true_beam_grand_daughter_IDs.clear();
  true_beam_grand_daughter_ParIDs.clear();
  true_beam_grand_daughter_PDGs.clear();
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
  reco_daughter_deltaR.clear();

  reco_daughter_shower_startX.clear();
  reco_daughter_shower_startY.clear();
  reco_daughter_shower_startZ.clear();

  reco_daughter_shower_len.clear(); 

  resRange.clear();

  reco_daughter_dQdX.clear();
  reco_daughter_dEdX.clear();
  reco_daughter_resRange.clear();
  reco_daughter_shower_dQdX.clear();
  reco_daughter_shower_dEdX.clear();
  reco_daughter_shower_resRange.clear();
  reco_daughter_len.clear();

  beamTrackID = -1;
  reco_daughter_trackID.clear();
  reco_daughter_completeness.clear();
  reco_daughter_truth_PDG.clear();
  reco_daughter_truth_ID.clear();
  reco_daughter_truth_Origin.clear();
  reco_daughter_truth_ParID.clear();
  reco_daughter_truth_Process.clear();

  reco_daughter_showerID.clear();
  reco_daughter_shower_truth_PDG.clear();
  reco_daughter_shower_truth_ID.clear();
  reco_daughter_shower_truth_Origin.clear();
  reco_daughter_shower_truth_ParID.clear();

}

DEFINE_ART_MODULE(pionana::PionAnalyzerMC)
