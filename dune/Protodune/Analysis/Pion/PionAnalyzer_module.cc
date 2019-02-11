////////////////////////////////////////////////////////////////////////
// Class:       PionAnalyzer
// Plugin Type: analyzer (art v3_00_00)
// File:        PionAnalyzer_module.cc
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
  class PionAnalyzer;
}


class pionana::PionAnalyzer : public art::EDAnalyzer {
public:
  explicit PionAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PionAnalyzer(PionAnalyzer const&) = delete;
  PionAnalyzer(PionAnalyzer&&) = delete;
  PionAnalyzer& operator=(PionAnalyzer const&) = delete;
  PionAnalyzer& operator=(PionAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();

private:

  // Declare member data here.
  const art::InputTag fBeamModuleLabel;
  const art::InputTag fTrackModuleLabel;

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;

  int MC;

  int nPiPlus_truth, nPiMinus_truth, nPi0_truth;
  int nProton_truth, nNeutron_truth;

  double startX, startY, startZ;
  double endX, endY, endZ;
  double len;
  double stitch_len;
  double combined_len;
  double beam_costheta;
  double new_beam_costheta;
  double beamDirX, beamDirY, beamDirZ;
  double trackDirX, trackDirY, trackDirZ;
  double newDirX, newDirY, newDirZ;
  int beamTrackID;
  std::vector< int >  stitchTrackID;
  std::vector< double > stitch_cos_theta;
  std::vector< double > stitch_dQdX;
  std::vector< double > stitch_dEdX;
  std::vector< double > stitch_resRange;
  std::vector< double > combined_dQdX;
  std::vector< double > combined_dEdX;
  std::vector< double > combined_resRange;

  std::vector< double > dEdX, dQdX, resRange;
  std::vector< std::vector< double > > daughter_dEdX, daughter_dQdX, daughter_resRange;
  std::vector< double > daughter_len;

  int type;

  std::string fCalorimetryTag;
  std::string fTrackerTag;    
  std::string fShowerTag;     
  std::string fPFParticleTag; 
  std::string fGeneratorTag;
  bool fVerbose;             
  fhicl::ParameterSet dataUtil;
  fhicl::ParameterSet beamlineUtil;
  double fBrokenTrackZ_low, fBrokenTrackZ_high;
  double fStitchTrackZ_low, fStitchTrackZ_high;
  double fStitchXTol, fStitchYTol;

};


pionana::PionAnalyzer::PionAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),

  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  beamlineUtil( p.get<fhicl::ParameterSet>("BeamlineUtils")),

  fBrokenTrackZ_low( p.get<double>("BrokenTrackZ_low") ),
  fBrokenTrackZ_high( p.get<double>("BrokenTrackZ_high") ),
  fStitchTrackZ_low( p.get<double>("StitchTrackZ_low") ),
  fStitchTrackZ_high( p.get<double>("StitchTrackZ_high") ),
  fStitchXTol( p.get<double>("StitchXTol") ),
  fStitchYTol( p.get<double>("StitchYTol") )

{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pionana::PionAnalyzer::analyze(art::Event const& evt)
{
  //reset containers
  reset();  


  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

  //Instantiate this here. Fill later if real data
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  
  //Beamline utils 
  protoana::ProtoDUNEBeamlineUtils BLUtil(beamlineUtil);
   
  // Get the PFParticle utility
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);

  if(beamParticles.size() == 0){
    std::cerr << "We found no beam particles for this event... moving on" << std::endl;
    return;
  }

  // We can now look at these particles
  for(const recob::PFParticle* particle : beamParticles){

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
      len  = thisTrack->Length();    
      beamTrackID = thisTrack->ID();
      
      std::cout << "Start: " << startX << " " << startY << " " << startZ << std::endl;
      std::cout << "End: " << endX << " " << endY << " " << endZ << std::endl;
      std::cout << "len: " << len << std::endl;

      if( evt.isRealData() ){
        //Get Beamline info
        auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamModuleLabel);
        if( beamHandle.isValid()){
          art::fill_ptr_vector(beamVec, beamHandle);
        }
        //Should just have one
        const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0));
      
        std::vector< recob::Track > newTracks = BLUtil.MakeTracks( evt );
        const std::vector< recob::Track > & beamTracks = beamEvent.GetBeamTracks();
        std::cout << "Beamline tracks: " << beamTracks.size() << std::endl;

        if( beamTracks.size() == 1 ){
          auto trackDir = thisTrack->StartDirection();
          auto beamDir = beamTracks.at(0).StartDirection();
          
          double flip = 1.;
          if( beamDir.Z() < 0. ) flip = -1.;
          beamDirX = flip * beamDir.X(); 
          beamDirY = flip * beamDir.Y(); 
          beamDirZ = flip * beamDir.Z(); 

          if( trackDir.Z() < 0. ) flip = -1.;
          else flip = 1.;
          trackDirX = flip * trackDir.X(); 
          trackDirY = flip * trackDir.Y(); 
          trackDirZ = flip * trackDir.Z(); 

          std::cout << "beamDirX: " << beamDir.X() << std::endl;
          std::cout << "beamDirY: " << beamDir.Y() << std::endl;
          std::cout << "beamDirZ: " << beamDir.Z() << std::endl;
          std::cout << "trackDirX: " << trackDir.X() << std::endl;
          std::cout << "trackDirY: " << trackDir.Y() << std::endl;
          std::cout << "trackDirZ: " << trackDir.Z() << std::endl;

          beam_costheta = beamDirX*trackDirX + beamDirY*trackDirY + beamDirZ*trackDirZ;
        }

        if( newTracks.size() == 1 ){
          auto newDir = newTracks.at(0).StartDirection();

          double flip = 1.;
          if( newDir.Z() < 0. ) flip = -1.;

          newDirX = flip * newDir.X(); 
          newDirY = flip * newDir.Y(); 
          newDirZ = flip * newDir.Z(); 

          std::cout << "newDirX: " << newDir.X() << std::endl;
          std::cout << "newDirY: " << newDir.Y() << std::endl;
          std::cout << "newDirZ: " << newDir.Z() << std::endl;
          new_beam_costheta = newDirX*trackDirX + newDirY*trackDirY + newDirZ*trackDirZ;
        }
      }
      else{

        art::ServiceHandle<cheat::BackTrackerService> bt_serv;

        MC = 1;
        // Get the truth utility to help us out
        protoana::ProtoDUNETruthUtils truthUtil;
        // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
        // simulation has fGeneratorTag = "generator"
        auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
        // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
        // of these, so we pass the first element into the function to get the good particle
        const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
        if(geantGoodParticle != 0x0){
          std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() << std::endl;
          std::cout << "ID: " << geantGoodParticle->TrackId() << std::endl;
          std::cout << "Mother: " << geantGoodParticle->Mother() << std::endl;
        }

        if( geantGoodParticle->PdgCode() != 211 ) return;
        std::cout << "Found Pi+" << std::endl;
        
        std::cout << "Has " << geantGoodParticle->NumberDaughters() << " daughters" << std::endl;

        art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
        const sim::ParticleList & plist = pi_serv->ParticleList(); 


        for( int i = 0; i < geantGoodParticle->NumberDaughters(); ++i ){
          int daughterID = geantGoodParticle->Daughter(i);

          //Skip photons, neutrons, the nucleus
          if( plist[ daughterID ]->PdgCode() == 22 || plist[ daughterID ]->PdgCode() == 2112 || plist[ daughterID ]->PdgCode() > 1000000000  ) continue;

          std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
          auto part = plist[ daughterID ];
          int pid = part->PdgCode();
          std::cout << "PID: " << pid << std::endl;
          std::cout << "Start: " << part->Position(0).X() << " " << part->Position(0).Y() << " " << part->Position(0).Z() << std::endl;
          std::cout << "End: " << part->EndPosition().X() << " " << part->EndPosition().Y() << " " << part->EndPosition().Z() << std::endl;
          std::cout << "Len: " << part->Trajectory().TotalLength() << std::endl;

          if( pid == 211 )  nPiPlus_truth++;
          if( pid == -211 ) nPiMinus_truth++;
          if( pid == 111 )  nPi0_truth++;
          if( pid == 2212 ) nProton_truth++;

        }

        //Want to see all of the truth particles that contributed to this track
        std::vector< std::pair< const simb::MCParticle*, double > > contribParts = truthUtil.GetAllMCParticlesFromRecoTrack(*thisTrack, evt, fTrackerTag);
        std::cout << contribParts.size() << " Truth Particles Contributed to this track" << std::endl;
        for( size_t ip = 0; ip < contribParts.size(); ++ip ){
          auto part = contribParts.at( ip ).first;
          double energy = contribParts.at( ip ).second;
          std::cout << ip << " " << part->TrackId() << " " << part->PdgCode() << " " << energy << std::endl;
        }
      }

      // Now we can look for the interaction point of the particle if one exists, i.e where the particle
      // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
      // check this by making sure we get a sensible position
      const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
      std::cout << "Interaction Vertex: " << interactionVtx.X() << " " <<  interactionVtx.Y() << " " <<  interactionVtx.Z() << std::endl;

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
        dQdX.push_back( calo_dQdX[i] );
        dEdX.push_back( calo_dEdX[i] );
        resRange.push_back( calo_range[i] );
      }


      //Go through the track-like daughters and save their calorimetry
      for( size_t i = 0; i < trackDaughters.size(); ++i ){
        auto daughterTrack = trackDaughters.at(i);
        std::vector< anab::Calorimetry > dummy_calo = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
        auto dummy_dQdx = dummy_calo[0].dQdx();
        auto dummy_dEdx = dummy_calo[0].dEdx();
        auto dummy_Range = dummy_calo[0].ResidualRange();
 
        daughter_dQdX.push_back( std::vector<double>() );   
        daughter_dEdX.push_back( std::vector<double>() );
        daughter_resRange.push_back( std::vector<double>() );

        for( size_t i = 0; i < dummy_dQdx.size(); ++i ){
          daughter_dQdX.back().push_back( dummy_dQdx[i] );
          daughter_dEdX.back().push_back( dummy_dEdx[i] );
          daughter_resRange.back().push_back( dummy_Range[i] );
        }
      }

      //Looking for tracks ending at APA3-2 boundary
      //
      if( fBrokenTrackZ_low < endZ && endZ < fBrokenTrackZ_high ){
        std::cout << "Possibly broken track " << std::endl; 
        std::cout << "Track: " <<  beamTrackID << " at " << endZ << std::endl;
        const auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackModuleLabel);

        //Look through all of the tracks
//        for( size_t i = 0; i < recoTracks->size(); ++i ){
          
//          recob::Track tr = recoTracks[i];
        for( auto const & tr : *recoTracks ){          

          //Skip the beam particle
          if( tr.ID() == beamTrackID ) continue;

          std::cout << "Checking track: " << tr.ID() << std::endl;
          
          //Check if the track is close enough to the APA boundary
          double stitchStartZ = tr.Trajectory().Start().Z();
          std::cout << "Start: " << stitchStartZ << std::endl;
          if( fStitchTrackZ_low < stitchStartZ && stitchStartZ < fStitchTrackZ_high ){

            double deltaX = fabs(endX - tr.Trajectory().Start().X());
            double deltaY = fabs(endY - tr.Trajectory().Start().Y());
            double deltaZ = tr.Trajectory().Start().Z() - endZ;

            if( deltaX < fStitchXTol && deltaY < fStitchYTol ){
              std::cout << "Found track " << tr.ID() << " within (dx, dy, dz): " << deltaX << " " <<  deltaY << " " << deltaZ 
                        << " of the beam track" << std::endl;

              stitch_len = tr.Length();
              //Possibly add more requirements?
              stitchTrackID.push_back(tr.ID());

              //Check the cosine of the angle between them
              auto stitchDir = tr.StartDirection();
              stitch_cos_theta.push_back( stitchDir.X()*thisTrack->EndDirection().X() + stitchDir.Y()*thisTrack->EndDirection().Y() + stitchDir.Z()*thisTrack->EndDirection().Z() );
              std::cout << "cosine between possible stitch: " << stitch_cos_theta.back() << std::endl;
              

              //Get calorimetry of the stitched track
              std::vector< anab::Calorimetry> stitch_calo = trackUtil.GetRecoTrackCalorimetry(tr, evt, fTrackerTag, fCalorimetryTag);
              auto temp_dQdX  = stitch_calo[0].dQdx();
              auto temp_dEdX  = stitch_calo[0].dEdx();
              auto temp_range = stitch_calo[0].ResidualRange();
              for( size_t i = 0; i < temp_dQdX.size(); ++i ){
                stitch_dQdX.push_back( temp_dQdX[i] );
                stitch_dEdX.push_back( temp_dEdX[i] );
                stitch_resRange.push_back( temp_range[i] );
              }
       
              ///This will possibly be out of order
              combined_resRange = stitch_resRange;
              if( stitch_resRange[0] > stitch_resRange.back() ){      
                for( size_t i = 0; i < resRange.size(); ++i ){
                  combined_resRange.push_back( resRange[i] + stitch_resRange[0] );
                }
              }
              else{
                for( size_t i = 0; i < resRange.size(); ++i ){
                  combined_resRange.push_back( resRange[i] + stitch_resRange.back() );
                }
              }
              combined_dQdX = stitch_dQdX;
              combined_dQdX.insert(combined_dQdX.end(), dQdX.begin(), dQdX.end() );
              combined_dEdX = stitch_dEdX;
              combined_dEdX.insert(combined_dEdX.end(), dEdX.begin(), dEdX.end() );

              combined_len = stitch_len + len;

            }
          }
        }
      }
    }  
  }
  

  fTree->Fill();
}

void pionana::PionAnalyzer::beginJob()
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
  fTree->Branch("stitch_len", &stitch_len);
  fTree->Branch("combined_len", &combined_len);
  fTree->Branch("run", &run);
  fTree->Branch("event", &event);
  fTree->Branch("type", &type);

  fTree->Branch("beam_costheta", &beam_costheta);
  fTree->Branch("new_beam_costheta", &new_beam_costheta);
  fTree->Branch("beamDirX", &beamDirX);
  fTree->Branch("beamDirY", &beamDirY);
  fTree->Branch("beamDirZ", &beamDirZ);
  fTree->Branch("trackDirX", &trackDirX);
  fTree->Branch("trackDirY", &trackDirY);
  fTree->Branch("trackDirZ", &trackDirZ);
  fTree->Branch("newDirX", &newDirX);
  fTree->Branch("newDirY", &newDirY);
  fTree->Branch("newDirZ", &newDirZ);
  fTree->Branch("beamTrackID", &beamTrackID);
  fTree->Branch("stitchTrackID", &stitchTrackID);
  fTree->Branch("stitch_cos_theta", &stitch_cos_theta);
  fTree->Branch("stitch_dQdX", &stitch_dQdX);
  fTree->Branch("stitch_dEdX", &stitch_dEdX);
  fTree->Branch("stitch_resRange", &stitch_resRange);
  fTree->Branch("combined_dQdX", &combined_dQdX);
  fTree->Branch("combined_dEdX", &combined_dEdX);
  fTree->Branch("combined_resRange", &combined_resRange);

  fTree->Branch("MC", &MC);
  fTree->Branch("dQdX", &dQdX);
  fTree->Branch("dEdX", &dEdX);
  fTree->Branch("resRange", &resRange);
  fTree->Branch("daughter_dQdX", &daughter_dQdX);
  fTree->Branch("daughter_dEdX", &daughter_dEdX);
  fTree->Branch("daughter_resRange", &daughter_resRange);
  fTree->Branch("daughter_len", &daughter_len);
  fTree->Branch("nProton_truth", &nProton_truth);
  fTree->Branch("nPi0_truth", &nPi0_truth);
  fTree->Branch("nPiPlus_truth", &nPiPlus_truth);
  fTree->Branch("nPiMinus_truth", &nPiMinus_truth);

  
}

void pionana::PionAnalyzer::endJob()
{

}

void pionana::PionAnalyzer::reset()
{
  startX = -1;
  startY = -1;
  startZ = -1;
  endX = -1;
  endY = -1;
  endZ = -1;

  len = -1;
  stitch_len = -1;
  combined_len = -1;
  type = -1;
  beam_costheta = -100;
  new_beam_costheta = -100;

  MC = 0;
  nProton_truth = 0;
  nPi0_truth = 0;
  nPiPlus_truth = 0;
  nPiMinus_truth = 0;

  dQdX.clear();
  dEdX.clear();
  resRange.clear();

  daughter_dQdX.clear();
  daughter_dEdX.clear();
  daughter_resRange.clear();
  daughter_len.clear();

  beamTrackID = -1;

  if( stitchTrackID.size() ){
    stitchTrackID.clear();
    stitch_cos_theta.clear();

    stitch_dQdX.clear();
    stitch_dEdX.clear();
    stitch_resRange.clear();
    combined_dQdX.clear();
    combined_dEdX.clear();
    combined_resRange.clear();
  }

}

DEFINE_ART_MODULE(pionana::PionAnalyzer)
