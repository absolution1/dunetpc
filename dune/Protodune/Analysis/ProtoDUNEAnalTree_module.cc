////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEAnalTree
// File:        ProtoDUNEAnalTree_module.cc
//
// Extract protoDUNE useful information, do a firs tpre-selection and save output to a flat tree
// 
// Some parts are copied from the beam module example
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
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
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEShowerUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"

// C++ Includes
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// Maximum number of beam particles to save
const int NMAXBEAMPARTICLES = 10;
const int NMAXDAUGTHERS = 20;
const int NMAXGRANDDAUGTHERS = 20;
const int NMAXT0S = 10;

namespace protoana {
  class ProtoDUNEAnalTree;
}


class protoana::ProtoDUNEAnalTree : public art::EDAnalyzer {
public:

  explicit ProtoDUNEAnalTree(fhicl::ParameterSet const & p);

  ProtoDUNEAnalTree(ProtoDUNEAnalTree const &) = delete;
  ProtoDUNEAnalTree(ProtoDUNEAnalTree &&) = delete;
  ProtoDUNEAnalTree & operator = (ProtoDUNEAnalTree const &) = delete;
  ProtoDUNEAnalTree & operator = (ProtoDUNEAnalTree &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & evt) override;

private:

  // Initialise tree variables
  void Initialise();

  // Fill cosmics tree
  void FillCosmicsTree(art::Event const & evt, std::string pfParticleTag);

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  bool fVerbose;

  TTree *fPandoraBeam;
  TTree *fPandoraCosmics;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
  int fNactivefembs[6];

  // Beam tracks
  int fNBEAMPARTICLES;
  int fbeamtrigger[NMAXBEAMPARTICLES];
  double ftof[NMAXBEAMPARTICLES];
  int fcerenkov1[NMAXBEAMPARTICLES];
  int fcerenkov2[NMAXBEAMPARTICLES];
  double fbeamtrackMomentum[NMAXBEAMPARTICLES];
  double fbeamtrackPos[NMAXBEAMPARTICLES][3];
  double fbeamtrackDir[NMAXBEAMPARTICLES][3];

  // Reconstructed tracks/showers
  int fNPRIMARYPARTICLES;
  double fvertex[3];
  double fsecvertex[3];

  int fisprimarytrack;
  int fisprimaryshower;
  double fprimaryBDTScore;
  int fprimaryNHits;
  double fprimaryTheta;
  double fprimaryPhi;
  double fprimaryLength;
  double fprimaryMomentum;
  double fprimaryEndMomentum;
  double fprimaryEndPosition[3];
  double fprimaryStartPosition[3];
  double fprimaryEndDirection[3];
  double fprimaryStartDirection[3];
  double fprimaryOpeningAngle;
  int fprimaryShowerBestPlane;
  double fprimaryShowerEnergy;
  double fprimaryShowerMIPEnergy;
  double fprimaryShowerdEdx;
  double fprimaryMomentumByRangeProton;
  double fprimaryMomentumByRangeMuon;
  double fprimaryKineticEnergy[3];
  double fprimaryRange[3];
  int fprimaryID;
  int fNPRIMARYT0S;
  double fprimaryT0s[NMAXT0S];

  // Daughters from primary
  int fNDAUGHTERS;
  int fisdaughtertrack[NMAXDAUGTHERS];
  int fisdaughtershower[NMAXDAUGTHERS];
  int fdaughterNHits[NMAXDAUGTHERS];
  double fdaughterTheta[NMAXDAUGTHERS];
  double fdaughterPhi[NMAXDAUGTHERS];
  double fdaughterLength[NMAXDAUGTHERS];
  double fdaughterMomentum[NMAXDAUGTHERS];
  double fdaughterEndMomentum[NMAXDAUGTHERS];
  double fdaughterEndPosition[NMAXDAUGTHERS][3];
  double fdaughterStartPosition[NMAXDAUGTHERS][3];
  double fdaughterEndDirection[NMAXDAUGTHERS][3];
  double fdaughterStartDirection[NMAXDAUGTHERS][3];
  double fdaughterOpeningAngle[NMAXDAUGTHERS];
  double fdaughterShowerEnergy[NMAXDAUGTHERS];
  double fdaughterShowerMIPEnergy[NMAXDAUGTHERS];
  double fdaughterShowerdEdx[NMAXDAUGTHERS];
  int fdaughterShowerBestPlane[NMAXDAUGTHERS];
  double fdaughterMomentumByRangeProton[NMAXDAUGTHERS];
  double fdaughterMomentumByRangeMuon[NMAXDAUGTHERS];
  double fdaughterKineticEnergy[NMAXDAUGTHERS][3];
  double fdaughterRange[NMAXDAUGTHERS][3];
  int fdaughterID[NMAXDAUGTHERS];

};


protoana::ProtoDUNEAnalTree::ProtoDUNEAnalTree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose"))
{

}

void protoana::ProtoDUNEAnalTree::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fPandoraBeam = tfs->make<TTree>("PandoraBeam", "Beam events reconstructed with Pandora");
  fPandoraBeam->Branch("run",                           &fRun,                          "run/I");
  fPandoraBeam->Branch("subrun",                        &fSubRun,                       "subrun/I");
  fPandoraBeam->Branch("event",                         &fevent,                        "event/I");
  fPandoraBeam->Branch("timestamp",                     &fTimeStamp,                    "timestamp/D");
  fPandoraBeam->Branch("Nactivefembs",                  &fNactivefembs,                 "Nactivefembs[5]/I");

  fPandoraBeam->Branch("NBEAMPARTICLES",                &fNBEAMPARTICLES,               "NBEAMPARTICLES/I");
  fPandoraBeam->Branch("beamtrigger",                   &fbeamtrigger,                  "beamtrigger[NBEAMPARTICLES]/I");
  fPandoraBeam->Branch("tof",                           &ftof,                          "tof[NBEAMPARTICLES]/D");
  fPandoraBeam->Branch("cerenkov1",                     &fcerenkov1,                    "cerenkov1[NBEAMPARTICLES]/I");
  fPandoraBeam->Branch("cerenkov2",                     &fcerenkov2,                    "cerenkov2[NBEAMPARTICLES]/I");
  fPandoraBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum[NBEAMPARTICLES]/D");
  fPandoraBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[NBEAMPARTICLES][3]/D");
  fPandoraBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[NBEAMPARTICLES][3]/D");

  fPandoraBeam->Branch("NPRIMARYPARTICLES",             &fNPRIMARYPARTICLES,            "NPRIMARYPARTICLES/I");
  fPandoraBeam->Branch("vertex",                        &fvertex,                       "vertex[3]/D");
  fPandoraBeam->Branch("secvertex",                     &fsecvertex,                    "secvertex[3]/D");
  fPandoraBeam->Branch("isprimarytrack",                &fisprimarytrack,               "isprimarytrack/I");
  fPandoraBeam->Branch("isprimaryshower",               &fisprimaryshower,              "isprimaryshower/I");
  fPandoraBeam->Branch("primaryBDTScore",               &fprimaryBDTScore,              "primaryBDTScore/D");
  fPandoraBeam->Branch("primaryNHits",                  &fprimaryNHits,                 "fprimaryNHits/I");
  fPandoraBeam->Branch("primaryTheta",                  &fprimaryTheta,                 "primaryTheta/D");
  fPandoraBeam->Branch("primaryPhi",                    &fprimaryPhi,                   "primaryPhi/D");
  fPandoraBeam->Branch("primaryLength",                 &fprimaryLength,                "primaryLength/D");
  fPandoraBeam->Branch("primaryMomentum",               &fprimaryMomentum,              "primaryMomentum/D");
  fPandoraBeam->Branch("primaryEndMomentum",            &fprimaryEndMomentum,           "primaryEndMomentum/D");
  fPandoraBeam->Branch("primaryEndPosition",            &fprimaryEndPosition,           "primaryEndPosition[3]/D");
  fPandoraBeam->Branch("primaryStartPosition",          &fprimaryStartPosition,         "primaryStartPosition[3]/D");
  fPandoraBeam->Branch("primaryEndDirection",           &fprimaryEndDirection,          "primaryEndDirection[3]/D");
  fPandoraBeam->Branch("primaryStartDirection",         &fprimaryStartDirection,        "primaryStartDirection[3]/D");
  fPandoraBeam->Branch("primaryOpeningAngle",           &fprimaryOpeningAngle,          "primaryOpeningAngle/D");
  fPandoraBeam->Branch("primaryID",                     &fprimaryID,                    "primaryID/I");
  fPandoraBeam->Branch("primaryShowerBestPlane",        &fprimaryShowerBestPlane,       "primaryShowerBestPlane/I");
  fPandoraBeam->Branch("primaryShowerEnergy",           &fprimaryShowerEnergy,          "primaryShowerEnergy/I");
  fPandoraBeam->Branch("primaryShowerMIPEnergy",        &fprimaryShowerMIPEnergy,       "primaryShowerMIPEnergy/I");
  fPandoraBeam->Branch("primaryShowerdEdx",             &fprimaryShowerdEdx,            "primaryShowerdEdx/I");
  fPandoraBeam->Branch("primaryMomentumByRangeProton",  &fprimaryMomentumByRangeProton, "primaryMomentumByRangeProton/D");
  fPandoraBeam->Branch("primaryMomentumByRangeMuon",    &fprimaryMomentumByRangeMuon,   "primaryMomentumByRangeMuon/D");
  fPandoraBeam->Branch("primaryKineticEnergy",          &fprimaryKineticEnergy,         "primaryKineticEnergy[3]/D");
  fPandoraBeam->Branch("primaryRange",                  &fprimaryRange,                 "primaryRange[3]/D");
  fPandoraBeam->Branch("NPRIMARYT0S",                   &fNPRIMARYT0S,                  "NPRIMARYT0S/I");
  fPandoraBeam->Branch("primaryT0s",                    &fprimaryT0s,                   "primaryT0s[NPRIMARYT0S]/D");

  fPandoraBeam->Branch("NDAUGHTERS",                    &fNDAUGHTERS,                   "NDAUGHTERS/I");
  fPandoraBeam->Branch("isdaughtertrack",               &fisdaughtertrack,              "isdaughtertrack[NDAUGHTERS]/I");
  fPandoraBeam->Branch("isdaughtershower",              &fisdaughtershower,             "isdaughtershower[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterNHits",                 &fdaughterNHits,                "daughterNHits[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterTheta",                 &fdaughterTheta,                "daughterTheta[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterPhi",                   &fdaughterPhi,                  "daughterPhi[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterLength",                &fdaughterLength,               "daughterLength[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentum",              &fdaughterMomentum,             "daughterMomentum[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndMomentum",           &fdaughterEndMomentum,          "daughterEndMomentum[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndPosition",           &fdaughterEndPosition,          "daughterEndPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartPosition",         &fdaughterStartPosition,        "daughterStartPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartDirection",        &fdaughterStartDirection,       "daughterStartDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterEndDirection",          &fdaughterEndDirection,         "daughterEndDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterOpeningAngle",          &fdaughterOpeningAngle,         "daughterOpeningAngle[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerBestPlane",       &fdaughterShowerBestPlane,      "daughterShowerBestPlane[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerEnergy",          &fdaughterShowerEnergy,         "daughterShowerEnergy[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerMIPEnergy",       &fdaughterShowerMIPEnergy,      "daughterShowerMIPEnergy[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerdEdx",            &fdaughterShowerdEdx,           "daughterShowerdEdx[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentumByRangeProton", &fdaughterMomentumByRangeProton,"daughterMomentumByRangeProton[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentumByRangeMuon",   &fdaughterMomentumByRangeMuon,  "daughterMomentumByRangeMuon[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterKineticEnergy",         &fdaughterKineticEnergy,        "daughterKineticEnergy[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterRange",                 &fdaughterRange,                "daughterRange[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterID",                    &fdaughterID,                   "daughterID[NDAUGHTERS]/I");

  fPandoraCosmics = tfs->make<TTree>("PandoraCosmics", "Cosmic tracks reconstructed with Pandora");
  fPandoraCosmics->Branch("run",                 &fRun,                "run/I");
  fPandoraCosmics->Branch("subrun",              &fSubRun,             "subrun/I");
  fPandoraCosmics->Branch("event",               &fevent,              "event/I");
  fPandoraCosmics->Branch("timestamp",           &fTimeStamp,          "timestamp/D");
  fPandoraCosmics->Branch("Nactivefembs",        &fNactivefembs,       "Nactivefembs[5]/I");
  fPandoraCosmics->Branch("NBEAMPARTICLES",      &fNBEAMPARTICLES,     "NBEAMPARTICLES/I");
  fPandoraCosmics->Branch("beamtrigger",         &fbeamtrigger,        "beamtrigger[NBEAMPARTICLES]/I");
  fPandoraCosmics->Branch("tof",                 &ftof,                "tof[NBEAMPARTICLES]/D");
  fPandoraCosmics->Branch("cerenkov1",           &fcerenkov1,          "cerenkov1[NBEAMPARTICLES]/I");
  fPandoraCosmics->Branch("cerenkov2",           &fcerenkov2,          "cerenkov2[NBEAMPARTICLES]/I");
  fPandoraCosmics->Branch("beamtrackMomentum",   &fbeamtrackMomentum,  "beamtrackMomentum[NBEAMPARTICLES]/D");
  fPandoraCosmics->Branch("beamtrackPos",        &fbeamtrackPos,       "beamtrackPos[NBEAMPARTICLES][3]/D");
  fPandoraCosmics->Branch("beamtrackDir",        &fbeamtrackDir,       "beamtrackDir[NBEAMPARTICLES][3]/D");

}

void protoana::ProtoDUNEAnalTree::analyze(art::Event const & evt){

  // Initialise tree parameters
  Initialise();

  fRun = evt.run();
  fSubRun = evt.subRun();
  fevent  = evt.id().event(); 
  art::Timestamp ts = evt.time();
  if (ts.timeHigh() == 0){
    TTimeStamp ts2(ts.timeLow());
    fTimeStamp = ts2.AsDouble();
  }
  else{
    TTimeStamp ts2(ts.timeHigh(), ts.timeLow());
    fTimeStamp = ts2.AsDouble();
  }

  // Helper utility functions
  protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trackUtil;

  // Get number of active fembs
  if(!evt.isRealData()){
    for(int k=0; k < 6; k++)
      fNactivefembs[0] = 20;
  }
  else{
    for(int k=0; k < 6; k++)
     fNactivefembs[k] = dataUtil.GetNActiveFembsForAPA(evt, k);
  }

  bool beamTriggerEvent = false;
  // If this event is MC then we can check what the true beam particle is
  if(!evt.isRealData()){
    // For MC always have a beam trigger
    beamTriggerEvent = true;

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
    }
  }
  else{
    // For data we can see if this event comes from a beam trigger
    beamTriggerEvent = dataUtil.IsBeamTrigger(evt);
    //if(beamTriggerEvent){
    //std::cout << "This data event has a beam trigger" << std::endl;
    //}

    art::Handle< std::vector<beam::ProtoDUNEBeamEvent> > pdbeamHandle;
    std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
    if(evt.getByLabel(fBeamModuleLabel, pdbeamHandle))
      art::fill_ptr_vector(beaminfo, pdbeamHandle);
    
    fNBEAMPARTICLES = beaminfo.size();
    if(fNBEAMPARTICLES > NMAXBEAMPARTICLES){
      fNBEAMPARTICLES = NMAXBEAMPARTICLES;
      std::cout << "INFO::Large number of beam particles " << beaminfo.size() << ". Only the first " << NMAXBEAMPARTICLES << "beam particles are processed." << std::endl;
    }

    for(int i = 0; i < fNBEAMPARTICLES; ++i){
      fbeamtrigger[i] = beaminfo[i]->GetTimingTrigger();
      //if(!beaminfo[i]->CheckIsMatched()) continue;

      // If ToF is 0-3 there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
      if(beaminfo[i]->GetTOFChan() >= 0)
	ftof[i] =  beaminfo[i]->GetTOF();

      // Get Cerenkov
      if(beaminfo[i]->GetBITrigger() == 1){
	fcerenkov1[i] = beaminfo[i]->GetCKov0Status();
	fcerenkov2[i] = beaminfo[i]->GetCKov1Status();
      }

      // Beam particle could have more than one tracks - for now take the first one, need to do this properly
      auto & tracks = beaminfo[i]->GetBeamTracks();
      if(!tracks.empty()){
	fbeamtrackPos[0][0] = tracks[0].End().X();
	fbeamtrackPos[0][1] = tracks[0].End().Y();
	fbeamtrackPos[0][2] = tracks[0].End().Z();
	fbeamtrackDir[0][0] = tracks[0].StartDirection().X();
	fbeamtrackDir[0][1] = tracks[0].StartDirection().Y();
	fbeamtrackDir[0][2] = tracks[0].StartDirection().Z();
      }
  
      // Beam momentum
      auto & beammom = beaminfo[i]->GetRecoBeamMomenta();
      if(!beammom.empty())
	fbeamtrackMomentum[i] = beammom[0];
 
    }

  }

  /*
  // Now we want to access the output from Pandora. This comes in the form of particle flow objects (recob::PFParticle).
  // The primary PFParticles are those we want to consider and these PFParticles then have a hierarchy of daughters that
  // describe the whole interaction of a given primary particle
  //
  //                     / daughter track
  //                    /
  //  primary track    /   
  //  ---------------- ---- daughter track
  //                   \
  //                   /\-
  //                   /\\-- daughter shower
  //
  // The above primary PFParticle will have links to three daughter particles, two track-like and one shower-like
  */

  // Track momentum algorithm calculates momentum based on track range
  trkf::TrackMomentumCalculator trmom;
  //trmom.SetMinLength(100);

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  std::cout << "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag) << std::endl;

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  fNPRIMARYPARTICLES = pfParticles.size();

  // We can now look at these particles
  for(const recob::PFParticle* particle : pfParticles){

    // Pandora's BDT beam-cosmic score
    fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);
    
    // NHits associated with this pfParticle
    fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();

    // Get the T0 for this pfParticle
    //std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
    //for(const anab::T0 pfT0 : pfT0vec){
    //fprimaryT0s[fNPRIMARYT0S] = pfT0.Time();
    //fNPRIMARYT0S++;
    //}

    //std::cout << "Pdg Code = " << particle->PdgCode() << std::endl;
    // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
    // of this particle might be more helpful. These return null pointers if not track-like / shower-like
    const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
    if(thisTrack != 0x0){
      fisprimarytrack               = 1;

      fprimaryID                    = thisTrack->ParticleId();
      fprimaryTheta                 = thisTrack->Theta();
      fprimaryPhi                   = thisTrack->Phi();
      fprimaryLength                = thisTrack->Length();
      fprimaryMomentum              = thisTrack->StartMomentum();
      fprimaryEndMomentum           = thisTrack->EndMomentum();

      fprimaryEndPosition[0]        = thisTrack->Trajectory().End().X();
      fprimaryEndPosition[1]        = thisTrack->Trajectory().End().Y();
      fprimaryEndPosition[2]        = thisTrack->Trajectory().End().Z();
      fprimaryStartPosition[0]      = thisTrack->Trajectory().Start().X();
      fprimaryStartPosition[1]      = thisTrack->Trajectory().Start().Y();
      fprimaryStartPosition[2]      = thisTrack->Trajectory().Start().Z();
      fprimaryEndDirection[0]       = thisTrack->Trajectory().EndDirection().X();
      fprimaryEndDirection[1]       = thisTrack->Trajectory().EndDirection().Y();
      fprimaryEndDirection[2]       = thisTrack->Trajectory().EndDirection().Z();
      fprimaryStartDirection[0]     = thisTrack->Trajectory().StartDirection().X();
      fprimaryStartDirection[1]     = thisTrack->Trajectory().StartDirection().Y();
      fprimaryStartDirection[2]     = thisTrack->Trajectory().StartDirection().Z();

      fprimaryMomentumByRangeMuon   = trmom.GetTrackMomentum(thisTrack->Length(),13);
      fprimaryMomentumByRangeProton = trmom.GetTrackMomentum(thisTrack->Length(),2212);      
      
      std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
      if(calovector.size() != 3)
	std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;

      for(size_t k = 0; k < calovector.size() && k<3; k++){
	fprimaryKineticEnergy[k] = calovector[k].KineticEnergy();
	fprimaryRange[k] = calovector[k].Range();
	//const std::vector< double > & dedxvec = calovector[k].dEdx();
	//for(size_t l = 0; l<dedxvec.size(); l++){
	//std::cout << "dE/dx = " << dedxvec[l] << " , Z = " << calovector[k].XYZ()[l].Z() << " , R = " << calovector[k].ResidualRange()[l] << std::endl;
      }

      //std::vector<anab::T0> T0vector = trackUtil.GetRecoTrackT0(*thisTrack, evt, fTrackerTag);
      //std::cout << "t0vector size = " << T0vector.size() << std::endl;

      //std::cout << "Track Length = " << thisTrack->Length() << " , theta = " << thisTrack->Theta() << " , phi = " << thisTrack->Phi() << " , mom = " << thisTrack->VertexMomentum() << std::endl;
      //std::cout << "Trajectory Length = " << thisTrack->Trajectory().Length() << " , theta = " << thisTrack->Trajectory().Theta() << " , phi = " << thisTrack->Trajectory().Phi() << " , mom = " << thisTrack->Trajectory().StartMomentum() << " , end = " << thisTrack->Trajectory().End().X() << std::endl;
      //std::cout << "Beam particle is track-like" << std::endl;
    }
    else if(thisShower != 0x0){
      fisprimaryshower              = 1;

      fprimaryID                    = thisShower->ID();
      fprimaryLength                = thisShower->Length();
      fprimaryShowerBestPlane       = thisShower->best_plane();
      fprimaryOpeningAngle          = thisShower->OpenAngle();
      fprimaryStartPosition[0]      = thisShower->ShowerStart().X();
      fprimaryStartPosition[1]      = thisShower->ShowerStart().Y();
      fprimaryStartPosition[2]      = thisShower->ShowerStart().Z();
      fprimaryStartDirection[0]     = thisShower->Direction().X();
      fprimaryStartDirection[1]     = thisShower->Direction().Y();
      fprimaryStartDirection[2]     = thisShower->Direction().Z();
      if( (thisShower->Energy()).size() > 0 )
	fprimaryShowerEnergy = thisShower->Energy()[0]; // thisShower->best_plane()
      if( (thisShower->MIPEnergy()).size() > 0 )
	fprimaryShowerMIPEnergy = thisShower->MIPEnergy()[0];
      if( (thisShower->dEdx()).size() > 0 )
	fprimaryShowerdEdx = thisShower->dEdx()[0];

      //std::cout << "Shower Length = " << thisShower->Length() << " , open angle = " << thisShower->OpenAngle() << " , start X-Y-Z = " << thisShower->ShowerStart().X() << "-" << thisShower->ShowerStart().Y() << "-" << thisShower->ShowerStart().Z() << " , energy size = " << thisShower->dEdx().size() <<  std::endl;
      //std::cout << "Beam particle is shower-like" << std::endl;
    }
    
    // Find the particle vertex. We need the tracker tag here because we need to do a bit of
    // additional work if the PFParticle is track-like to find the vertex. 
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fvertex[0] = vtx.X(); fvertex[1] = vtx.Y(); fvertex[2] = vtx.Z();

    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fsecvertex[0] = interactionVtx.X(); fsecvertex[1] = interactionVtx.Y(); fsecvertex[2] = interactionVtx.Z();

    // Let's get the daughter PFParticles... we can do this simply without the utility
    for(const int daughterID : particle->Daughters()){
      // Daughter ID is the element of the original recoParticle vector
      const recob::PFParticle *daughterParticle      = &(recoParticles->at(daughterID));
      std::cout << "Daughter " << daughterID << " has " << daughterParticle->NumDaughters() << " daughters" << std::endl;
      
      const recob::Track* daughterTrack              = pfpUtil.GetPFParticleTrack(*daughterParticle,evt,fPFParticleTag,fTrackerTag);
      const recob::Shower* daughterShower            = pfpUtil.GetPFParticleShower(*daughterParticle,evt,fPFParticleTag,fShowerTag);
      
      fdaughterID[fNDAUGHTERS]                       = daughterID;
      // NHits associated with this pfParticle
      fdaughterNHits[fNDAUGHTERS]                    = (pfpUtil.GetPFParticleHits(*daughterParticle,evt,fPFParticleTag)).size();

      if(daughterTrack != 0x0){
	fisdaughtertrack[fNDAUGHTERS]                = 1;
	fdaughterTheta[fNDAUGHTERS]                  = daughterTrack->Theta();
	fdaughterPhi[fNDAUGHTERS]                    = daughterTrack->Phi();
	fdaughterLength[fNDAUGHTERS]                 = daughterTrack->Length();
	fdaughterMomentum[fNDAUGHTERS]               = daughterTrack->StartMomentum();
	fdaughterEndMomentum[fNDAUGHTERS]            = daughterTrack->EndMomentum();
	fdaughterStartPosition[fNDAUGHTERS][0]       = daughterTrack->Trajectory().Start().X();
	fdaughterStartPosition[fNDAUGHTERS][1]       = daughterTrack->Trajectory().Start().Y();
	fdaughterStartPosition[fNDAUGHTERS][2]       = daughterTrack->Trajectory().Start().Z();
	fdaughterEndPosition[fNDAUGHTERS][0]         = daughterTrack->Trajectory().End().X();
	fdaughterEndPosition[fNDAUGHTERS][1]         = daughterTrack->Trajectory().End().Y();
	fdaughterEndPosition[fNDAUGHTERS][2]         = daughterTrack->Trajectory().End().Z();
	fdaughterStartDirection[fNDAUGHTERS][0]      = daughterTrack->Trajectory().StartDirection().X();
	fdaughterStartDirection[fNDAUGHTERS][1]      = daughterTrack->Trajectory().StartDirection().Y();
	fdaughterStartDirection[fNDAUGHTERS][2]      = daughterTrack->Trajectory().StartDirection().Z();
	fdaughterEndDirection[fNDAUGHTERS][0]        = daughterTrack->Trajectory().EndDirection().X();
	fdaughterEndDirection[fNDAUGHTERS][1]        = daughterTrack->Trajectory().EndDirection().Y();
	fdaughterEndDirection[fNDAUGHTERS][2]        = daughterTrack->Trajectory().EndDirection().Z();

	fdaughterMomentumByRangeMuon[fNDAUGHTERS]    = trmom.GetTrackMomentum(daughterTrack->Length(),13);
	fdaughterMomentumByRangeProton[fNDAUGHTERS]  = trmom.GetTrackMomentum(daughterTrack->Length(),2212);

	std::vector<anab::Calorimetry> daughtercalovector = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
	if(daughtercalovector.size() != 3)
	  std::cerr << "WARNING::Calorimetry vector size for primary is = " << daughtercalovector.size() << std::endl;

	for(size_t k = 0; k < daughtercalovector.size() && k<3; k++){
	  fdaughterKineticEnergy[fNDAUGHTERS][k] = daughtercalovector[k].KineticEnergy();
	  fdaughterRange[fNDAUGHTERS][k] = daughtercalovector[k].Range();
	}

      }
      else if(daughterShower != 0x0){
	fisdaughtershower[fNDAUGHTERS]               = 1;
	fdaughterLength[fNDAUGHTERS]                 = daughterShower->Length();
	fdaughterShowerBestPlane[fNDAUGHTERS]        = daughterShower->best_plane();
	fdaughterOpeningAngle[fNDAUGHTERS]           = daughterShower->OpenAngle();
	fdaughterStartPosition[fNDAUGHTERS][0]       = daughterShower->ShowerStart().X();
	fdaughterStartPosition[fNDAUGHTERS][1]       = daughterShower->ShowerStart().Y();
	fdaughterStartPosition[fNDAUGHTERS][2]       = daughterShower->ShowerStart().Z();
	fdaughterStartDirection[fNDAUGHTERS][0]      = daughterShower->Direction().X();
	fdaughterStartDirection[fNDAUGHTERS][1]      = daughterShower->Direction().Y();
	fdaughterStartDirection[fNDAUGHTERS][2]      = daughterShower->Direction().Z();
	if( (daughterShower->Energy()).size() > 0 )
	  fdaughterShowerEnergy[fNDAUGHTERS] = daughterShower->Energy()[0]; // thisShower->best_plane()
	if( (daughterShower->MIPEnergy()).size() > 0 )
	  fdaughterShowerMIPEnergy[fNDAUGHTERS] = daughterShower->MIPEnergy()[0];
	if( (daughterShower->dEdx()).size() > 0 )
	  fdaughterShowerdEdx[fNDAUGHTERS] = daughterShower->dEdx()[0];
      }

      fNDAUGHTERS++;

    }
 
    // For actually studying the objects, it is easier to have the daughters in their track and shower forms.
    // We can use the utility to get a vector of track-like and a vector of shower-like daughters
    //const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
    //const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);  
    //fNDAUGHTERS = trackDaughters.size() + showerDaughters.size();
    //for(const recob::Track* DaughterParticleTrack : trackDaughters){
      
    //}
    //for(unsigned int k=0; k < trackDaughters.size(); k++){
    //std::cout << "Daughter start pos X = " << trackDaughters[k]->Trajectory().Start().X() << std::endl;
    //}
    //std::cout << "Beam particle has " << trackDaughters.size() << " track-like daughters and " << showerDaughters.size() << " shower-like daughters." << std::endl;

    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found
    break;
  } 

  // Fill trees
  if(beamTriggerEvent)
    fPandoraBeam->Fill();

  fPandoraCosmics->Fill();

}

void protoana::ProtoDUNEAnalTree::endJob(){

}

void protoana::ProtoDUNEAnalTree::FillCosmicsTree(art::Event const & evt, std::string pfParticleTag){

  // To fill

}

void protoana::ProtoDUNEAnalTree::Initialise(){
  
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  for(int k=0; k < 3; k++){
    fvertex[k] = -999.0;
    fsecvertex[k] = -999.0;
    fprimaryEndPosition[k] = -999.0;
    fprimaryStartPosition[k] = -999.0;
    fprimaryEndDirection[k] = -999.0;
    fprimaryStartDirection[k] = -999.0;
    fprimaryKineticEnergy[k] = -999.0;
    fprimaryRange[k] = -999.0;
  }

  fNBEAMPARTICLES = 0;
  for(int k=0; k < NMAXBEAMPARTICLES; k++){
    fbeamtrigger[k] = -999;
    ftof[k] = -999.0;
    fcerenkov1[k] = -999;
    fcerenkov2[k] = -999;
    fbeamtrackMomentum[k] = -999.0;
    for(int l=0; l < 3; l++){
      fbeamtrackPos[k][l] = -999.0;
      fbeamtrackDir[k][l] = -999.0;
    }
  }

  fNPRIMARYPARTICLES = 0;
  fNPRIMARYT0S = 0;
  fisprimarytrack = 0;
  fisprimaryshower = 0;
  fNDAUGHTERS = 0;

  fprimaryBDTScore = -999.0;
  fprimaryNHits = -999;
  fprimaryTheta = -999.0;
  fprimaryPhi = -999.0;
  fprimaryLength = -999.0;
  fprimaryMomentum = -999.0;
  fprimaryEndMomentum = -999.0;
  fprimaryOpeningAngle = -999.0;
  fprimaryShowerBestPlane = -999;
  fprimaryShowerEnergy = -999;
  fprimaryShowerMIPEnergy = -999;
  fprimaryShowerdEdx = -999;
  fprimaryID = -999;
  fprimaryMomentumByRangeProton = -999.0;
  fprimaryMomentumByRangeMuon = -999.0;
  
  for(int k=0; k < NMAXT0S; k++){
    fprimaryT0s[k] = -999.0;
  }

  for(int k=0; k < NMAXDAUGTHERS; k++){
    fisdaughtertrack[k] = -999;
    fisdaughtershower[k] = -999;
    fdaughterNHits[k] = -999;
    fdaughterTheta[k] = -999.0;
    fdaughterPhi[k] = -999.0;
    fdaughterLength[k] = -999.0;
    fdaughterMomentum[k] = -999.0;
    fdaughterEndMomentum[k] = -999.0;
    for(int l=0; l < 3; l++){
      fdaughterEndPosition[k][l] = -999.0;
      fdaughterStartPosition[k][l] = -999.0;
      fdaughterEndDirection[k][l] = -999.0;
      fdaughterStartDirection[k][l] = -999.0;
      fdaughterKineticEnergy[k][l] = -999.0;
      fdaughterRange[k][l] = -999.0;
    }
    fdaughterOpeningAngle[k] = -999.0;
    fdaughterShowerBestPlane[k] = -999;
    fdaughterShowerEnergy[k] = -999;
    fdaughterShowerMIPEnergy[k] = -999;
    fdaughterShowerdEdx[k] = -999;
    fdaughterMomentumByRangeProton[k] = -999.0;
    fdaughterMomentumByRangeMuon[k] = -999.0;
    fdaughterID[k] = -999;
  }


}

DEFINE_ART_MODULE(protoana::ProtoDUNEAnalTree)

