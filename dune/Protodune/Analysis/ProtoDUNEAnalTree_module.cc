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
//#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
//#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"

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
#include <stdio.h>
#include <stdlib.h> 
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// Maximum number of beam particles to save
const int NMAXDAUGTHERS = 15;

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

  // Helper utility functions
  protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNETruthUtils truthUtil;

  // Track momentum algorithm calculates momentum based on track range
  trkf::TrackMomentumCalculator trmom;

  // Initialise tree variables
  void Initialise();

  void FillPrimaryPFParticle(art::Event const & evt, const recob::PFParticle* particle);
  void FillPrimaryDaughterPFParticle(art::Event const & evt, const recob::PFParticle* daughterParticle, int daughterID);
  void FillPrimaryGrandDaughterPFParticle(art::Event const & evt, const recob::PFParticle* gdaughterParticle, int daughterID, int gdaughterID);

  // Fill cosmics tree
  //void FillCosmicsTree(art::Event const & evt, std::string pfParticleTag);

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fCalorimetryTag;
  std::string fParticleIDTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  bool fVerbose;

  TTree *fPandoraBeam;
  //TTree *fPandoraCosmics;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
  int fNactivefembs[6];

  // Beam track
  int fbeamtrigger;
  double ftof;
  int fcerenkovStatus[2];
  double fcerenkovTime[2];
  double fcerenkovPressure[2];
  double fbeamtrackMomentum;
  double fbeamtrackEnergy;
  double fbeamtrackPos[3];
  double fbeamtrackDir[3];
  double fbeamtrackTime;
  int fbeamtrackPdg;
  int fbeamtrackID;

  // Reconstructed tracks/showers
  double fprimaryVertex[3];
  int fprimaryIstrack;
  int fprimaryIsshower;
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
  double fprimaryTrkPitchC[3];
  int fprimaryID;
  double fprimaryT0;

  int fprimaryPID_Pdg[3];
  int fprimaryPID_Ndf[3];
  double fprimaryPID_MinChi2[3];
  double fprimaryPID_DeltaChi2[3];
  double fprimaryPID_Chi2Proton[3];
  double fprimaryPID_Chi2Kaon[3];
  double fprimaryPID_Chi2Pion[3];
  double fprimaryPID_Chi2Muon[3];
  double fprimaryPID_MissingE[3];
  double fprimaryPID_MissingEavg[3];
  double fprimaryPID_PIDA[3];

  int fprimary_truth_TrackId;
  int fprimary_truth_Pdg;
  int fprimary_truth_Mother;
  double fprimary_truth_StartPosition[4];
  double fprimary_truth_EndPosition[4];
  double fprimary_truth_P;
  double fprimary_truth_Momentum[4];
  double fprimary_truth_EndMomentum[4];
  double fprimary_truth_Pt;
  double fprimary_truth_Mass;
  double fprimary_truth_Theta;
  double fprimary_truth_Phi;
  double fprimary_truth_TotalLength;
  int fprimary_truth_Process;
  int fprimary_truth_EndProcess;
  int fprimary_truth_Isbeammatched;
  int fprimary_truth_NDaughters;
  
  // Daughters from primary
  int fNDAUGHTERS;
  double fdaughterVertex[3];
  int fdaughterIstrack[NMAXDAUGTHERS];
  int fdaughterIsshower[NMAXDAUGTHERS];
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
  double fdaughterTrkPitchC[NMAXDAUGTHERS][3];
  int fdaughterID[NMAXDAUGTHERS];
  double fdaughterT0[NMAXDAUGTHERS];

  int fdaughterPID_Pdg[NMAXDAUGTHERS][3];
  int fdaughterPID_Ndf[NMAXDAUGTHERS][3];
  double fdaughterPID_MinChi2[NMAXDAUGTHERS][3];
  double fdaughterPID_DeltaChi2[NMAXDAUGTHERS][3];
  double fdaughterPID_Chi2Proton[NMAXDAUGTHERS][3];
  double fdaughterPID_Chi2Kaon[NMAXDAUGTHERS][3];
  double fdaughterPID_Chi2Pion[NMAXDAUGTHERS][3];
  double fdaughterPID_Chi2Muon[NMAXDAUGTHERS][3];
  double fdaughterPID_MissingE[NMAXDAUGTHERS][3];
  double fdaughterPID_MissingEavg[NMAXDAUGTHERS][3];
  double fdaughterPID_PIDA[NMAXDAUGTHERS][3];

  int fdaughter_truth_TrackId[NMAXDAUGTHERS];
  int fdaughter_truth_Pdg[NMAXDAUGTHERS];
  int fdaughter_truth_Mother[NMAXDAUGTHERS];
  double fdaughter_truth_StartPosition[NMAXDAUGTHERS][4];
  double fdaughter_truth_EndPosition[NMAXDAUGTHERS][4];
  double fdaughter_truth_P[NMAXDAUGTHERS];
  double fdaughter_truth_Momentum[NMAXDAUGTHERS][4];
  double fdaughter_truth_EndMomentum[NMAXDAUGTHERS][4];
  double fdaughter_truth_Pt[NMAXDAUGTHERS];
  double fdaughter_truth_Mass[NMAXDAUGTHERS];
  double fdaughter_truth_Theta[NMAXDAUGTHERS];
  double fdaughter_truth_Phi[NMAXDAUGTHERS];
  double fdaughter_truth_TotalLength[NMAXDAUGTHERS];
  int fdaughter_truth_Process[NMAXDAUGTHERS];
  int fdaughter_truth_EndProcess[NMAXDAUGTHERS];

  int fNGRANDDAUGHTERS;
  double fgranddaughterVertex[NMAXDAUGTHERS][3];
  int fgranddaughterIstrack[NMAXDAUGTHERS];
  int fgranddaughterIsshower[NMAXDAUGTHERS];
  int fgranddaughterNHits[NMAXDAUGTHERS];
  double fgranddaughterTheta[NMAXDAUGTHERS];
  double fgranddaughterPhi[NMAXDAUGTHERS];
  double fgranddaughterLength[NMAXDAUGTHERS];
  double fgranddaughterMomentum[NMAXDAUGTHERS];
  double fgranddaughterEndMomentum[NMAXDAUGTHERS];
  double fgranddaughterEndPosition[NMAXDAUGTHERS][3];
  double fgranddaughterStartPosition[NMAXDAUGTHERS][3];
  double fgranddaughterEndDirection[NMAXDAUGTHERS][3];
  double fgranddaughterStartDirection[NMAXDAUGTHERS][3];
  double fgranddaughterOpeningAngle[NMAXDAUGTHERS];
  double fgranddaughterShowerEnergy[NMAXDAUGTHERS];
  double fgranddaughterShowerMIPEnergy[NMAXDAUGTHERS];
  double fgranddaughterShowerdEdx[NMAXDAUGTHERS];
  int fgranddaughterShowerBestPlane[NMAXDAUGTHERS];
  double fgranddaughterMomentumByRangeProton[NMAXDAUGTHERS];
  double fgranddaughterMomentumByRangeMuon[NMAXDAUGTHERS];
  double fgranddaughterKineticEnergy[NMAXDAUGTHERS][3];
  double fgranddaughterRange[NMAXDAUGTHERS][3];
  double fgranddaughterTrkPitchC[NMAXDAUGTHERS][3];
  int fgranddaughterID[NMAXDAUGTHERS];
  int fgranddaughterMotherID[NMAXDAUGTHERS];
  double fgranddaughterT0[NMAXDAUGTHERS];

  int fgranddaughterPID_Pdg[NMAXDAUGTHERS][3];
  int fgranddaughterPID_Ndf[NMAXDAUGTHERS][3];
  double fgranddaughterPID_MinChi2[NMAXDAUGTHERS][3];
  double fgranddaughterPID_DeltaChi2[NMAXDAUGTHERS][3];
  double fgranddaughterPID_Chi2Proton[NMAXDAUGTHERS][3];
  double fgranddaughterPID_Chi2Kaon[NMAXDAUGTHERS][3];
  double fgranddaughterPID_Chi2Pion[NMAXDAUGTHERS][3];
  double fgranddaughterPID_Chi2Muon[NMAXDAUGTHERS][3];
  double fgranddaughterPID_MissingE[NMAXDAUGTHERS][3];
  double fgranddaughterPID_MissingEavg[NMAXDAUGTHERS][3];
  double fgranddaughterPID_PIDA[NMAXDAUGTHERS][3];

  int fgranddaughter_truth_TrackId[NMAXDAUGTHERS];
  int fgranddaughter_truth_Pdg[NMAXDAUGTHERS];
  int fgranddaughter_truth_Mother[NMAXDAUGTHERS];
  double fgranddaughter_truth_StartPosition[NMAXDAUGTHERS][4];
  double fgranddaughter_truth_EndPosition[NMAXDAUGTHERS][4];
  double fgranddaughter_truth_P[NMAXDAUGTHERS];
  double fgranddaughter_truth_Momentum[NMAXDAUGTHERS][4];
  double fgranddaughter_truth_EndMomentum[NMAXDAUGTHERS][4];
  double fgranddaughter_truth_Pt[NMAXDAUGTHERS];
  double fgranddaughter_truth_Mass[NMAXDAUGTHERS];
  double fgranddaughter_truth_Theta[NMAXDAUGTHERS];
  double fgranddaughter_truth_Phi[NMAXDAUGTHERS];
  double fgranddaughter_truth_TotalLength[NMAXDAUGTHERS];
  int fgranddaughter_truth_Process[NMAXDAUGTHERS];
  int fgranddaughter_truth_EndProcess[NMAXDAUGTHERS];

};


protoana::ProtoDUNEAnalTree::ProtoDUNEAnalTree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fParticleIDTag(p.get<std::string>("ParticleIDTag")),
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

  fPandoraBeam->Branch("beamtrigger",                   &fbeamtrigger,                  "beamtrigger/I");
  fPandoraBeam->Branch("tof",                           &ftof,                          "tof/D");
  fPandoraBeam->Branch("cerenkovStatus",                &fcerenkovStatus,               "cerenkovStatus[2]/I");
  fPandoraBeam->Branch("cerenkovTime",                  &fcerenkovTime,                 "cerenkovTime[2]/D");
  fPandoraBeam->Branch("cerenkovPressure",              &fcerenkovPressure,             "cerenkovPressure[2]/D");
  fPandoraBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum/D");
  fPandoraBeam->Branch("beamtrackEnergy",               &fbeamtrackEnergy,              "beamtrackEnergy/D");
  fPandoraBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[3]/D");
  fPandoraBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[3]/D");
  fPandoraBeam->Branch("beamtrackTime",                 &fbeamtrackTime,                "beamtrackTime/D");
  fPandoraBeam->Branch("beamtrackPdg",                  &fbeamtrackPdg,                 "beamtrackPdg/I");
  fPandoraBeam->Branch("beamtrackID",                   &fbeamtrackID,                  "beamtrackID/I");

  fPandoraBeam->Branch("primaryVertex",                 &fprimaryVertex,                "primaryVertex[3]/D");
  fPandoraBeam->Branch("primaryIstrack",                &fprimaryIstrack,               "primaryIstrack/I");
  fPandoraBeam->Branch("primaryIsshower",               &fprimaryIsshower,              "primaryIsshower/I");
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
  fPandoraBeam->Branch("primaryTrkPitchC",              &fprimaryTrkPitchC,             "primaryTrkPitchC[3]/D");
  fPandoraBeam->Branch("primaryT0",                     &fprimaryT0,                    "primaryT0/D");
  fPandoraBeam->Branch("primaryPID_Pdg",                &fprimaryPID_Pdg,               "primaryPID_Pdg[3]/I");
  fPandoraBeam->Branch("primaryPID_Ndf",                &fprimaryPID_Ndf,               "primaryPID_Ndf[3]/I");
  fPandoraBeam->Branch("primaryPID_MinChi2",            &fprimaryPID_MinChi2,           "primaryPID_MinChi2[3]/D");
  fPandoraBeam->Branch("primaryPID_DeltaChi2",          &fprimaryPID_DeltaChi2,         "primaryPID_DeltaChi2[3]/D");
  fPandoraBeam->Branch("primaryPID_Chi2Proton",         &fprimaryPID_Chi2Proton,        "primaryPID_Chi2Proton[3]/D");
  fPandoraBeam->Branch("primaryPID_Chi2Kaon",           &fprimaryPID_Chi2Kaon,          "primaryPID_Chi2Kaon[3]/D");
  fPandoraBeam->Branch("primaryPID_Chi2Pion",           &fprimaryPID_Chi2Pion,          "primaryPID_Chi2Pion[3]/D");
  fPandoraBeam->Branch("primaryPID_Chi2Muon",           &fprimaryPID_Chi2Muon,          "primaryPID_Chi2Muon[3]/D");
  fPandoraBeam->Branch("primaryPID_MissingE",           &fprimaryPID_MissingE,          "primaryPID_MissingE[3]/D");
  fPandoraBeam->Branch("primaryPID_MissingEavg",        &fprimaryPID_MissingEavg,       "primaryPID_MissingEavg[3]/D");
  fPandoraBeam->Branch("primaryPID_PIDA",               &fprimaryPID_PIDA,              "primaryPID_PIDA[3]/D");

  fPandoraBeam->Branch("primary_truth_TrackId",         &fprimary_truth_TrackId,         "primary_truth_TrackId/I");
  fPandoraBeam->Branch("primary_truth_Pdg",             &fprimary_truth_Pdg,             "primary_truth_Pdg/I");
  fPandoraBeam->Branch("primary_truth_Mother",          &fprimary_truth_Mother,          "primary_truth_Mother/I");
  fPandoraBeam->Branch("primary_truth_StartPosition",   &fprimary_truth_StartPosition,   "primary_truth_StartPosition[4]/D");
  fPandoraBeam->Branch("primary_truth_EndPosition",     &fprimary_truth_EndPosition,     "primary_truth_EndPosition[4]/D");
  fPandoraBeam->Branch("primary_truth_Momentum",        &fprimary_truth_Momentum,        "primary_truth_Momentum[4]/D");
  fPandoraBeam->Branch("primary_truth_EndMomentum",     &fprimary_truth_EndMomentum,     "primary_truth_EndMomentum[4]/D");
  fPandoraBeam->Branch("primary_truth_P",               &fprimary_truth_P,               "primary_truth_P/D");
  fPandoraBeam->Branch("primary_truth_Pt",              &fprimary_truth_Pt,              "primary_truth_Pt/D");
  fPandoraBeam->Branch("primary_truth_Mass",            &fprimary_truth_Mass,            "primary_truth_Mass/D");
  fPandoraBeam->Branch("primary_truth_Theta",           &fprimary_truth_Theta,           "primary_truth_Theta/D");
  fPandoraBeam->Branch("primary_truth_Phi",             &fprimary_truth_Phi,             "primary_truth_Phi/D");
  fPandoraBeam->Branch("primary_truth_TotalLength",     &fprimary_truth_TotalLength,     "primary_truth_TotalLength/D");
  fPandoraBeam->Branch("primary_truth_Process",         &fprimary_truth_Process,         "primary_truth_Process/I");
  fPandoraBeam->Branch("primary_truth_EndProcess",      &fprimary_truth_EndProcess,      "primary_truth_EndProcess/I");
  fPandoraBeam->Branch("primary_truth_Isbeammatched",   &fprimary_truth_Isbeammatched,   "primary_truth_Isbeammatched/I");
  fPandoraBeam->Branch("primary_truth_NDaughters",      &fprimary_truth_NDaughters,      "primary_truth_NDaughters/I");

  fPandoraBeam->Branch("NDAUGHTERS",                    &fNDAUGHTERS,                    "NDAUGHTERS/I");
  fPandoraBeam->Branch("daughterVertex",                &fdaughterVertex,                "daughterVertex[3]/D");
  fPandoraBeam->Branch("daughterIstrack",               &fdaughterIstrack,               "daughterIstrack[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterIsshower",              &fdaughterIsshower,              "daughterIsshower[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterNHits",                 &fdaughterNHits,                 "daughterNHits[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterTheta",                 &fdaughterTheta,                 "daughterTheta[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterPhi",                   &fdaughterPhi,                   "daughterPhi[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterLength",                &fdaughterLength,                "daughterLength[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentum",              &fdaughterMomentum,              "daughterMomentum[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndMomentum",           &fdaughterEndMomentum,           "daughterEndMomentum[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndPosition",           &fdaughterEndPosition,           "daughterEndPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartPosition",         &fdaughterStartPosition,         "daughterStartPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartDirection",        &fdaughterStartDirection,        "daughterStartDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterEndDirection",          &fdaughterEndDirection,          "daughterEndDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterOpeningAngle",          &fdaughterOpeningAngle,          "daughterOpeningAngle[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerBestPlane",       &fdaughterShowerBestPlane,       "daughterShowerBestPlane[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerEnergy",          &fdaughterShowerEnergy,          "daughterShowerEnergy[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerMIPEnergy",       &fdaughterShowerMIPEnergy,       "daughterShowerMIPEnergy[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerdEdx",            &fdaughterShowerdEdx,            "daughterShowerdEdx[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentumByRangeProton", &fdaughterMomentumByRangeProton, "daughterMomentumByRangeProton[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentumByRangeMuon",   &fdaughterMomentumByRangeMuon,   "daughterMomentumByRangeMuon[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterKineticEnergy",         &fdaughterKineticEnergy,         "daughterKineticEnergy[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterRange",                 &fdaughterRange,                 "daughterRange[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterTrkPitchC",             &fdaughterTrkPitchC,             "daughterTrkPitchC[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterID",                    &fdaughterID,                    "daughterID[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterT0",                    &fdaughterT0,                    "daughterT0[NDAUGHTERS]/D");

  fPandoraBeam->Branch("daughterPID_Pdg",               &fdaughterPID_Pdg,               "daughterPID_Pdg[NDAUGHTERS][3]/I");
  fPandoraBeam->Branch("daughterPID_Ndf",               &fdaughterPID_Ndf,               "daughterPID_Ndf[NDAUGHTERS][3]/I");
  fPandoraBeam->Branch("daughterPID_MinChi2",           &fdaughterPID_MinChi2,           "daughterPID_MinChi2[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterPID_DeltaChi2",         &fdaughterPID_DeltaChi2,         "daughterPID_DeltaChi2[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterPID_Chi2Proton",        &fdaughterPID_Chi2Proton,        "daughterPID_Chi2Proton[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterPID_Chi2Kaon",          &fdaughterPID_Chi2Kaon,          "daughterPID_Chi2Kaon[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterPID_Chi2Pion",          &fdaughterPID_Chi2Pion,          "daughterPID_Chi2Pion[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterPID_Chi2Muon",          &fdaughterPID_Chi2Muon,          "daughterPID_Chi2Muon[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterPID_MissingE",          &fdaughterPID_MissingE,          "daughterPID_MissingE[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterPID_MissingEavg",       &fdaughterPID_MissingEavg,       "daughterPID_MissingEavg[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterPID_PIDA",              &fdaughterPID_PIDA,              "daughterPID_PIDA[NDAUGHTERS][3]/D");

  fPandoraBeam->Branch("daughter_truth_TrackId",        &fdaughter_truth_TrackId,        "daughter_truth_TrackId[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughter_truth_Pdg",            &fdaughter_truth_Pdg,            "daughter_truth_Pdg[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughter_truth_Mother",         &fdaughter_truth_Mother,         "daughter_truth_Mother[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughter_truth_StartPosition",  &fdaughter_truth_StartPosition,  "daughter_truth_StartPosition[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_EndPosition",    &fdaughter_truth_EndPosition,    "daughter_truth_EndPosition[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_Momentum",       &fdaughter_truth_Momentum,       "daughter_truth_Momentum[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_EndMomentum",    &fdaughter_truth_EndMomentum,    "daughter_truth_EndMomentum[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_P",              &fdaughter_truth_P,              "daughter_truth_P[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Pt",             &fdaughter_truth_Pt,             "daughter_truth_Pt[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Mass",           &fdaughter_truth_Mass,           "daughter_truth_Mass[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Theta",          &fdaughter_truth_Theta,          "daughter_truth_Theta[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Phi",            &fdaughter_truth_Phi,            "daughter_truth_Phi[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_TotalLength",    &fdaughter_truth_TotalLength,    "daughter_truth_TotalLength[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Process",        &fdaughter_truth_Process,        "daughter_truth_Process[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughter_truth_EndProcess",     &fdaughter_truth_EndProcess,     "daughter_truth_EndProcess[NDAUGHTERS]/I");

  fPandoraBeam->Branch("granddaughterVertex",           &fgranddaughterVertex,           "granddaughterVertex[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("NGRANDDAUGHTERS",               &fNGRANDDAUGHTERS,               "NGRANDDAUGHTERS/I");
  fPandoraBeam->Branch("granddaughterIstrack",               &fgranddaughterIstrack,               "granddaughterIstrack[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughterIsshower",              &fgranddaughterIsshower,              "granddaughterIsshower[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughterNHits",                 &fgranddaughterNHits,                 "granddaughterNHits[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughterTheta",                 &fgranddaughterTheta,                 "granddaughterTheta[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterPhi",                   &fgranddaughterPhi,                   "granddaughterPhi[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterLength",                &fgranddaughterLength,                "granddaughterLength[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterMomentum",              &fgranddaughterMomentum,              "granddaughterMomentum[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterEndMomentum",           &fgranddaughterEndMomentum,           "granddaughterEndMomentum[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterEndPosition",           &fgranddaughterEndPosition,           "granddaughterEndPosition[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterStartPosition",         &fgranddaughterStartPosition,         "granddaughterStartPosition[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterStartDirection",        &fgranddaughterStartDirection,        "granddaughterStartDirection[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterEndDirection",          &fgranddaughterEndDirection,          "granddaughterEndDirection[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterOpeningAngle",          &fgranddaughterOpeningAngle,          "granddaughterOpeningAngle[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterShowerBestPlane",       &fgranddaughterShowerBestPlane,       "granddaughterShowerBestPlane[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterShowerEnergy",          &fgranddaughterShowerEnergy,          "granddaughterShowerEnergy[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterShowerMIPEnergy",       &fgranddaughterShowerMIPEnergy,       "granddaughterShowerMIPEnergy[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterShowerdEdx",            &fgranddaughterShowerdEdx,            "granddaughterShowerdEdx[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterMomentumByRangeProton", &fgranddaughterMomentumByRangeProton, "granddaughterMomentumByRangeProton[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterMomentumByRangeMuon",   &fgranddaughterMomentumByRangeMuon,   "granddaughterMomentumByRangeMuon[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughterKineticEnergy",         &fgranddaughterKineticEnergy,         "granddaughterKineticEnergy[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterRange",                 &fgranddaughterRange,                 "granddaughterRange[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterTrkPitchC",             &fgranddaughterTrkPitchC,             "granddaughterTrkPitchC[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterID",                    &fgranddaughterID,                    "granddaughterID[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughterMotherID",              &fgranddaughterMotherID,              "granddaughterMotherID[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughterT0",                    &fgranddaughterT0,                    "granddaughterT0[NGRANDDAUGHTERS]/D");

  fPandoraBeam->Branch("granddaughterPID_Pdg",               &fgranddaughterPID_Pdg,               "granddaughterPID_Pdg[NGRANDDAUGHTERS][3]/I");
  fPandoraBeam->Branch("granddaughterPID_Ndf",               &fgranddaughterPID_Ndf,               "granddaughterPID_Ndf[NGRANDDAUGHTERS][3]/I");
  fPandoraBeam->Branch("granddaughterPID_MinChi2",           &fgranddaughterPID_MinChi2,           "granddaughterPID_MinChi2[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterPID_DeltaChi2",         &fgranddaughterPID_DeltaChi2,         "granddaughterPID_DeltaChi2[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterPID_Chi2Proton",        &fgranddaughterPID_Chi2Proton,        "granddaughterPID_Chi2Proton[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterPID_Chi2Kaon",          &fgranddaughterPID_Chi2Kaon,          "granddaughterPID_Chi2Kaon[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterPID_Chi2Pion",          &fgranddaughterPID_Chi2Pion,          "granddaughterPID_Chi2Pion[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterPID_Chi2Muon",          &fgranddaughterPID_Chi2Muon,          "granddaughterPID_Chi2Muon[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterPID_MissingE",          &fgranddaughterPID_MissingE,          "granddaughterPID_MissingE[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterPID_MissingEavg",       &fgranddaughterPID_MissingEavg,       "granddaughterPID_MissingEavg[NGRANDDAUGHTERS][3]/D");
  fPandoraBeam->Branch("granddaughterPID_PIDA",              &fgranddaughterPID_PIDA,              "granddaughterPID_PIDA[NGRANDDAUGHTERS][3]/D");

  fPandoraBeam->Branch("granddaughter_truth_TrackId",        &fgranddaughter_truth_TrackId,        "granddaughter_truth_TrackId[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughter_truth_Pdg",            &fgranddaughter_truth_Pdg,            "granddaughter_truth_Pdg[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughter_truth_Mother",         &fgranddaughter_truth_Mother,         "granddaughter_truth_Mother[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughter_truth_StartPosition",  &fgranddaughter_truth_StartPosition,  "granddaughter_truth_StartPosition[NGRANDDAUGHTERS][4]/D");
  fPandoraBeam->Branch("granddaughter_truth_EndPosition",    &fgranddaughter_truth_EndPosition,    "granddaughter_truth_EndPosition[NGRANDDAUGHTERS][4]/D");
  fPandoraBeam->Branch("granddaughter_truth_Momentum",       &fgranddaughter_truth_Momentum,       "granddaughter_truth_Momentum[NGRANDDAUGHTERS][4]/D");
  fPandoraBeam->Branch("granddaughter_truth_EndMomentum",    &fgranddaughter_truth_EndMomentum,    "granddaughter_truth_EndMomentum[NGRANDDAUGHTERS][4]/D");
  fPandoraBeam->Branch("granddaughter_truth_P",              &fgranddaughter_truth_P,              "granddaughter_truth_P[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughter_truth_Pt",             &fgranddaughter_truth_Pt,             "granddaughter_truth_Pt[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughter_truth_Mass",           &fgranddaughter_truth_Mass,           "granddaughter_truth_Mass[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughter_truth_Theta",          &fgranddaughter_truth_Theta,          "granddaughter_truth_Theta[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughter_truth_Phi",            &fgranddaughter_truth_Phi,            "granddaughter_truth_Phi[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughter_truth_TotalLength",    &fgranddaughter_truth_TotalLength,    "granddaughter_truth_TotalLength[NGRANDDAUGHTERS]/D");
  fPandoraBeam->Branch("granddaughter_truth_Process",        &fgranddaughter_truth_Process,        "granddaughter_truth_Process[NGRANDDAUGHTERS]/I");
  fPandoraBeam->Branch("granddaughter_truth_EndProcess",     &fgranddaughter_truth_EndProcess,     "granddaughter_truth_EndProcess[NGRANDDAUGHTERS]/I");


  //fPandoraCosmics = tfs->make<TTree>("PandoraCosmics", "Cosmic tracks reconstructed with Pandora");
  //fPandoraCosmics->Branch("run",                 &fRun,                "run/I");
  //fPandoraCosmics->Branch("subrun",              &fSubRun,             "subrun/I");
  //fPandoraCosmics->Branch("event",               &fevent,              "event/I");
  //fPandoraCosmics->Branch("timestamp",           &fTimeStamp,          "timestamp/D");
  //fPandoraCosmics->Branch("Nactivefembs",        &fNactivefembs,       "Nactivefembs[5]/I");
  //fPandoraCosmics->Branch("beamtrigger",         &fbeamtrigger,        "beamtrigger/I");
  //fPandoraCosmics->Branch("tof",                 &ftof,                "tof/D");
  //fPandoraCosmics->Branch("cerenkov1",           &fcerenkov1,          "cerenkov1/I");
  //fPandoraCosmics->Branch("cerenkov2",           &fcerenkov2,          "cerenkov2/I");
  //fPandoraCosmics->Branch("beamtrackMomentum",   &fbeamtrackMomentum,  "beamtrackMomentum/D");
  //fPandoraCosmics->Branch("beamtrackPos",        &fbeamtrackPos,       "beamtrackPos[3]/D");
  //fPandoraCosmics->Branch("beamtrackDir",        &fbeamtrackDir,       "beamtrackDir[3]/D");

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

  // Get number of active fembs
  if(!evt.isRealData()){
    for(int k=0; k < 6; k++)
      fNactivefembs[k] = 20;
  }
  else{
    for(int k=0; k < 6; k++)
     fNactivefembs[k] = dataUtil.GetNActiveFembsForAPA(evt, k);
  }

  bool beamTriggerEvent = false;
  // If this event is MC then we can check what the true beam particle is
  if(!evt.isRealData()){
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    // simulation has fGeneratorTag = "generator"
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);

    // Also get the reconstructed beam information in the MC - TO DO
    //auto beamsim = evt.getValidHandle<std::vector<sim::ProtoDUNEbeamsim> >(fGeneratorTag);
    //const sim::ProtoDUNEbeamsim beamsimobj = (*beamsim)[0];
    //std::cout << beamsimobj.NInstruments() << std::endl;
    //sim::ProtoDUNEBeamInstrument beamsim_tof1 = beamsimobj.GetInstrument("TOF1");
    //sim::ProtoDUNEBeamInstrument beamsim_trig2 = beamsimobj.GetInstrument("TRIG2");
    //std::cout << beamsim_trig2.GetT() - beamsim_tof1.GetT() << " , " << beamsim_trig2.GetSmearedVar1() - beamsim_tof1.GetSmearedVar1() << std::endl;
    //std::cout << beamsimobj.GetInstrument("TRIG2").GetT() - beamsimobj.GetInstrument("TOF1").GetT() << " , " << beamsimobj.GetInstrument("TRIG2").GetSmearedVar1() - beamsimobj.GetInstrument("TOF1").GetSmearedVar1() << std::endl;
  
    // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
    // of these, so we pass the first element into the function to get the good particle
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);

    if(geantGoodParticle != 0x0){
      std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() 
		<< " , track id = " << geantGoodParticle->TrackId()
		<< std::endl;

      beamTriggerEvent = true;
      fbeamtrigger       = 12;
      fbeamtrackPos[0]   = geantGoodParticle->Vx();
      fbeamtrackPos[1]   = geantGoodParticle->Vy();
      fbeamtrackPos[2]   = geantGoodParticle->Vz();
      fbeamtrackMomentum = geantGoodParticle->P();
      fbeamtrackEnergy   = geantGoodParticle->E();
      fbeamtrackPdg      = geantGoodParticle->PdgCode();
      fbeamtrackTime     = geantGoodParticle->T();
      fbeamtrackID       = geantGoodParticle->TrackId();
    }
  }
  else{
    // For data we can see if this event comes from a beam trigger
    beamTriggerEvent = dataUtil.IsBeamTrigger(evt);

    art::Handle< std::vector<beam::ProtoDUNEBeamEvent> > pdbeamHandle;
    std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
    if(evt.getByLabel(fBeamModuleLabel, pdbeamHandle))
      art::fill_ptr_vector(beaminfo, pdbeamHandle);
  
    for(unsigned int i = 0; i < beaminfo.size(); ++i){
      //if(!beaminfo[i]->CheckIsMatched()) continue;
      fbeamtrigger = beaminfo[i]->GetTimingTrigger();
      fbeamtrackTime = (double)beaminfo[i]->GetRDTimestamp();

      // If ToF is 0-3 there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
      if(beaminfo[i]->GetTOFChan() >= 0)
	ftof =  beaminfo[i]->GetTOF();

      // Get Cerenkov
      if(beaminfo[i]->GetBITrigger() == 1){
	fcerenkovStatus[0]   = beaminfo[i]->GetCKov0Status();
	fcerenkovStatus[1]   = beaminfo[i]->GetCKov1Status();
	fcerenkovTime[0]     = beaminfo[i]->GetCKov0Time();
	fcerenkovTime[1]     = beaminfo[i]->GetCKov1Time();
	fcerenkovPressure[0] = beaminfo[i]->GetCKov0Pressure();
	fcerenkovPressure[1] = beaminfo[i]->GetCKov1Pressure();
      }

      // Beam particle could have more than one tracks - for now take the first one, need to do this properly
      auto & tracks = beaminfo[i]->GetBeamTracks();
      if(!tracks.empty()){
	fbeamtrackPos[0] = tracks[0].End().X();
	fbeamtrackPos[1] = tracks[0].End().Y();
	fbeamtrackPos[2] = tracks[0].End().Z();
	fbeamtrackDir[0] = tracks[0].StartDirection().X();
	fbeamtrackDir[1] = tracks[0].StartDirection().Y();
	fbeamtrackDir[2] = tracks[0].StartDirection().Z();
      }
  
      // Beam momentum
      auto & beammom = beaminfo[i]->GetRecoBeamMomenta();
      if(!beammom.empty())
	fbeamtrackMomentum = beammom[0];

      // For now only take the first beam particle - need to add some criteria if more than one are found
      break;
 
    }
  }

  /*
  // Now we want to access the output from Pandora. This comes in the form of particle flow objects (recob::PFParticle).
  // The primary PFParticles are those we want to consider and these PFParticles then have a hierarchy of daughters that
  // describe the whole interaction of a given primary particle
  //
  //                     / daughter track    / grand-daughter track
  //                    /                   /
  //  primary track    /                   /
  //  ---------------- ---- daughter track ----- grand-daughter track
  //                   \                   
  //                   /\-
  //                   /\\-- daughter shower
  //
  // The above primary PFParticle will have links to three daughter particles, two track-like and one shower-like
  */

  // Track momentum algorithm calculates momentum based on track range
  //trkf::TrackMomentumCalculator trmom;
  //trmom.SetMinLength(100);

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  //std::cout << "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag) << std::endl;

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);

  // We can now look at these particles
  for(const recob::PFParticle* particle : pfParticles){

    FillPrimaryPFParticle(evt, particle);

    // Find the particle vertex. We need the tracker tag here because we need to do a bit of
    // additional work if the PFParticle is track-like to find the vertex. 
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fprimaryVertex[0] = vtx.X(); fprimaryVertex[1] = vtx.Y(); fprimaryVertex[2] = vtx.Z();

    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fdaughterVertex[0] = interactionVtx.X(); fdaughterVertex[1] = interactionVtx.Y(); fdaughterVertex[2] = interactionVtx.Z();

    // Maximum number of daugthers to be processed
    if(particle->NumDaughters() > NMAXDAUGTHERS)
      std::cout << "INFO::Number of daughters is " << particle->NumDaughters() << ". Only the first NMAXDAUGTHERS are processed." << std::endl;

    // Let's get the daughter PFParticles... we can do this simply without the utility
    for(const int daughterID : particle->Daughters()){
      // Daughter ID is the element of the original recoParticle vector
      const recob::PFParticle *daughterParticle      = &(recoParticles->at(daughterID));
      std::cout << "Daughter " << daughterID << " has " << daughterParticle->NumDaughters() << " daughters" << std::endl;

      // Fill tree with daughter info
      FillPrimaryDaughterPFParticle(evt, daughterParticle, daughterID);

      // Get the secondary vertex from daughter interactions
      const TVector3 secinteractionVtx = pfpUtil.GetPFParticleSecondaryVertex(*daughterParticle,evt,fPFParticleTag,fTrackerTag);
      fgranddaughterVertex[fNDAUGHTERS][0] = secinteractionVtx.X(); fgranddaughterVertex[fNDAUGHTERS][1] = secinteractionVtx.Y(); fgranddaughterVertex[fNDAUGHTERS][2] = secinteractionVtx.Z();

      for(const int gdaughterID : daughterParticle->Daughters()){
	FillPrimaryGrandDaughterPFParticle(evt, daughterParticle, gdaughterID, daughterID);
	
	// Only process NMAXDAUGTHERS
	if(fNGRANDDAUGHTERS > NMAXDAUGTHERS) break;
      }

      // Only process NMAXDAUGTHERS
      if(fNDAUGHTERS > NMAXDAUGTHERS) break;
      
    }

    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
    break;
  } 

  // Fill trees
  if(beamTriggerEvent)
    fPandoraBeam->Fill();

  //fPandoraCosmics->Fill();

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalTree::endJob(){

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalTree::FillPrimaryPFParticle(art::Event const & evt, const recob::PFParticle* particle){

  // Pandora's BDT beam-cosmic score
  fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);
  
  // NHits associated with this pfParticle
  fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();
  
  // Get the T0 for this pfParticle
  std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
  if(!pfT0vec.empty())
    fprimaryT0 = pfT0vec[0].Time();
    
  // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
  // of this particle might be more helpful. These return null pointers if not track-like / shower-like
  const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle, evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);

  if(thisTrack != 0x0){
    fprimaryIstrack                    = 1;
    fprimaryIsshower                   = 0;
    fprimaryID                         = thisTrack->ParticleId();
    fprimaryTheta                      = thisTrack->Theta();
    fprimaryPhi                        = thisTrack->Phi();
    fprimaryLength                     = thisTrack->Length();
    fprimaryMomentum                   = thisTrack->StartMomentum();
    fprimaryEndMomentum                = thisTrack->EndMomentum();
    fprimaryEndPosition[0]             = thisTrack->Trajectory().End().X();
    fprimaryEndPosition[1]             = thisTrack->Trajectory().End().Y();
    fprimaryEndPosition[2]             = thisTrack->Trajectory().End().Z();
    fprimaryStartPosition[0]           = thisTrack->Trajectory().Start().X();
    fprimaryStartPosition[1]           = thisTrack->Trajectory().Start().Y();
    fprimaryStartPosition[2]           = thisTrack->Trajectory().Start().Z();
    fprimaryEndDirection[0]            = thisTrack->Trajectory().EndDirection().X();
    fprimaryEndDirection[1]            = thisTrack->Trajectory().EndDirection().Y();
    fprimaryEndDirection[2]            = thisTrack->Trajectory().EndDirection().Z();
    fprimaryStartDirection[0]          = thisTrack->Trajectory().StartDirection().X();
    fprimaryStartDirection[1]          = thisTrack->Trajectory().StartDirection().Y();
    fprimaryStartDirection[2]          = thisTrack->Trajectory().StartDirection().Z();
    
    fprimaryMomentumByRangeMuon        = trmom.GetTrackMomentum(thisTrack->Length(),13);
    fprimaryMomentumByRangeProton      = trmom.GetTrackMomentum(thisTrack->Length(),2212);      
    
    // Calorimetry
    std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
    if(calovector.size() != 3)
      std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;
    
    for(size_t k = 0; k < calovector.size() && k<3; k++){
      int plane = calovector[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fprimaryKineticEnergy[plane]      = calovector[k].KineticEnergy();
      fprimaryRange[plane]              = calovector[k].Range();
      fprimaryTrkPitchC[plane]          = calovector[k].TrkPitchC();
    }
    
    // PID
    std::vector<anab::ParticleID> pids = trackUtil.GetRecoTrackPID(*thisTrack, evt, fTrackerTag, fParticleIDTag);
    if(pids.size() != 3)
      std::cerr << "WARNING::PID vector size for primary is = " << pids.size() << std::endl;
    
    for(size_t k = 0; k < pids.size() && k<3; k++){
      int plane = pids[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fprimaryPID_Pdg[plane]            = pids[plane].Pdg();
      fprimaryPID_Ndf[plane]            = pids[plane].Ndf();
      fprimaryPID_MinChi2[plane]        = pids[plane].MinChi2();
      fprimaryPID_DeltaChi2[plane]      = pids[plane].DeltaChi2();
      fprimaryPID_Chi2Proton[plane]     = pids[plane].Chi2Proton();
      fprimaryPID_Chi2Kaon[plane]       = pids[plane].Chi2Kaon();
      fprimaryPID_Chi2Pion[plane]       = pids[plane].Chi2Pion();
      fprimaryPID_Chi2Muon[plane]       = pids[plane].Chi2Muon();
      fprimaryPID_MissingE[plane]       = pids[plane].MissingE();
      fprimaryPID_MissingEavg[plane]    = pids[plane].MissingEavg();
      fprimaryPID_PIDA[plane]           = pids[plane].PIDA();
    }
    
    // Get the true mc particle
    const simb::MCParticle* mcparticle = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackerTag);
    if(mcparticle != 0x0){
      fprimary_truth_TrackId            = mcparticle->TrackId();
      fprimary_truth_Pdg                = mcparticle->PdgCode();
      fprimary_truth_Mother             = mcparticle->Mother();
      fprimary_truth_StartPosition[0]   = mcparticle->Vx();
      fprimary_truth_StartPosition[1]   = mcparticle->Vy();
      fprimary_truth_StartPosition[2]   = mcparticle->Vz();
      fprimary_truth_StartPosition[3]   = mcparticle->T();
      fprimary_truth_EndPosition[0]     = mcparticle->EndX();
      fprimary_truth_EndPosition[1]     = mcparticle->EndY();
      fprimary_truth_EndPosition[2]     = mcparticle->EndZ();
      fprimary_truth_EndPosition[3]     = mcparticle->EndT();
      fprimary_truth_P                  = mcparticle->P();
      fprimary_truth_Momentum[0]        = mcparticle->Px();
      fprimary_truth_Momentum[1]        = mcparticle->Py();
      fprimary_truth_Momentum[2]        = mcparticle->Pz();
      fprimary_truth_Momentum[3]        = mcparticle->E();
      fprimary_truth_Pt                 = mcparticle->Pt();
      fprimary_truth_Mass               = mcparticle->Mass();
      fprimary_truth_EndMomentum[0]     = mcparticle->EndPx();
      fprimary_truth_EndMomentum[1]     = mcparticle->EndPy();
      fprimary_truth_EndMomentum[2]     = mcparticle->EndPz();
      fprimary_truth_EndMomentum[3]     = mcparticle->EndE();
      fprimary_truth_Theta              = mcparticle->Momentum().Theta();
      fprimary_truth_Phi                = mcparticle->Momentum().Phi();
      fprimary_truth_NDaughters         = mcparticle->NumberDaughters();
      fprimary_truth_TotalLength        = mcparticle->Trajectory().TotalLength();
      fprimary_truth_Process            = int(mcparticle->Trajectory().ProcessToKey(mcparticle->Process()));
      fprimary_truth_EndProcess         = int(mcparticle->Trajectory().ProcessToKey(mcparticle->EndProcess()));
      if(fbeamtrackID != -999 && fprimary_truth_TrackId == fbeamtrackID)
	fprimary_truth_Isbeammatched    = 1;
      else
	fprimary_truth_Isbeammatched    = 0;
      
      std::cout << "Process = " << (mcparticle->Process()).c_str() << " , End process = " << (mcparticle->EndProcess()).c_str()
		<< " , track ID = " << mcparticle->TrackId()
		<< std::endl;
    }
  } // end is track
  else if(thisShower != 0x0){
    fprimaryIstrack                     = 0;
    fprimaryIsshower                    = 1;
    
    fprimaryID                          = thisShower->ID();
    fprimaryLength                      = thisShower->Length();
    fprimaryShowerBestPlane             = thisShower->best_plane();
    fprimaryOpeningAngle                = thisShower->OpenAngle();
    fprimaryStartPosition[0]            = thisShower->ShowerStart().X();
    fprimaryStartPosition[1]            = thisShower->ShowerStart().Y();
    fprimaryStartPosition[2]            = thisShower->ShowerStart().Z();
    fprimaryStartDirection[0]           = thisShower->Direction().X();
    fprimaryStartDirection[1]           = thisShower->Direction().Y();
    fprimaryStartDirection[2]           = thisShower->Direction().Z();
    if( (thisShower->Energy()).size() > 0 )
      fprimaryShowerEnergy = thisShower->Energy()[0]; // thisShower->best_plane()
    if( (thisShower->MIPEnergy()).size() > 0 )
      fprimaryShowerMIPEnergy = thisShower->MIPEnergy()[0];
    if( (thisShower->dEdx()).size() > 0 )
      fprimaryShowerdEdx = thisShower->dEdx()[0];
  } // end is shower
  else{
    std::cout << "INFO::Primary pfParticle is not track or shower!" << std::endl;
    return;
  }

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalTree::FillPrimaryDaughterPFParticle(art::Event const & evt, const recob::PFParticle* daughterParticle, int daughterID){

  const recob::Track* daughterTrack              = pfpUtil.GetPFParticleTrack(*daughterParticle,evt, fPFParticleTag,fTrackerTag);
  const recob::Shower* daughterShower            = pfpUtil.GetPFParticleShower(*daughterParticle,evt,fPFParticleTag,fShowerTag);
  
  if(daughterTrack != 0x0){
    fdaughterIstrack[fNDAUGHTERS]                = 1;
    fdaughterIsshower[fNDAUGHTERS]               = 0;
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
    
    // Calorimetry
    std::vector<anab::Calorimetry> daughtercalovector = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
    if(daughtercalovector.size() != 3)
      std::cerr << "WARNING::Calorimetry vector size for daughter is = " << daughtercalovector.size() << std::endl;
    
    for(size_t k = 0; k < daughtercalovector.size() && k<3; k++){
      int plane = daughtercalovector[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fdaughterKineticEnergy[fNDAUGHTERS][plane] = daughtercalovector[k].KineticEnergy();
      fdaughterRange[fNDAUGHTERS][plane]         = daughtercalovector[k].Range();
      fdaughterTrkPitchC[fNDAUGHTERS][plane]     = daughtercalovector[k].TrkPitchC();
    }

    // PID
    std::vector<anab::ParticleID> daughterpids = trackUtil.GetRecoTrackPID(*daughterTrack, evt, fTrackerTag, fParticleIDTag);
    if(daughterpids.size() != 3)
      std::cerr << "WARNING::PID vector size for daughter is = " << daughterpids.size() << std::endl;
    
    for(size_t k = 0; k < daughterpids.size() && k<3; k++){
      int plane = daughterpids[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fdaughterPID_Pdg[fNDAUGHTERS][plane]          = daughterpids[plane].Pdg();
      fdaughterPID_Ndf[fNDAUGHTERS][plane]          = daughterpids[plane].Ndf();
      fdaughterPID_MinChi2[fNDAUGHTERS][plane]      = daughterpids[plane].MinChi2();
      fdaughterPID_DeltaChi2[fNDAUGHTERS][plane]    = daughterpids[plane].DeltaChi2();
      fdaughterPID_Chi2Proton[fNDAUGHTERS][plane]   = daughterpids[plane].Chi2Proton();
      fdaughterPID_Chi2Kaon[fNDAUGHTERS][plane]     = daughterpids[plane].Chi2Kaon();
      fdaughterPID_Chi2Pion[fNDAUGHTERS][plane]     = daughterpids[plane].Chi2Pion();
      fdaughterPID_Chi2Muon[fNDAUGHTERS][plane]     = daughterpids[plane].Chi2Muon();
      fdaughterPID_MissingE[fNDAUGHTERS][plane]     = daughterpids[plane].MissingE();
      fdaughterPID_MissingEavg[fNDAUGHTERS][plane]  = daughterpids[plane].MissingEavg();
      fdaughterPID_PIDA[fNDAUGHTERS][plane]         = daughterpids[plane].PIDA();
    }
    
    // Get the true mc particle
    const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(*daughterTrack, evt, fTrackerTag);
    if(mcdaughterparticle != 0x0){
      fdaughter_truth_TrackId[fNDAUGHTERS]          = mcdaughterparticle->TrackId();
      fdaughter_truth_Pdg[fNDAUGHTERS]              = mcdaughterparticle->PdgCode();
      fdaughter_truth_Mother[fNDAUGHTERS]           = mcdaughterparticle->Mother();
      fdaughter_truth_StartPosition[fNDAUGHTERS][0] = mcdaughterparticle->Vx();
      fdaughter_truth_StartPosition[fNDAUGHTERS][1] = mcdaughterparticle->Vy();
      fdaughter_truth_StartPosition[fNDAUGHTERS][2] = mcdaughterparticle->Vz();
      fdaughter_truth_StartPosition[fNDAUGHTERS][3] = mcdaughterparticle->T();
      fdaughter_truth_EndPosition[fNDAUGHTERS][0]   = mcdaughterparticle->EndX();
      fdaughter_truth_EndPosition[fNDAUGHTERS][1]   = mcdaughterparticle->EndY();
      fdaughter_truth_EndPosition[fNDAUGHTERS][2]   = mcdaughterparticle->EndZ();
      fdaughter_truth_EndPosition[fNDAUGHTERS][3]   = mcdaughterparticle->EndT();
      fdaughter_truth_P[fNDAUGHTERS]                = mcdaughterparticle->P();
      fdaughter_truth_Momentum[fNDAUGHTERS][0]      = mcdaughterparticle->Px();
      fdaughter_truth_Momentum[fNDAUGHTERS][1]      = mcdaughterparticle->Py();
      fdaughter_truth_Momentum[fNDAUGHTERS][2]      = mcdaughterparticle->Pz();
      fdaughter_truth_Momentum[fNDAUGHTERS][3]      = mcdaughterparticle->E();
      fdaughter_truth_Pt[fNDAUGHTERS]               = mcdaughterparticle->Pt();
      fdaughter_truth_Mass[fNDAUGHTERS]             = mcdaughterparticle->Mass();
      fdaughter_truth_EndMomentum[fNDAUGHTERS][0]   = mcdaughterparticle->EndPx();
      fdaughter_truth_EndMomentum[fNDAUGHTERS][1]   = mcdaughterparticle->EndPy();
      fdaughter_truth_EndMomentum[fNDAUGHTERS][2]   = mcdaughterparticle->EndPz();
      fdaughter_truth_EndMomentum[fNDAUGHTERS][3]   = mcdaughterparticle->EndE();
      fdaughter_truth_Theta[fNDAUGHTERS]            = mcdaughterparticle->Momentum().Theta();
      fdaughter_truth_Phi[fNDAUGHTERS]              = mcdaughterparticle->Momentum().Phi();
      fdaughter_truth_TotalLength[fNDAUGHTERS]      = mcdaughterparticle->Trajectory().TotalLength();
      fdaughter_truth_Process[fNDAUGHTERS]          = int(mcdaughterparticle->Trajectory().ProcessToKey(mcdaughterparticle->Process()));
      fdaughter_truth_EndProcess[fNDAUGHTERS]       = int(mcdaughterparticle->Trajectory().ProcessToKey(mcdaughterparticle->EndProcess()));
      std::cout << "Daughter Process = " << (mcdaughterparticle->Process()).c_str() 
		<< " , mother = " << mcdaughterparticle->Mother() 
		<< std::endl;
    }
  }
  else if(daughterShower != 0x0){
    fdaughterIstrack[fNDAUGHTERS]                   = 0;
    fdaughterIsshower[fNDAUGHTERS]                  = 1;
    fdaughterLength[fNDAUGHTERS]                    = daughterShower->Length();
    fdaughterShowerBestPlane[fNDAUGHTERS]           = daughterShower->best_plane();
    fdaughterOpeningAngle[fNDAUGHTERS]              = daughterShower->OpenAngle();
    fdaughterStartPosition[fNDAUGHTERS][0]          = daughterShower->ShowerStart().X();
    fdaughterStartPosition[fNDAUGHTERS][1]          = daughterShower->ShowerStart().Y();
    fdaughterStartPosition[fNDAUGHTERS][2]          = daughterShower->ShowerStart().Z();
    fdaughterStartDirection[fNDAUGHTERS][0]         = daughterShower->Direction().X();
    fdaughterStartDirection[fNDAUGHTERS][1]         = daughterShower->Direction().Y();
    fdaughterStartDirection[fNDAUGHTERS][2]         = daughterShower->Direction().Z();
    if( (daughterShower->Energy()).size() > 0 )
      fdaughterShowerEnergy[fNDAUGHTERS] = daughterShower->Energy()[0]; // thisShower->best_plane()
    if( (daughterShower->MIPEnergy()).size() > 0 )
      fdaughterShowerMIPEnergy[fNDAUGHTERS] = daughterShower->MIPEnergy()[0];
    if( (daughterShower->dEdx()).size() > 0 )
      fdaughterShowerdEdx[fNDAUGHTERS] = daughterShower->dEdx()[0];
  }
  else{
    std::cout << "INFO::Daughter pfParticle is not track or shower!" << std::endl;
    //return;
  }
  
  fdaughterID[fNDAUGHTERS]                          = daughterID;
  // NHits associated with this pfParticle
  fdaughterNHits[fNDAUGHTERS]                       = (pfpUtil.GetPFParticleHits(*daughterParticle,evt,fPFParticleTag)).size();
  // T0
  std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*daughterParticle,evt,fPFParticleTag);
  if(!pfdaughterT0vec.empty())
    fdaughterT0[fNDAUGHTERS] = pfdaughterT0vec[0].Time();
  
  // Increment counter
  fNDAUGHTERS++;

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalTree::FillPrimaryGrandDaughterPFParticle(art::Event const & evt, const recob::PFParticle* gdaughterParticle, int daughterID, int gdaughterID){

  const recob::Track* gdaughterTrack              = pfpUtil.GetPFParticleTrack(*gdaughterParticle,evt, fPFParticleTag,fTrackerTag);
  const recob::Shower* gdaughterShower            = pfpUtil.GetPFParticleShower(*gdaughterParticle,evt,fPFParticleTag,fShowerTag);
  
  if(gdaughterTrack != 0x0){
    fgranddaughterIstrack[fNGRANDDAUGHTERS]                = 1;
    fgranddaughterIsshower[fNGRANDDAUGHTERS]               = 0;
    fgranddaughterTheta[fNGRANDDAUGHTERS]                  = gdaughterTrack->Theta();
    fgranddaughterPhi[fNGRANDDAUGHTERS]                    = gdaughterTrack->Phi();
    fgranddaughterLength[fNGRANDDAUGHTERS]                 = gdaughterTrack->Length();
    fgranddaughterMomentum[fNGRANDDAUGHTERS]               = gdaughterTrack->StartMomentum();
    fgranddaughterEndMomentum[fNGRANDDAUGHTERS]            = gdaughterTrack->EndMomentum();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][0]       = gdaughterTrack->Trajectory().Start().X();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][1]       = gdaughterTrack->Trajectory().Start().Y();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][2]       = gdaughterTrack->Trajectory().Start().Z();
    fgranddaughterEndPosition[fNGRANDDAUGHTERS][0]         = gdaughterTrack->Trajectory().End().X();
    fgranddaughterEndPosition[fNGRANDDAUGHTERS][1]         = gdaughterTrack->Trajectory().End().Y();
    fgranddaughterEndPosition[fNGRANDDAUGHTERS][2]         = gdaughterTrack->Trajectory().End().Z();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][0]      = gdaughterTrack->Trajectory().StartDirection().X();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][1]      = gdaughterTrack->Trajectory().StartDirection().Y();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][2]      = gdaughterTrack->Trajectory().StartDirection().Z();
    fgranddaughterEndDirection[fNGRANDDAUGHTERS][0]        = gdaughterTrack->Trajectory().EndDirection().X();
    fgranddaughterEndDirection[fNGRANDDAUGHTERS][1]        = gdaughterTrack->Trajectory().EndDirection().Y();
    fgranddaughterEndDirection[fNGRANDDAUGHTERS][2]        = gdaughterTrack->Trajectory().EndDirection().Z();
    
    fgranddaughterMomentumByRangeMuon[fNGRANDDAUGHTERS]    = trmom.GetTrackMomentum(gdaughterTrack->Length(),13);
    fgranddaughterMomentumByRangeProton[fNGRANDDAUGHTERS]  = trmom.GetTrackMomentum(gdaughterTrack->Length(),2212);
    
    // Calorimetry
    std::vector<anab::Calorimetry> daughtercalovector = trackUtil.GetRecoTrackCalorimetry(*gdaughterTrack, evt, fTrackerTag, fCalorimetryTag);
    if(daughtercalovector.size() != 3)
      std::cerr << "WARNING::Calorimetry vector size for grand-daughter is = " << daughtercalovector.size() << std::endl;
    
    for(size_t k = 0; k < daughtercalovector.size() && k<3; k++){
      int plane = daughtercalovector[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fgranddaughterKineticEnergy[fNGRANDDAUGHTERS][plane] = daughtercalovector[k].KineticEnergy();
      fgranddaughterRange[fNGRANDDAUGHTERS][plane]         = daughtercalovector[k].Range();
      fgranddaughterTrkPitchC[fNGRANDDAUGHTERS][plane]     = daughtercalovector[k].TrkPitchC();
    }

    // PID
    std::vector<anab::ParticleID> daughterpids = trackUtil.GetRecoTrackPID(*gdaughterTrack, evt, fTrackerTag, fParticleIDTag);
    if(daughterpids.size() != 3)
      std::cerr << "WARNING::PID vector size for grand-daughter is = " << daughterpids.size() << std::endl;
    
    for(size_t k = 0; k < daughterpids.size() && k<3; k++){
      int plane = daughterpids[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fgranddaughterPID_Pdg[fNGRANDDAUGHTERS][plane]          = daughterpids[plane].Pdg();
      fgranddaughterPID_Ndf[fNGRANDDAUGHTERS][plane]          = daughterpids[plane].Ndf();
      fgranddaughterPID_MinChi2[fNGRANDDAUGHTERS][plane]      = daughterpids[plane].MinChi2();
      fgranddaughterPID_DeltaChi2[fNGRANDDAUGHTERS][plane]    = daughterpids[plane].DeltaChi2();
      fgranddaughterPID_Chi2Proton[fNGRANDDAUGHTERS][plane]   = daughterpids[plane].Chi2Proton();
      fgranddaughterPID_Chi2Kaon[fNGRANDDAUGHTERS][plane]     = daughterpids[plane].Chi2Kaon();
      fgranddaughterPID_Chi2Pion[fNGRANDDAUGHTERS][plane]     = daughterpids[plane].Chi2Pion();
      fgranddaughterPID_Chi2Muon[fNGRANDDAUGHTERS][plane]     = daughterpids[plane].Chi2Muon();
      fgranddaughterPID_MissingE[fNGRANDDAUGHTERS][plane]     = daughterpids[plane].MissingE();
      fgranddaughterPID_MissingEavg[fNGRANDDAUGHTERS][plane]  = daughterpids[plane].MissingEavg();
      fgranddaughterPID_PIDA[fNGRANDDAUGHTERS][plane]         = daughterpids[plane].PIDA();
    }
    
    // Get the true mc particle
    const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(*gdaughterTrack, evt, fTrackerTag);
    if(mcdaughterparticle != 0x0){
      fgranddaughter_truth_TrackId[fNGRANDDAUGHTERS]          = mcdaughterparticle->TrackId();
      fgranddaughter_truth_Pdg[fNGRANDDAUGHTERS]              = mcdaughterparticle->PdgCode();
      fgranddaughter_truth_Mother[fNGRANDDAUGHTERS]           = mcdaughterparticle->Mother();
      fgranddaughter_truth_StartPosition[fNGRANDDAUGHTERS][0] = mcdaughterparticle->Vx();
      fgranddaughter_truth_StartPosition[fNGRANDDAUGHTERS][1] = mcdaughterparticle->Vy();
      fgranddaughter_truth_StartPosition[fNGRANDDAUGHTERS][2] = mcdaughterparticle->Vz();
      fgranddaughter_truth_StartPosition[fNGRANDDAUGHTERS][3] = mcdaughterparticle->T();
      fgranddaughter_truth_EndPosition[fNGRANDDAUGHTERS][0]   = mcdaughterparticle->EndX();
      fgranddaughter_truth_EndPosition[fNGRANDDAUGHTERS][1]   = mcdaughterparticle->EndY();
      fgranddaughter_truth_EndPosition[fNGRANDDAUGHTERS][2]   = mcdaughterparticle->EndZ();
      fgranddaughter_truth_EndPosition[fNGRANDDAUGHTERS][3]   = mcdaughterparticle->EndT();
      fgranddaughter_truth_P[fNGRANDDAUGHTERS]                = mcdaughterparticle->P();
      fgranddaughter_truth_Momentum[fNGRANDDAUGHTERS][0]      = mcdaughterparticle->Px();
      fgranddaughter_truth_Momentum[fNGRANDDAUGHTERS][1]      = mcdaughterparticle->Py();
      fgranddaughter_truth_Momentum[fNGRANDDAUGHTERS][2]      = mcdaughterparticle->Pz();
      fgranddaughter_truth_Momentum[fNGRANDDAUGHTERS][3]      = mcdaughterparticle->E();
      fgranddaughter_truth_Pt[fNGRANDDAUGHTERS]               = mcdaughterparticle->Pt();
      fgranddaughter_truth_Mass[fNGRANDDAUGHTERS]             = mcdaughterparticle->Mass();
      fgranddaughter_truth_EndMomentum[fNGRANDDAUGHTERS][0]   = mcdaughterparticle->EndPx();
      fgranddaughter_truth_EndMomentum[fNGRANDDAUGHTERS][1]   = mcdaughterparticle->EndPy();
      fgranddaughter_truth_EndMomentum[fNGRANDDAUGHTERS][2]   = mcdaughterparticle->EndPz();
      fgranddaughter_truth_EndMomentum[fNGRANDDAUGHTERS][3]   = mcdaughterparticle->EndE();
      fgranddaughter_truth_Theta[fNGRANDDAUGHTERS]            = mcdaughterparticle->Momentum().Theta();
      fgranddaughter_truth_Phi[fNGRANDDAUGHTERS]              = mcdaughterparticle->Momentum().Phi();
      fgranddaughter_truth_TotalLength[fNGRANDDAUGHTERS]      = mcdaughterparticle->Trajectory().TotalLength();
      fgranddaughter_truth_Process[fNGRANDDAUGHTERS]          = int(mcdaughterparticle->Trajectory().ProcessToKey(mcdaughterparticle->Process()));
      fgranddaughter_truth_EndProcess[fNGRANDDAUGHTERS]       = int(mcdaughterparticle->Trajectory().ProcessToKey(mcdaughterparticle->EndProcess()));
      std::cout << "GrandDaughter Process = " << (mcdaughterparticle->Process()).c_str() 
		<< " , mother = " << mcdaughterparticle->Mother() 
		<< std::endl;
    }
  }
  else if(gdaughterShower != 0x0){
    fgranddaughterIstrack[fNGRANDDAUGHTERS]                   = 0;
    fgranddaughterIsshower[fNGRANDDAUGHTERS]                  = 1;
    fgranddaughterLength[fNGRANDDAUGHTERS]                    = gdaughterShower->Length();
    fgranddaughterShowerBestPlane[fNGRANDDAUGHTERS]           = gdaughterShower->best_plane();
    fgranddaughterOpeningAngle[fNGRANDDAUGHTERS]              = gdaughterShower->OpenAngle();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][0]          = gdaughterShower->ShowerStart().X();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][1]          = gdaughterShower->ShowerStart().Y();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][2]          = gdaughterShower->ShowerStart().Z();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][0]         = gdaughterShower->Direction().X();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][1]         = gdaughterShower->Direction().Y();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][2]         = gdaughterShower->Direction().Z();
    if( (gdaughterShower->Energy()).size() > 0 )
      fgranddaughterShowerEnergy[fNGRANDDAUGHTERS] = gdaughterShower->Energy()[0]; // thisShower->best_plane()
    if( (gdaughterShower->MIPEnergy()).size() > 0 )
      fgranddaughterShowerMIPEnergy[fNGRANDDAUGHTERS] = gdaughterShower->MIPEnergy()[0];
    if( (gdaughterShower->dEdx()).size() > 0 )
      fgranddaughterShowerdEdx[fNGRANDDAUGHTERS] = gdaughterShower->dEdx()[0];
  }
  else{
    std::cout << "INFO::GrandDaughter pfParticle is not track or shower!" << std::endl;
    //return;
  }
  
  fgranddaughterID[fNGRANDDAUGHTERS]                          = gdaughterID;
  fgranddaughterMotherID[fNGRANDDAUGHTERS]                    = daughterID;
  // NHits associated with this pfParticle
  fgranddaughterNHits[fNGRANDDAUGHTERS]                       = (pfpUtil.GetPFParticleHits(*gdaughterParticle,evt,fPFParticleTag)).size();
  // T0
  std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*gdaughterParticle,evt,fPFParticleTag);
  if(!pfdaughterT0vec.empty())
    fgranddaughterT0[fNGRANDDAUGHTERS] = pfdaughterT0vec[0].Time();
  
  // Increment counter
  fNGRANDDAUGHTERS++;

}

// -----------------------------------------------------------------------------
//void protoana::ProtoDUNEAnalTree::FillCosmicsTree(art::Event const & evt, std::string pfParticleTag){

  // To fill

//}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalTree::Initialise(){
  
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  for(int k=0; k < 3; k++){
    fprimaryVertex[k] = -999.0;
    fdaughterVertex[k] = -999.0;
    fprimaryEndPosition[k] = -999.0;
    fprimaryStartPosition[k] = -999.0;
    fprimaryEndDirection[k] = -999.0;
    fprimaryStartDirection[k] = -999.0;
    fprimaryKineticEnergy[k] = -999.0;
    fprimaryRange[k] = -999.0;
    fprimaryTrkPitchC[k] = -999.0;

    fprimaryPID_Pdg[k] = -999.0;
    fprimaryPID_Ndf[k] = -999.0;
    fprimaryPID_MinChi2[k] = -999.0;
    fprimaryPID_DeltaChi2[k] = -999.0;
    fprimaryPID_Chi2Proton[k] = -999.0;
    fprimaryPID_Chi2Kaon[k] = -999.0;
    fprimaryPID_Chi2Pion[k] = -999.0;
    fprimaryPID_Chi2Muon[k] = -999.0;
    fprimaryPID_MissingE[k] = -999.0;
    fprimaryPID_MissingEavg[k] = -999.0;
    fprimaryPID_PIDA[k] = -999.0;
  }

  fbeamtrigger = -999;
  ftof = -999.0;
  for(int k=0; k < 2; k++){
    fcerenkovPressure[k] = -999.0;
    fcerenkovTime[k] = -999.0; 
    fcerenkovStatus[k] = -999;
  }
  fbeamtrackMomentum = -999.0;
  fbeamtrackEnergy = 999.0;
  fbeamtrackPdg = -999;
  fbeamtrackTime = -999.0;
  fbeamtrackID = -999;
  for(int l=0; l < 3; l++){
    fbeamtrackPos[l] = -999.0;
    fbeamtrackDir[l] = -999.0;
  }
 
  fprimaryIstrack = 0;
  fprimaryIsshower = 0;

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
  fprimaryT0 = -999.0;

  fprimary_truth_TrackId = -999;
  fprimary_truth_Pdg = -999;
  fprimary_truth_Mother = -999;
  fprimary_truth_P = -999.0;
  fprimary_truth_Pt = -999.0;
  fprimary_truth_Mass = -999.0;
  fprimary_truth_Theta = -999.0;
  fprimary_truth_Phi = -999.0;
  fprimary_truth_Process = -999;
  fprimary_truth_Isbeammatched = -999;
  fprimary_truth_NDaughters = -999;
  fprimary_truth_TotalLength = -999;
  for(int l=0; l < 4; l++){
    fprimary_truth_StartPosition[l] = -999.0;
    fprimary_truth_EndPosition[l] = -999.0;
    fprimary_truth_Momentum[l] = -999.0;
    fprimary_truth_EndMomentum[l] = -999.0;
  }

  fNDAUGHTERS = 0;
  fNGRANDDAUGHTERS = 0;
  for(int k=0; k < NMAXDAUGTHERS; k++){
    fdaughterIstrack[k] = -999;
    fdaughterIsshower[k] = -999;
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
      fdaughterTrkPitchC[k][l] = -999.0;

      fdaughterPID_Pdg[k][l] = -999.0;
      fdaughterPID_Ndf[k][l] = -999.0;
      fdaughterPID_MinChi2[k][l] = -999.0;
      fdaughterPID_DeltaChi2[k][l] = -999.0;
      fdaughterPID_Chi2Proton[k][l] = -999.0;
      fdaughterPID_Chi2Kaon[k][l] = -999.0;
      fdaughterPID_Chi2Pion[k][l] = -999.0;
      fdaughterPID_Chi2Muon[k][l] = -999.0;
      fdaughterPID_MissingE[k][l] = -999.0;
      fdaughterPID_MissingEavg[k][l] = -999.0;
      fdaughterPID_PIDA[k][l] = -999.0;
      fgranddaughterVertex[k][l] = -999.0;
    }
    fdaughterOpeningAngle[k] = -999.0;
    fdaughterShowerBestPlane[k] = -999;
    fdaughterShowerEnergy[k] = -999;
    fdaughterShowerMIPEnergy[k] = -999;
    fdaughterShowerdEdx[k] = -999;
    fdaughterMomentumByRangeProton[k] = -999.0;
    fdaughterMomentumByRangeMuon[k] = -999.0;
    fdaughterID[k] = -999;
    fdaughterT0[k] = -999;

    fdaughter_truth_TrackId[k] = -999;
    fdaughter_truth_Pdg[k] = -999;
    fdaughter_truth_Mother[k] = -999;
    fdaughter_truth_P[k] = -999.0;
    fdaughter_truth_Pt[k] = -999.0;
    fdaughter_truth_Mass[k] = -999.0;
    fdaughter_truth_Theta[k] = -999.0;
    fdaughter_truth_Phi[k] = -999.0;
    fdaughter_truth_Process[k] = -999;
    fdaughter_truth_TotalLength[k] = -999.0;
    for(int l=0; l < 4; l++){
      fdaughter_truth_StartPosition[k][l] = -999.0;
      fdaughter_truth_EndPosition[k][l] = -999.0;
      fdaughter_truth_Momentum[k][l] = -999.0;
      fdaughter_truth_EndMomentum[k][l] = -999.0;
    }

    // Grand-daugthers
    fgranddaughterIstrack[k] = -999;
    fgranddaughterIsshower[k] = -999;
    fgranddaughterNHits[k] = -999;
    fgranddaughterTheta[k] = -999.0;
    fgranddaughterPhi[k] = -999.0;
    fgranddaughterLength[k] = -999.0;
    fgranddaughterMomentum[k] = -999.0;
    fgranddaughterEndMomentum[k] = -999.0;
    for(int l=0; l < 3; l++){
      fgranddaughterEndPosition[k][l] = -999.0;
      fgranddaughterStartPosition[k][l] = -999.0;
      fgranddaughterEndDirection[k][l] = -999.0;
      fgranddaughterStartDirection[k][l] = -999.0;
      fgranddaughterKineticEnergy[k][l] = -999.0;
      fgranddaughterRange[k][l] = -999.0;
      fgranddaughterTrkPitchC[k][l] = -999.0;

      fgranddaughterPID_Pdg[k][l] = -999.0;
      fgranddaughterPID_Ndf[k][l] = -999.0;
      fgranddaughterPID_MinChi2[k][l] = -999.0;
      fgranddaughterPID_DeltaChi2[k][l] = -999.0;
      fgranddaughterPID_Chi2Proton[k][l] = -999.0;
      fgranddaughterPID_Chi2Kaon[k][l] = -999.0;
      fgranddaughterPID_Chi2Pion[k][l] = -999.0;
      fgranddaughterPID_Chi2Muon[k][l] = -999.0;
      fgranddaughterPID_MissingE[k][l] = -999.0;
      fgranddaughterPID_MissingEavg[k][l] = -999.0;
      fgranddaughterPID_PIDA[k][l] = -999.0;
    }
    fgranddaughterOpeningAngle[k] = -999.0;
    fgranddaughterShowerBestPlane[k] = -999;
    fgranddaughterShowerEnergy[k] = -999;
    fgranddaughterShowerMIPEnergy[k] = -999;
    fgranddaughterShowerdEdx[k] = -999;
    fgranddaughterMomentumByRangeProton[k] = -999.0;
    fgranddaughterMomentumByRangeMuon[k] = -999.0;
    fgranddaughterID[k] = -999;
    fgranddaughterMotherID[k] = -999;
    fgranddaughterT0[k] = -999;

    fgranddaughter_truth_TrackId[k] = -999;
    fgranddaughter_truth_Pdg[k] = -999;
    fgranddaughter_truth_Mother[k] = -999;
    fgranddaughter_truth_P[k] = -999.0;
    fgranddaughter_truth_Pt[k] = -999.0;
    fgranddaughter_truth_Mass[k] = -999.0;
    fgranddaughter_truth_Theta[k] = -999.0;
    fgranddaughter_truth_Phi[k] = -999.0;
    fgranddaughter_truth_Process[k] = -999;
    fgranddaughter_truth_TotalLength[k] = -999.0;
    for(int l=0; l < 4; l++){
      fgranddaughter_truth_StartPosition[k][l] = -999.0;
      fgranddaughter_truth_EndPosition[k][l] = -999.0;
      fgranddaughter_truth_Momentum[k][l] = -999.0;
      fgranddaughter_truth_EndMomentum[k][l] = -999.0;
    }

  }


}

DEFINE_ART_MODULE(protoana::ProtoDUNEAnalTree)

