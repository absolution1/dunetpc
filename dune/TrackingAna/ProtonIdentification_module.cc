////////////////////////////////////////////////////////////////////////
// Class:       ProtonIdentification
// Module Type: analyzer
// File:        ProtonIdentification_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
//
//  Using AnaTree as base, have written a module to compare reconstruction
//  algorithm efficiencies.
//
//  Thomas Karl Warburton
//  k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TLorentzVector.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TNtuple.h"
#include "TFile.h"

//standard library includes
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits> // std::numeric_limits<>

namespace ProtonIdentification {
  class ProtonIdentification;
}

class ProtonIdentification::ProtonIdentification : public art::EDAnalyzer {
public:
  explicit ProtonIdentification(fhicl::ParameterSet const & p);
  virtual ~ProtonIdentification();

  void analyze(art::Event const & e) override;

  void beginRun(art::Run& run);
  void beginJob() override;
  void endJob();
  void endRun();
  //void reconfigure(fhicl::ParameterSet const & p) override;
  
  // ----------- Declare my structs ----------
  struct CheatHists {
    //EfficHists();
    // EfficHists(const std::string& subdir);
    TH1D* dEdx;
    TH1D* ResRange;
    TH2D* dEdx_ResRange;
    TH1D* PIDA;
    TH1D* PIDA_Plane0;
    TH1D* PIDA_Plane1;
    TH1D* PIDA_Plane2;
    TH1D* DistFromEdge;
    TH2D* XYPlaneStart;
    TH2D* XZPlaneStart;
    TH2D* YZPlaneStart;
    TH2D* XYPlaneEnd;
    TH2D* XZPlaneEnd;
    TH2D* YZPlaneEnd;
    TH2D* PIDA_nHits;
} Cheat_Proton, Cheat_Muon, Cheat_All;
  // ----------- Declare my structs ----------

private:
  
  void ResetVars();
  void MCTruthInformation ( const simb::MCParticle *particle );

  void TrackBoundaries ( TVector3 larStart, TVector3 larEnd );

  double CalcPIDA ( std::vector<const anab::Calorimetry*> calos );
  
  void CheatPIDAFiller ( struct ProtonIdentification::ProtonIdentification::CheatHists *HistPtr, 
			 std::vector<const anab::Calorimetry*> calos, TVector3 larStart, TVector3 larEnd, int TrackHits );

  double boundaries[6]; // For use in finding detector boundaries.  
  TTree *fCheatTree, *fRecoTree;

  double MCTruthT0, MCTruthTickT0, MCTruthTrackID;
  double PhotonCounterT0, PhotonCounterTickT0, PhotonCounterID;
  double TrackTheta_XZ, TrackTheta_YZ, TrackEta_XY, TrackEta_ZY, TrackTheta, TrackPhi;
  double TrackLength; 
  std::vector<double> trackStart, trackEnd;

  double MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi;
  double MCTPCLength, MCEnergy, MCEnergyDeposited;
  bool   MCStartInTPC, MCEndInTPC, MCStops;
  int    MCPdgCode, MCTrackId, MatchedTrackID, MCContainment, MCParticleKept;

  double PIDA, PIDA_Plane0, PIDA_Plane1, PIDA_Plane2, TrackDistFromEdge, AvdEdx;
  double CorrectedStartX, CorrectedStartY, CorrectedStartZ, CorrectedEndX, CorrectedEndY, CorrectedEndZ;
  bool   startsonboundary, endsonboundary, trackstops;
  int    RecoContainment, CounterID, TrackHits;

  bool AllChargedTrack;

  int TrueMuons=0, TrueElectrons=0, TrueProtons=0, TrueGammas=0, TrueOthers=0;
  int AllBad = 0, MuonBad = 0, ProtonBad = 0, GammaBad = 0, ElectronBad=0, OtherBad = 0;
  int AllAll = 0, MuonAll = 0, ProtonAll = 0, GammaAll = 0, ElectronAll=0, OtherAll = 0;
  int AllRange = 0, MuonRange = 0, ProtonRange = 0, GammaRange = 0, ElectronRange=0, OtherRange = 0;
  // Handles
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTracker> bktrk;
  detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const *ts = lar::providerFrom<detinfo::DetectorClocksService>();
  
  double XDriftVelocity      = detprop->DriftVelocity()*1e-3; //cm/ns
  double WindowSize          = detprop->NumberTimeSamples() * ts->TPCClock().TickPeriod() * 1e3;
  // Parameter List
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fMCTruthT0ModuleLabel;
  std::string fPhotonT0ModuleLabel;
  std::string fCounterT0ModuleLabel;
  std::string fCalorimetryModuleLabel;
  double PIDApower, fBoundaryEdge;

  TH2D* PIDA_vs_Hits_good;
  TH2D* PIDA_vs_Hits_bad;
 };
// ********************************** Begin Run *******************************************************
void ProtonIdentification::ProtonIdentification::beginRun(art::Run& run) {

}
// *********************************** Begin Job ********************************************************
void ProtonIdentification::ProtonIdentification::beginJob()
{
  // --------- Working out detector dimensions ----------------
  double TempBounds[6];
  int NCryo = geom->Ncryostats();
  for ( int c1 = 0; c1 < NCryo; c1++ ) { // loop through TPC's
    geom->CryostatBoundaries(TempBounds, c1); // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )
    if ( boundaries[0] > TempBounds [0] ) boundaries[0] = TempBounds [0];
    if ( boundaries[1] < TempBounds [1] ) boundaries[1] = TempBounds [1];
    if ( boundaries[2] > TempBounds [2] ) boundaries[2] = TempBounds [2];
    if ( boundaries[3] < TempBounds [3] ) boundaries[3] = TempBounds [3];
    if ( boundaries[4] > TempBounds [4] ) boundaries[4] = TempBounds [4];
    if ( boundaries[5] < TempBounds [5] ) boundaries[5] = TempBounds [5];
  }
  //****************** The cryo boundaries are much larger than TPC boundaries.....
  boundaries[0] = -35.;
  boundaries[1] = 221.;
  boundaries[2] = -85.;
  boundaries[3] = 115.;
  boundaries[4] = -2.;
  boundaries[5] = 157.;
  std::cout << boundaries[0] << " " << boundaries[1] << " " << boundaries[2] << " " << boundaries[3] << " " <<boundaries[4] << " " << boundaries[5] << std::endl;
  
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fCheatTree = tfs->make<TTree>("CheatedTree","analysis tree");
  fCheatTree->Branch("MCTPCLength"      ,&MCTPCLength      ,"MCTPCLength/D"      );
  fCheatTree->Branch("MatchedTrackLength",&TrackLength     ,"MatchedTrackLength/D");  
  fCheatTree->Branch("MCEnergy"         ,&MCEnergy         ,"MCEnergy/D"         );
  fCheatTree->Branch("MCEnergyDeposited",&MCEnergyDeposited,"MCEnergyDeposited/D");
  fCheatTree->Branch("MCTheta"          ,&MCTheta          ,"MCTheta/D"          );
  fCheatTree->Branch("MCTheta_XZ"       ,&MCTheta_XZ       ,"MCTheta_XZ/D"       );
  fCheatTree->Branch("MCTheta_YZ"       ,&MCTheta_YZ       ,"MCTheta/D"          );
  fCheatTree->Branch("MCPhi"            ,&MCPhi            ,"MCPhi/D"            );
  fCheatTree->Branch("MCEta_XY"         ,&MCEta_XY         ,"MCEta_XY/D"         );
  fCheatTree->Branch("MCEta_ZY"         ,&MCEta_ZY         ,"MCEta_ZY/D"         );
  fCheatTree->Branch("MCPdgCode"        ,&MCPdgCode        ,"MCPdgCode/I"        );
  fCheatTree->Branch("MCTrackId"        ,&MCTrackId        ,"MCTrackId/I"        );
  fCheatTree->Branch("MatchedTrackID"   ,&MatchedTrackID   ,"MatchedTrackID/I"   );
  fCheatTree->Branch("MCStartInTPC"     ,&MCStartInTPC     ,"MCStartInTPC/B"     );
  fCheatTree->Branch("MCEndInTPC"       ,&MCEndInTPC       ,"MCEndInTPC/B"       );
  fCheatTree->Branch("MCStops"          ,&MCStops          ,"MCStops/B"          );
  fCheatTree->Branch("MCContainment"    ,&MCContainment    ,"MCContainment/I"    );
  fCheatTree->Branch("MCParticleKept"   ,&MCParticleKept   ,"MCParticleKept/I"   );
  
  fRecoTree = tfs->make<TTree>("ReconstructedTree","analysis tree");
  fRecoTree->Branch("MCPdgCode"        ,&MCPdgCode        ,"MCPdgCode/I"        );
  fRecoTree->Branch("PIDA"             ,&PIDA             ,"PIDA/D"             );
  fRecoTree->Branch("PIDA_Plane0"      ,&PIDA_Plane0      ,"PIDA_Plane0/D"      );
  fRecoTree->Branch("PIDA_Plane1"      ,&PIDA_Plane1      ,"PIDA_Plane1/D"      );
  fRecoTree->Branch("PIDA_Plane2"      ,&PIDA_Plane2      ,"PIDA_Plane2/D"      );
  fRecoTree->Branch("AvdEdx"           ,&AvdEdx           ,"AvdEdx/D"           );
  fRecoTree->Branch("startsonboundary" ,&startsonboundary ,"startsonboundary/B" );
  fRecoTree->Branch("endsonboundary"   ,&endsonboundary   ,"endsonboundary/B"   );
  fRecoTree->Branch("trackstops"       ,&trackstops       ,"trackstops/B"       );
  fRecoTree->Branch("RecoContainment"  ,&RecoContainment  ,"RecoContainment/I"  );
  fRecoTree->Branch("TrackLength"      ,&TrackLength      ,"TrackLength/D"      );
  fRecoTree->Branch("TrackHits"        ,&TrackHits        ,"TrackHits/I"        );
  fRecoTree->Branch("CorrectedStartX"  ,&CorrectedStartX  ,"CorrectedStartX/D"  );
  fRecoTree->Branch("CorrectedStartY"  ,&CorrectedStartY  ,"CorrectedStartY/D"  );
  fRecoTree->Branch("CorrectedStartZ"  ,&CorrectedStartZ  ,"CorrectedStartZ/D"  );
  fRecoTree->Branch("CorrectedEndX"    ,&CorrectedEndX    ,"CorrectedEndX/D"    );
  fRecoTree->Branch("CorrectedEndY"    ,&CorrectedEndY    ,"CorrectedEndY/D"    );
  fRecoTree->Branch("CorrectedEndZ"    ,&CorrectedEndZ    ,"CorrectedEndZ/D"    );
  fRecoTree->Branch("CounterID"        ,&CounterID        ,"CounterID/D"        );

  // ----- CHEATED HISTOGRAMS ------- Protons first -----
  Cheat_Proton.dEdx     = tfs->make<TH1D>("Proton_dEdx", "dEdx values for protons in the detector; dEdx (MeV cm); Number", 100, 0, 50);
  Cheat_Proton.ResRange = tfs->make<TH1D>("Proton_ResRange", "Residual range for protons in the detector; Residual Range (cm); Number", 100, 0, 50);
  Cheat_Proton.dEdx_ResRange = tfs->make<TH2D>("Proton_dEdx_ResRange", 
					       "Residual range against dEdx for protons in the detector; Residual Range (cm); dEdx (MeV cm)",
					       100, 0, 50, 100, 0, 50);
  Cheat_Proton.PIDA     = tfs->make<TH1D>("Proton_PIDA", "PIDA values for protons in the detector; PIDA, Number", 100, 0, 30);
  Cheat_Proton.PIDA_Plane0 = tfs->make<TH1D>("Proton_PIDA_Plane0", "PIDA values for protons in the detector on plane 0; PIDA; Number", 100, 0, 30);
  Cheat_Proton.PIDA_Plane1 = tfs->make<TH1D>("Proton_PIDA_Plane1", "PIDA values for protons in the detector on plane 1; PIDA; Number", 100, 0, 30);
  Cheat_Proton.PIDA_Plane2 = tfs->make<TH1D>("Proton_PIDA_Plane2", "PIDA values for protons in the detector on plane 2; PIDA; Number", 100, 0, 30);
  
  Cheat_Proton.DistFromEdge= tfs->make<TH1D>("Proton_DistFromEdge","Distance of reconstructed track from TPC edge; Distance (cm); Number", 100, 0, 80);
  Cheat_Proton.XYPlaneStart= tfs->make<TH2D>("Proton_XYPlaneStart","Starting position of track in XY plane; X (cm); Y (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[2], boundaries[3] );
  Cheat_Proton.XZPlaneStart= tfs->make<TH2D>("Proton_XZPlaneStart","Starting position of track in XZ plane; X (cm); Z (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[4], boundaries[5] );
  Cheat_Proton.YZPlaneStart= tfs->make<TH2D>("Proton_YZPlaneStart","Starting position of track in YZ plane; Y (cm); Z (cm);", 100, boundaries[2], boundaries[3], 100, boundaries[4], boundaries[5] );
  Cheat_Proton.XYPlaneEnd= tfs->make<TH2D>("Proton_XYPlaneEnd","Ending position of track in XY plane; X (cm); Y (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[2], boundaries[3] );
  Cheat_Proton.XZPlaneEnd= tfs->make<TH2D>("Proton_XZPlaneEnd","Ending position of track in XZ plane; X (cm); Z (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[4], boundaries[5] );
  Cheat_Proton.YZPlaneEnd= tfs->make<TH2D>("Proton_YZPlaneEnd","Ending position of track in YZ plane; Y (cm); Z (cm);", 100, boundaries[2], boundaries[3], 100, boundaries[4], boundaries[5] );
  Cheat_Proton.PIDA_nHits= tfs->make<TH2D>("Proton_PIDA_nHits","PIDA value versus number of hits; PIDA value; number of hits;", 100, 0, 30, 100, 0, 500 );

  // ------ Now for Muons ------
  Cheat_Muon.dEdx     = tfs->make<TH1D>("Muon_dEdx", "dEdx values for muons in the detector; dEdx (MeV cm); Number", 100, 0, 50);
  Cheat_Muon.ResRange = tfs->make<TH1D>("Muon_ResRange", "Residual range for muons in the detector; Residual Range (cm); Number", 100, 0, 50);
  Cheat_Muon.dEdx_ResRange = tfs->make<TH2D>("Muon_dEdx_ResRange", 
				       "Residual range against dEdx for muons in the detector; Residual Range (cm); dEdx (MeV cm)",
				       100, 0, 50, 100, 0, 50);
  Cheat_Muon.PIDA     = tfs->make<TH1D>("Muon_PIDA", "PIDA values for muons in the detector; PIDA, Number", 100, 0, 30);
  Cheat_Muon.PIDA_Plane0 = tfs->make<TH1D>("Muon_PIDA_Plane0", "PIDA values for protons in the detector on plane 0; PIDA; Number", 100, 0, 30);
  Cheat_Muon.PIDA_Plane1 = tfs->make<TH1D>("Muon_PIDA_Plane1", "PIDA values for protons in the detector on plane 1; PIDA; Number", 100, 0, 30);
  Cheat_Muon.PIDA_Plane2 = tfs->make<TH1D>("Muon_PIDA_Plane2", "PIDA values for protons in the detector on plane 2; PIDA; Number", 100, 0, 30);
  
  Cheat_Muon.DistFromEdge= tfs->make<TH1D>("Muon_DistFromEdge","Distance of reconstructed track from TPC edge; Distance (cm); Number", 100, 0, 150);
  Cheat_Muon.XYPlaneStart= tfs->make<TH2D>("Muon_XYPlaneStart","Starting position of track in XY plane; X (cm); Y (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[2], boundaries[3] );
  Cheat_Muon.XZPlaneStart= tfs->make<TH2D>("Muon_XZPlaneStart","Starting position of track in XZ plane; X (cm); Z (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[4], boundaries[5] );
  Cheat_Muon.YZPlaneStart= tfs->make<TH2D>("Muon_YZPlaneStart","Starting position of track in YZ plane; Y (cm); Z (cm);", 100, boundaries[2], boundaries[3], 100, boundaries[4], boundaries[5] );
  Cheat_Muon.XYPlaneEnd= tfs->make<TH2D>("Muon_XYPlaneEnd","Ending position of track in XY plane; X (cm); Y (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[2], boundaries[3] );
  Cheat_Muon.XZPlaneEnd= tfs->make<TH2D>("Muon_XZPlaneEnd","Ending position of track in XZ plane; X (cm); Z (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[4], boundaries[5] );
  Cheat_Muon.YZPlaneEnd= tfs->make<TH2D>("Muon_YZPlaneEnd","Ending position of track in YZ plane; Y (cm); Z (cm);", 100, boundaries[2], boundaries[3], 100, boundaries[4], boundaries[5] );
  Cheat_Muon.PIDA_nHits= tfs->make<TH2D>("Muon_PIDA_nHits","PIDA value versus number of hits; PIDA value; number of hits;", 100, 0, 30, 100, 0, 500 );
  // ------ One for all particle ------
  Cheat_All.dEdx     = tfs->make<TH1D>("All_dEdx", "dEdx values for muons in the detector; dEdx (MeV cm); Number", 100, 0, 50);
  Cheat_All.ResRange = tfs->make<TH1D>("All_ResRange", "Residual range for muons in the detector; Residual Range (cm); Number", 100, 0, 50);
  Cheat_All.dEdx_ResRange = tfs->make<TH2D>("All_dEdx_ResRange", 
				       "Residual range against dEdx for muons in the detector; Residual Range (cm); dEdx (MeV cm)",
				       100, 0, 50, 100, 0, 50);
  Cheat_All.PIDA     = tfs->make<TH1D>("All_PIDA", "PIDA values for muons in the detector; PIDA, Number", 100, 0, 30);
  Cheat_All.PIDA_Plane0 = tfs->make<TH1D>("All_PIDA_Plane0", "PIDA values for protons in the detector on plane 0; PIDA; Number", 100, 0, 30);
  Cheat_All.PIDA_Plane1 = tfs->make<TH1D>("All_PIDA_Plane1", "PIDA values for protons in the detector on plane 1; PIDA; Number", 100, 0, 30);
  Cheat_All.PIDA_Plane2 = tfs->make<TH1D>("All_PIDA_Plane2", "PIDA values for protons in the detector on plane 2; PIDA; Number", 100, 0, 30);
  
  Cheat_All.DistFromEdge= tfs->make<TH1D>("All_DistFromEdge","Distance of reconstructed track from TPC edge; Distance (cm); Number", 100, 0, 150);
  Cheat_All.XYPlaneStart= tfs->make<TH2D>("All_XYPlaneStart","Starting position of track in XY plane; X (cm); Y (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[2], boundaries[3] );
  Cheat_All.XZPlaneStart= tfs->make<TH2D>("All_XZPlaneStart","Starting position of track in XZ plane; X (cm); Z (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[4], boundaries[5] );
  Cheat_All.YZPlaneStart= tfs->make<TH2D>("All_YZPlaneStart","Starting position of track in YZ plane; Y (cm); Z (cm);", 100, boundaries[2], boundaries[3], 100, boundaries[4], boundaries[5] );
  Cheat_All.XYPlaneEnd= tfs->make<TH2D>("All_XYPlaneEnd","Ending position of track in XY plane; X (cm); Y (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[2], boundaries[3] );
  Cheat_All.XZPlaneEnd= tfs->make<TH2D>("All_XZPlaneEnd","Ending position of track in XZ plane; X (cm); Z (cm);", 100, boundaries[0], boundaries[1], 100, boundaries[4], boundaries[5] );
  Cheat_All.YZPlaneEnd= tfs->make<TH2D>("All_YZPlaneEnd","Ending position of track in YZ plane; Y (cm); Z (cm);", 100, boundaries[2], boundaries[3], 100, boundaries[4], boundaries[5] );
  Cheat_All.PIDA_nHits= tfs->make<TH2D>("All_PIDA_nHits","PIDA value versus number of hits; PIDA value; number of hits;", 100, 0, 30, 100, 0, 500 );
  
  PIDA_vs_Hits_good = tfs->make<TH2D>("PIDA_vs_Hits_good","PIDA vs Hits; PIDA (log 10); Number of hits", 100,0,1.5, 100,0,500);
  PIDA_vs_Hits_bad = tfs->make<TH2D>("PIDA_vs_Hits_bad","PIDA vs Hits; PIDA (log 10); Number of hits", 100,1.5,10, 100,0,500);

// ----- Non-Cheated Histograms ------- From here on in --------

}
// ************************************ End Job *********************************************************
void ProtonIdentification::ProtonIdentification::endJob() {
  std::cout << "\nFinished all the events..."
	    << "\nMonte Carlo's in the TPC had a total of " << TrueMuons << " Muons, " << TrueProtons << " Protons, " << TrueElectrons << " Electrons, " << TrueGammas << " Gammas, and " << TrueOthers << " others." 
	    << "\nHad a total of " << AllAll << " tracks. Comprised of Muons " << MuonAll << ", Protons " << ProtonAll << ", Gammas " << GammaAll << ", Electrons " << ElectronBad << ", Others " << OtherAll  
	    << "\n" << AllBad << " of which had high PIDA's, Muons " << MuonBad << ", Protons " << ProtonBad << ", Gammas " << GammaBad << ", Electrons " << ElectronBad << ", Others " << OtherBad
	    << "\n" << AllRange << " tracks were within PIDA range, Muons " << MuonRange << ", Protons " << ProtonRange << ", Gammas " << GammaRange << ", ElectronRange " << ElectronRange << ", Others " << OtherRange 
	    << std::endl;
}
// ************************************ End Run *********************************************************
void ProtonIdentification::ProtonIdentification::endRun() {
}

// ********************************** pset param *******************************************************
ProtonIdentification::ProtonIdentification::ProtonIdentification(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fHitsModuleLabel         ( pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel        ( pset.get< std::string >("TrackModuleLabel"))
  , fMCTruthT0ModuleLabel    ( pset.get< std::string >("MCTruthT0ModuleLabel"))
  , fPhotonT0ModuleLabel     ( pset.get< std::string >("PhotonT0ModuleLabel"))
  , fCounterT0ModuleLabel    ( pset.get< std::string >("CounterT0ModuleLabel"))
  , fCalorimetryModuleLabel  ( pset.get< std::string >("CalorimetryModuleLabel"))
  , PIDApower                ( pset.get< double      >("PIDApower"))
  , fBoundaryEdge            ( pset.get< double      >("BoundaryEdge"))
{

}
// ******************************************************************************************************
ProtonIdentification::ProtonIdentification::~ProtonIdentification()
{
  // Clean up dynamic memory and other resources here.

}
// ************************************ Analyse *********************************************************
void ProtonIdentification::ProtonIdentification::analyze(art::Event const & evt)
{
  //std::cout << "\n\n************* New Event / Module running *************\n\n" << std::endl;
  // Implementation of required member function here. 
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
  
  art::Handle< std::vector<raw::ExternalTrigger> > trigListHandle;
  std::vector<art::Ptr<raw::ExternalTrigger> > triglist;
  if (evt.getByLabel(fCounterT0ModuleLabel,trigListHandle))
    art::fill_ptr_vector(triglist, trigListHandle);  
  //std::cout << "There were " << triglist.size() << " counter hits in this event." << std::endl;
  
  const sim::ParticleList& plist = bktrk->ParticleList();
  // Quickly get total number of MC Particles...
  for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    simb::MCParticle *particle = ipar->second;
    bool InTPC = false;
    for ( int a=0; a<(int)particle->NumberTrajectoryPoints(); ++a ) {
      if (!InTPC) {
	const TLorentzVector& tmpPosition=particle->Position(a);
	double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
	geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
	if (tpcid.isValid) { 
	  InTPC = true;
	  if (abs(particle->PdgCode() == 13 )) ++TrueMuons;
	  else if (abs(particle->PdgCode() == 11 )) ++TrueElectrons;
	  else if (abs(particle->PdgCode() == 22 )) ++TrueGammas;
	  else if (abs(particle->PdgCode() == 2212 )) ++TrueProtons;
	  else ++TrueOthers;
	}
      }
    }
  }
  
  art::Handle< std::vector<recob::Track> > trackh;
  evt.getByLabel(fTrackModuleLabel, trackh);
  
  art::Handle< std::vector< art::PtrVector < recob::Track > > > trackvh;
  evt.getByLabel(fTrackModuleLabel, trackvh);

  if ( trackListHandle.isValid() ) { // Check that trackListHandle is valid.....
    art::FindManyP<recob::Hit>        fmht   (trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::T0>           fmt0   (trackListHandle, evt, fMCTruthT0ModuleLabel);
    art::FindMany<anab::T0>           fmphot (trackListHandle, evt, fPhotonT0ModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal  (trackListHandle, evt, fCalorimetryModuleLabel);
    int ntracks_reco=tracklist.size();      
    
    for(int Track=0; Track < ntracks_reco; ++Track){
      // *****************************************************************************
      // T0 stuff - So can correct X Start/End positions and identify MCParticle!!
      // *****************************************************************************
      double TickT0 = 0;
      if ( fmt0.isValid() ) {
	std::vector<const anab::T0*> T0s = fmt0.at(Track);
	for (size_t t0size =0; t0size < T0s.size(); t0size++) {
	  MCTruthT0      = T0s[t0size]->Time();
	  MCTruthTickT0  = MCTruthT0 / detprop->SamplingRate();
	  MCTruthTrackID = T0s[t0size]->TriggerBits();
	} // T0 size
	TickT0 = MCTruthTickT0;
      } // T0 valid
      if ( fmphot.isValid() ) {
	std::vector<const anab::T0*> PhotT0 = fmphot.at(Track);
	for (size_t T0it=0; T0it<PhotT0.size(); ++T0it) {
	  PhotonCounterT0     = PhotT0[T0it]->Time();
	  PhotonCounterTickT0 = PhotonCounterT0 / detprop->SamplingRate();
	  PhotonCounterID     = PhotT0[T0it]->TriggerBits();
	}
	TickT0 = PhotonCounterTickT0;
      } 
      // ************** END T0 stuff ***************
      
      // ******************************************************************************************
      // Correct X and get track length etc now that we have matched a Track with an MCParticle!!
      // ******************************************************************************************
      // ----- Hit Stuff ------
      std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(Track);
      double Hit_Size = allHits.size();
      art::Ptr<recob::Track> ptrack(trackh, Track);
      const recob::Track& track = *ptrack;
      
      CounterID = 0;
      int EarlyTick = std::min (allHits[0]->StartTick(), allHits[Hit_Size-1]->StartTick());
      for (size_t CLoop=0; CLoop < triglist.size(); ++CLoop) {
	if (triglist[CLoop]->GetTrigID() < 105 ) continue;
	int TickTime = triglist[CLoop]->GetTrigTime() * 15.625 / 500; // Conv from nova ticks ( 15.625 ns ) to TPC ticks ( 500 ns ).
	int TickDiff = EarlyTick - TickTime;
	if (TickDiff < 3200 && TickDiff > 0) CounterID = triglist[CLoop]->GetTrigID();
      }
      // ---- Get lengths and angles.
      TrackLength = track.Length();
      TrackHits = track.NumberTrajectoryPoints();
      if(TrackHits > 0) {
	TVector3 dir = track.VertexDirection();
	TrackTheta_XZ = std::atan2(dir.X(), dir.Z());
	TrackTheta_YZ = std::atan2(dir.Y(), dir.Z());
	TrackEta_XY   = std::atan2(dir.X(), dir.Y());
	TrackEta_ZY   = std::atan2(dir.Z(), dir.Y());
	TrackTheta    = dir.Theta();
	TrackPhi      = dir.Phi();
      }
      // ---- Correct X positions!
      std::vector < TVector3 > CorrectedLocations;
      for ( int point=0; point < TrackHits; ++point ) {
	const TVector3 ThisLoc = track.LocationAtPoint(point);
	TVector3 CorrectLoc = ThisLoc;
	CorrectLoc[0] = CorrectLoc[0] - detprop->ConvertTicksToX( TickT0, allHits[Hit_Size-(1+point)]->WireID().Plane, allHits[Hit_Size-(1+point)]->WireID().TPC, allHits[Hit_Size-(1+point)]->WireID().Cryostat );
	CorrectedLocations.push_back(CorrectLoc);
      }
      CorrectedStartX = CorrectedLocations[0][0];
      CorrectedStartY = CorrectedLocations[0][1];
      CorrectedStartZ = CorrectedLocations[0][2];
      CorrectedEndX   = CorrectedLocations[TrackHits-1][0];
      CorrectedEndY   = CorrectedLocations[TrackHits-1][1];
      CorrectedEndZ   = CorrectedLocations[TrackHits-1][2];
      // **************************************************************
      // Determine what kind of particle actually caused the track.....
      for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
	simb::MCParticle *particle = ipar->second;
	if ( particle->TrackId() != MCTruthTrackID ) continue;
	MatchedTrackID = Track;
	MCParticleKept = 0;
	// ---- Get MCTruth Information and check that MCParticle goes in TPC ---
	MCTruthInformation ( particle );
	
	++AllAll;
	if (fabs(MCPdgCode) == 13  ) ++MuonAll;
	else if (fabs(MCPdgCode) == 11 ) ++ElectronAll;
	else if (MCPdgCode == 2212 ) ++ProtonAll;
	else if (MCPdgCode == 22   ) ++GammaAll;
	else ++OtherAll;
      } // Check what particle caused this track....
      // ****************************************************************

      // Want to select only tracks which stop in the detector.
      TrackBoundaries( CorrectedLocations[0], CorrectedLocations[TrackHits-1]);
      //if (RecoContainment != 0) continue;
     
      // ---- Make a vector of calorimetry objects...
      std::vector<const anab::Calorimetry*> calos = fmcal.at(Track);
      PIDA = CalcPIDA ( calos );
      int QuickThresh = 25;
      if (PIDA < 25 ) PIDA_vs_Hits_good -> Fill (log10(PIDA), TrackHits);
      else {
	PIDA_vs_Hits_bad  -> Fill (log10(PIDA), TrackHits);
      }
     
      // ************************************************************
      // Now I've done my reco cuts, see what particles I have left.
      // ************************************************************
      MCParticleKept = 1;
      if (PIDA > QuickThresh ) {
	std::cout << "\nTrack " << Track << " has a PIDA value of " << PIDA << ", " << TrackHits << " hits, PdGCode " << MCPdgCode << " , MCTrackID " << MCTrackId << std::endl;
	++AllBad;
	if (fabs(MCPdgCode) == 13  ) ++MuonBad;
	else if (fabs(MCPdgCode) == 11) ++ElectronBad;
	else if (MCPdgCode == 2212 ) ++ProtonBad;
	else if (MCPdgCode == 22   ) ++GammaBad;
	else ++OtherBad;
      }
      if (PIDA > 15 && PIDA < 25 ) {
	++AllRange;
	if (fabs(MCPdgCode) == 13  ) ++MuonRange;
	else if (fabs(MCPdgCode) == 11) ++ElectronRange;
	else if (MCPdgCode == 2212 ) ++ProtonRange;
	else if (MCPdgCode == 22   ) ++GammaRange;
	  else ++OtherRange;
      }
      
      CheatPIDAFiller( &Cheat_All, calos, CorrectedLocations[0], CorrectedLocations[TrackHits-1], TrackHits );
      if ( fabs(MCPdgCode) == 13 ) CheatPIDAFiller( &Cheat_Muon, calos, CorrectedLocations[0], CorrectedLocations[TrackHits-1], TrackHits );
      else if ( MCPdgCode == 2212     ) {
	//if ( particle -> NumberDaughters() != 0 ) continue;
	CheatPIDAFiller( &Cheat_Proton, calos, CorrectedLocations[0], CorrectedLocations[TrackHits-1], TrackHits );
      }
      // ******** Fill Tree for each MCParticle **********
      //std::cout << "Average dEdx " << AvdEdx << std::endl;
      fCheatTree->Fill();
      fRecoTree->Fill();
    } // Loop over Tracks    
  } // if trackListHandle.isValid()
} // Analyse

// ******************************** Calc PIDA  ****************************************************
double ProtonIdentification::ProtonIdentification::CalcPIDA ( std::vector<const anab::Calorimetry*> calos ) {
  double PIDA  = 0, dEdxSum = 0;
  int UsedHits = 0, TotHits = 0;
  for ( int Plane=0; Plane < (int)calos.size(); ++Plane ) { // Loop through planes
    double PlanePIDA=0; int PlaneHits=0;
    for ( int PlaneHit=0; PlaneHit < (int)calos[Plane]->dEdx().size(); ++PlaneHit ) { // loop through hits on each plane
      dEdxSum += calos[Plane]->dEdx()[PlaneHit];
      ++TotHits;
      if ( calos[Plane]->ResidualRange()[PlaneHit] < 30 ) { // Only want PIDA for last 30 cm
	PIDA += calos[Plane]->dEdx()[PlaneHit] * pow(calos[Plane]->ResidualRange()[PlaneHit], PIDApower ); 
	++UsedHits;
	PlanePIDA += calos[Plane]->dEdx()[PlaneHit] * pow(calos[Plane]->ResidualRange()[PlaneHit], PIDApower );
	++PlaneHits;
      } // If ResRange < 30 cm
    } // Loop over hits on each plane
    PlanePIDA = PlanePIDA/PlaneHits;
    if (Plane == 0 ) PIDA_Plane0 = PlanePIDA;
    else if (Plane == 1 ) PIDA_Plane1 = PlanePIDA;
    else if (Plane == 2 ) PIDA_Plane2 = PlanePIDA;
  } // Loop over planes
  
  if ( UsedHits != 0 ) // If had any hits, work out PIDA and calculate
    PIDA = PIDA / UsedHits;
  AvdEdx = dEdxSum / TotHits;
  return PIDA;
} // CalcPIDA

// ******************************** Cheat PIDA Filler ****************************************************
void ProtonIdentification::ProtonIdentification::CheatPIDAFiller ( struct ProtonIdentification::ProtonIdentification::CheatHists *HistPtr, 
								   std::vector<const anab::Calorimetry*> calos, TVector3 larStart, TVector3 larEnd, int TrackHits ) {
  for ( int Plane=0; Plane < (int)calos.size(); ++Plane ) { // Loop through planes
    double TempPIDA = 0; int UsedHits = 0;
    for ( int PlaneHit=0; PlaneHit < (int)calos[Plane]->dEdx().size(); ++PlaneHit ) { // loop through hits on each plane
      HistPtr->dEdx          -> Fill( calos[Plane]->dEdx()[PlaneHit]          );
      HistPtr->ResRange      -> Fill( calos[Plane]->ResidualRange()[PlaneHit] );
      HistPtr->dEdx_ResRange -> Fill( calos[Plane]->ResidualRange()[PlaneHit], calos[Plane]->dEdx()[PlaneHit] );
      if ( calos[Plane]->ResidualRange()[PlaneHit] < 30 ) { // Only want PIDA for last 30 cm
	TempPIDA += calos[Plane]->dEdx()[PlaneHit] * pow(calos[Plane]->ResidualRange()[PlaneHit], PIDApower ); 
	++UsedHits;
      } // If ResRange < 30 cm
    } // Loop over hits on each plane
    if ( UsedHits != 0 ) {
      TempPIDA = TempPIDA / UsedHits;
      if ( Plane == 0 )      HistPtr->PIDA_Plane0 -> Fill ( TempPIDA );
      else if ( Plane == 1 ) HistPtr->PIDA_Plane1 -> Fill ( TempPIDA );
      else if ( Plane == 2 ) HistPtr->PIDA_Plane2 -> Fill ( TempPIDA );
    } // If had any hits on this plane, work out PIDA and calculate
  } // Loop over planes
  HistPtr->PIDA -> Fill( PIDA );
  HistPtr->DistFromEdge -> Fill (TrackDistFromEdge);
  HistPtr->XYPlaneStart -> Fill (larStart[0], larStart[1]);
  HistPtr->XZPlaneStart -> Fill (larStart[0], larStart[2]);
  HistPtr->YZPlaneStart -> Fill (larStart[1], larStart[2]);
  HistPtr->XYPlaneEnd   -> Fill (larEnd[0], larEnd[1]);
  HistPtr->XZPlaneEnd   -> Fill (larEnd[0], larEnd[2]);
  HistPtr->YZPlaneEnd   -> Fill (larEnd[1], larEnd[2]);
  HistPtr->PIDA_nHits   -> Fill (PIDA, TrackHits );
  return;
} // CheatPIDAFiller

// ******************************** Track Boundaries ****************************************************
void ProtonIdentification::ProtonIdentification::TrackBoundaries ( TVector3 larStart, TVector3 larEnd ) {
  
  double XStartdiff = std::min( abs(boundaries[0]-larStart[0]), abs(boundaries[1]-larStart[0]) );
  double YStartdiff = std::min( abs(boundaries[2]-larStart[1]), abs(boundaries[3]-larStart[1]) );
  double ZStartdiff = std::min( abs(boundaries[4]-larStart[2]), abs(boundaries[5]-larStart[2]) );
  if ( XStartdiff < fBoundaryEdge || YStartdiff < fBoundaryEdge || ZStartdiff < fBoundaryEdge )
    startsonboundary = true;
  else  startsonboundary = false;

  TrackDistFromEdge = std::min( XStartdiff, std::min(YStartdiff,ZStartdiff) );
  
    if (  std::min(    abs(boundaries[0]-larEnd[0]), abs(boundaries[1]-larEnd[0]) ) < fBoundaryEdge  
	  || std::min( abs(boundaries[2]-larEnd[1]), abs(boundaries[3]-larEnd[1]) ) < fBoundaryEdge 
	  || std::min( abs(boundaries[4]-larEnd[2]), abs(boundaries[5]-larEnd[2]) ) < fBoundaryEdge )
      endsonboundary = true;
    else  endsonboundary = false;
    
  // *** What to do if the track is reconstructed backwards though....
  // *** Do I want to look into that here or elsewhere?

  if(!startsonboundary && !endsonboundary ) RecoContainment = 0;  // contained track
  if(!startsonboundary &&  endsonboundary ) RecoContainment = 1;  // escaping track
  if( startsonboundary && !endsonboundary ) RecoContainment = 2;  // entering track
  if( startsonboundary &&  endsonboundary ) RecoContainment = 3;  // through going track
  
  trackstops = (RecoContainment == 0 || RecoContainment == 2);
 
} // TrackBoundaries
// *********************************** Monte Carlo Truth Extraction ********************************************************
void ProtonIdentification::ProtonIdentification::MCTruthInformation ( const simb::MCParticle *particle ) {
  int numberTrajectoryPoints = particle->NumberTrajectoryPoints(); // Looking at each MC hit
  double TPCLengthHits[numberTrajectoryPoints];
  bool BeenInVolume = false;
  int FirstHit=0, LastHit=0;
  MCTPCLength = 0;
  MCStartInTPC = MCEndInTPC = MCStops = false;

  MCPdgCode       = particle->PdgCode();
  MCTrackId       = particle->TrackId();
  
  for(int MCHit=0; MCHit < numberTrajectoryPoints; ++MCHit) {
    const TLorentzVector& tmpPosition=particle->Position(MCHit);
    double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
    
    if (MCHit!=0) TPCLengthHits[MCHit] = pow ( pow( (particle->Vx(MCHit-1)-particle->Vx(MCHit)),2)
					       + pow( (particle->Vy(MCHit-1)-particle->Vy(MCHit)),2)
					       + pow( (particle->Vz(MCHit-1)-particle->Vz(MCHit)),2)
					       , 0.5 );
    // --- Check if hit is in TPC...
    geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
    if (tpcid.isValid) { 
      if (MCHit == 0 ) MCStartInTPC = true;
      if (MCHit == numberTrajectoryPoints-1 ) MCEndInTPC = true;
      // -- Check if hit is within drift window...
      geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat);
      geo::TPCGeo      const& tpc  = cryo.TPC(tpcid.TPC); 
      double XPlanePosition      = tpc.PlaneLocation(0)[0];
      double DriftTimeCorrection = fabs( tmpPosition[0] - XPlanePosition ) / XDriftVelocity;
      double TimeAtPlane         = particle->T() + DriftTimeCorrection;
      if ( TimeAtPlane < detprop->TriggerOffset() 
	   || TimeAtPlane > detprop->TriggerOffset() + WindowSize 
	   ) continue;
      // -- Good hit in TPC
      LastHit = MCHit;
      if ( !BeenInVolume ) {
	BeenInVolume = true;
	FirstHit = MCHit;
      }
    } // TPC.valid	    
  } // MCTrajPoints
  if(!MCStartInTPC && !MCEndInTPC ) MCContainment = 0;  // through track
  if(!MCStartInTPC &&  MCEndInTPC ) MCContainment = 1;  // entering track
  if( MCStartInTPC && !MCEndInTPC ) MCContainment = 2;  // escaping track
  if( MCStartInTPC &&  MCEndInTPC ) MCContainment = 3;  // contained track
  
  MCStops = ( MCContainment == 0 || MCContainment == 2 );
  
  MCEnergy          = particle->E(FirstHit);
  MCEnergyDeposited = particle->E(FirstHit) - particle->E(LastHit);
  for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) {
    MCTPCLength += TPCLengthHits[Hit];
  }
  TLorentzVector& momentumStart  = (TLorentzVector&)particle->Momentum(FirstHit);   
  TVector3 mcstartmom = particle->Momentum(FirstHit).Vect();
  MCTheta_XZ = std::atan2(momentumStart.Px(), momentumStart.Pz());
  MCTheta_YZ = std::atan2(momentumStart.Py(), momentumStart.Pz());
  MCEta_XY   = std::atan2(momentumStart.Px(), momentumStart.Py());
  MCEta_ZY   = std::atan2(momentumStart.Pz(), momentumStart.Py());
  MCTheta    = mcstartmom.Theta();
  MCPhi      = mcstartmom.Phi();
  
  return;
} // MCTruthInformation
// ******************************** Reset Variables *****************************************************
void ProtonIdentification::ProtonIdentification::ResetVars() {
}
// ******************************** Define Module *****************************************************
DEFINE_ART_MODULE(ProtonIdentification::ProtonIdentification)
