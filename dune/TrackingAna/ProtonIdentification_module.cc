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
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/OpFlash.h"
#include "DetectorInfoServices/DetectorPropertiesService.h"
#include "Utilities/AssociationUtil.h"
#include "DetectorInfoServices/DetectorClocksService.h"
#include "RawData/ExternalTrigger.h"
#include "MCCheater/BackTracker.h"
#include "AnalysisBase/Calorimetry.h"
#include "AnalysisBase/T0.h"
#include "AnalysisBase/ParticleID.h"

#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"

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
  } Proton, Muon;
  struct TEfficiencies {
    //TEfficiencies();
    //TEfficiencies(const std::string& subdir);
    TEfficiency *Effic_Length;
    TEfficiency *Effic_Energy;
    TEfficiency *Effic_Energy_Depos;
    TEfficiency *Effic_Theta;
    TEfficiency *Effic_Theta_XZ;
    TEfficiency *Effic_Theta_YZ;
    TEfficiency *Effic_Phi;
    TEfficiency *Effic_Phi_v_Theta;
    TEfficiency *Effic_Eta_XY;
    TEfficiency *Effic_Eta_ZY;
  } AllEfficiencies, ProtonEfficiencies;
  // ----------- Declare my structs ----------

private:
  
  void ResetVars();
  void MCTruthInformation ( const simb::MCParticle *particle, double &Energy, double &EnergyDeposited, double &TPCLength, 
			    double &Theta_XZ, double &Theta_YZ, double &Eta_XY, double &Eta_ZY, double &Theta, double &Phi,
			    int &MCPdgCode, int &MCTrackId, int &StartInTPC, int &EndInTPC, int &Contained );
  void CheatPIDAFiller ( struct ProtonIdentification::ProtonIdentification::CheatHists *HistPtr, 
			 std::vector<const anab::Calorimetry*> calos, double &PIDA );
  void Make_Efficiencies ( struct ProtonIdentification::ProtonIdentification::TEfficiencies *AllHistPtr, bool matched, 
			   double Energy, double EnergyDeposited, double TPCLength, 
			   double Theta_XZ, double Theta_YZ, double Eta_XY, double Eta_ZY, double Theta, double Phi );
  void EfficPlot_1D ( TH1D *AllHist, TH1D *MatchHist, TEfficiency *EfficHist );
  void EfficPlot_2D ( TH2D *AllHist, TH2D *MatchHist, TEfficiency *EfficHist );

  TTree* fTree;

  double MCTruthT0, MCTruthTickT0, MCTruthTrackID;
  double PhotonCounterT0, PhotonCounterTickT0, PhotonCounterID;
  double TrackTheta_XZ, TrackTheta_YZ, TrackEta_XY, TrackEta_ZY, TrackTheta, TrackPhi;
  double TrackLength; 
  std::vector<double> trackStart, trackEnd;

  double MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi;
  double MCTPCLength, MCEnergy, MCEnergyDeposited;
  double PIDA;
  int MCPdgCode, MCTrackId, StartInTPC, EndInTPC, Contained;
  int MatchedTrackID;
  
  bool AllChargedTrack;
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
  std::string fCalorimetryModuleLabel;

  double PIDApower;
 };
// ********************************** Begin Run *******************************************************
void ProtonIdentification::ProtonIdentification::beginRun(art::Run& run) {

}
// *********************************** Begin Job ********************************************************
void ProtonIdentification::ProtonIdentification::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("CheatedTree","analysis tree");
  fTree->Branch("MCTPCLength"      ,&MCTPCLength      ,"MCTPCLength/D"      );
  fTree->Branch("MCEnergy"         ,&MCEnergy         ,"MCEnergy/D"         );
  fTree->Branch("MCEnergyDeposited",&MCEnergyDeposited,"MCEnergyDeposited/D");
  fTree->Branch("MCTheta"          ,&MCTheta          ,"MCTheta/D"          );
  fTree->Branch("MCTheta_XZ"       ,&MCTheta_XZ       ,"MCTheta_XZ/D"       );
  fTree->Branch("MCTheta_YZ"       ,&MCTheta_YZ       ,"MCTheta/D"          );
  fTree->Branch("MCPhi"            ,&MCPhi            ,"MCPhi/D"            );
  fTree->Branch("MCEta_XY"         ,&MCEta_XY         ,"MCEta_XY/D"         );
  fTree->Branch("MCEta_ZY"         ,&MCEta_ZY         ,"MCEta_ZY/D"         );
  fTree->Branch("MCPdgCode"        ,&MCPdgCode        ,"MCPdgCode/I"        );
  fTree->Branch("MCTrackId"        ,&MCTrackId        ,"MCTrackId/I"        );
  fTree->Branch("MatchedTrackID"   ,&MatchedTrackID   ,"MatchedTrackID/I"   );
  fTree->Branch("MatchedTrackLength",&TrackLength    ,"MatchedTrackLength/D");
  fTree->Branch("PIDA"             ,&PIDA             ,"PIDA/D"             );

  // ----- CHEATED HISTOGRAMS ------- Protons first -----
  Proton.dEdx     = tfs->make<TH1D>("Proton_dEdx", "dEdx values for protons in the detector; dEdx (MeV cm); Number", 100, 0, 50);
  Proton.ResRange = tfs->make<TH1D>("Proton_ResRange", "Residual range for protons in the detector; Residual Range (cm); Number", 100, 0, 50);
  Proton.dEdx_ResRange = tfs->make<TH2D>("Proton_dEdx_ResRange", 
					 "Residual range against dEdx for protons in the detector; Residual Range (cm); dEdx (MeV cm)",
					 100, 0, 50, 100, 0, 50);
  Proton.PIDA     = tfs->make<TH1D>("Proton_PIDA", "PIDA values for protons in the detector; PIDA, Number", 100, 0, 30);
  // ------ Now for Muons ------
  Muon.dEdx     = tfs->make<TH1D>("Muon_dEdx", "dEdx values for muons in the detector; dEdx (MeV cm); Number", 100, 0, 50);
  Muon.ResRange = tfs->make<TH1D>("Muon_ResRange", "Residual range for muons in the detector; Residual Range (cm); Number", 100, 0, 50);
  Muon.dEdx_ResRange = tfs->make<TH2D>("Muon_dEdx_ResRange", 
				       "Residual range against dEdx for muons in the detector; Residual Range (cm); dEdx (MeV cm)",
				       100, 0, 50, 100, 0, 50);
  Muon.PIDA     = tfs->make<TH1D>("Muon_PIDA", "PIDA values for muons in the detector; PIDA, Number", 100, 0, 30);

  // ----- Non-Cheated Histograms ------- From here on in --------
}
// ************************************ End Job *********************************************************
void ProtonIdentification::ProtonIdentification::endJob() {
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
  , fCalorimetryModuleLabel  ( pset.get< std::string >("CalorimetryModuleLabel"))
  , PIDApower                ( pset.get< double      >("PIDApower"))
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
  
  const sim::ParticleList& plist = bktrk->ParticleList();
  
  art::Handle< std::vector<recob::Track> > trackh;
  evt.getByLabel(fTrackModuleLabel, trackh);
  
  art::Handle< std::vector< art::PtrVector < recob::Track > > > trackvh;
  evt.getByLabel(fTrackModuleLabel, trackvh);

  int NPart = 0;

  std::cout << "Looking at a new event / module in my proton identification module!!" << std::endl;

  if ( trackListHandle.isValid() ) { // Check that trackListHandle is valid.....
    art::FindManyP<recob::Hit>        fmht  (trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::T0>           fmt0  (trackListHandle, evt, fMCTruthT0ModuleLabel);
    art::FindMany<anab::T0>           fmphot(trackListHandle, evt, fPhotonT0ModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal (trackListHandle, evt, fCalorimetryModuleLabel);
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
	//TickT0 = PhotonCounterTickT0;
      } 
      // ************** END T0 stuff ***************
      
      // ******************************************************************************************
      // Correct X and get track length etc now that we have matched a Track with an MCParticle!!
      // ******************************************************************************************
      // ----- Hit Stuff ------
      std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(Track);
      double Hit_Size = allHits.size();
      // ---- Correct X positions and get track positions.
      trackStart.clear();
      trackEnd.clear();
      tracklist[Track]->Extent(trackStart,trackEnd);
      trackStart[0] = trackStart[0] - detprop->ConvertTicksToX( TickT0, allHits[Hit_Size-1]->WireID().Plane, allHits[Hit_Size-1]->WireID().TPC, allHits[Hit_Size-1]->WireID().Cryostat ); // Correct X, last entry is first 'hit'
      trackEnd[0]   = trackEnd[0] - detprop->ConvertTicksToX( TickT0, allHits[0]->WireID().Plane, allHits[0]->WireID().TPC, allHits[0]->WireID().Cryostat ); // Correct X, first entry is last 'hit'
      // ---- Get lengths and angles.
      art::Ptr<recob::Track> ptrack(trackh, Track);
      const recob::Track& track = *ptrack;
      TrackLength = track.Length();
      int ntraj = track.NumberTrajectoryPoints();
      if(ntraj > 0) {
	TVector3 dir = track.VertexDirection();
	TrackTheta_XZ = std::atan2(dir.X(), dir.Z());
	TrackTheta_YZ = std::atan2(dir.Y(), dir.Z());
	TrackEta_XY   = std::atan2(dir.X(), dir.Y());
	TrackEta_ZY   = std::atan2(dir.Z(), dir.Y());
	TrackTheta    = dir.Theta();
	TrackPhi      = dir.Phi();
      }
      
      // *****************************************
      // IF CHEATING now want to loop through MC
      // *****************************************
      for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
	simb::MCParticle *particle = ipar->second;
	
	// ***** Check whether MCParticle made this track *****
	if ( particle->TrackId() != MCTruthTrackID ) continue;
	MatchedTrackID = Track;
	
	AllChargedTrack = false; // MCParticle not charged by default
	  
	// ----- Check that have a charged particle -----
	if ( fabs(particle->PdgCode() ) == 13 )  AllChargedTrack = true;
	else if ( particle->PdgCode() == 2212 ) AllChargedTrack = true;
	/*
	  else if ( fabs(particle->PdgCode() == 11) ) AllChargedTrack = true;  // Electron
	  else if ( fabs(particle->PdgCode() == 211) ) AllChargedTrack = true; // Pions
	  else if ( fabs(particle->PdgCode() == 321) ) AllChargedTrack = true; // Kaons
	//*/
	// ******** ONLY WANT CHARGED PARTICLES - MUONS AND PROTONS BY DEFAULT *********
	if (!AllChargedTrack) continue;

	// ---- Get MCTruth Information and check that MCParticle goes in TPC ---
	MCTruthInformation ( particle, MCEnergy, MCEnergyDeposited, MCTPCLength, 
			     MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi,
			     MCPdgCode, MCTrackId, StartInTPC, EndInTPC, Contained );
	
	// ******** Select only particles with some properties we are interested in *********
	if ( MCTPCLength == 0 ) continue; // Only want particles which enter the TPC
	if ( EndInTPC    == 0 ) continue; // Only want particles which stop in the TPC, so can use residual range! 
	
	++NPart;
	//std::cout << "Looking at MCParticle " << MCTrackId << ", with PdgCode " << MCPdgCode << ", MCTPCLength " << MCTPCLength << std::endl;
	
	// *********** Corrected X and track information **************
	// ****** Call PIDA cheat for this track/particle combination *************
	// ---- Make a vector of calorimetry objects...
	std::vector<const anab::Calorimetry*> calos = fmcal.at(Track);
	if ( fabs(MCPdgCode) == 13 ) CheatPIDAFiller( &Muon, calos, PIDA );
	else if ( MCPdgCode == 2212     ) {
	  if ( particle -> NumberDaughters() != 0 ) continue;
	  CheatPIDAFiller( &Proton, calos, PIDA );
	}
	// ****** Calculated PIDA for this track/particle combination *************
	// ******** Fill Tree for each MCParticle **********
	fTree->Fill();
      } // Loop through MCParticles....
      
      // ***********************************************************
      // IF NOT CHEATIING now want to work out PIDA with NO MC INFO
      // ***********************************************************
      
    } // Loop over Tracks    
  } // if trackListHandle.isValid()
  //std::cout << "\nThis event had " << NPart << " particles that were 'interesting'." << std::endl;
}
// *********************************** Monte Carlo Truth Extraction ********************************************************
void ProtonIdentification::ProtonIdentification::MCTruthInformation ( const simb::MCParticle *particle, double &Energy, double &EnergyDeposited, double &TPCLength, 
								      double &Theta_XZ, double &Theta_YZ, double &Eta_XY, double &Eta_ZY, double &Theta, double &Phi,
								      int &MCPdgCode, int &MCTrackId, int &StartInTPC, int &EndInTPC, int &Contained ) {
  int numberTrajectoryPoints = particle->NumberTrajectoryPoints(); // Looking at each MC hit
  double TPCLengthHits[numberTrajectoryPoints];
  bool BeenInVolume = false;
  int FirstHit=0, LastHit=0;
  TPCLength = StartInTPC = EndInTPC = Contained = 0;

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
      if (MCHit == 0 ) StartInTPC = 1;
      if (MCHit == numberTrajectoryPoints-1 ) EndInTPC = 1;
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
  if ( StartInTPC == 1 && EndInTPC == 1 ) Contained = 1;
  
  Energy          = particle->E(FirstHit);
  EnergyDeposited = particle->E(FirstHit) - particle->E(LastHit);
  for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) {
    TPCLength += TPCLengthHits[Hit];
  }
  TLorentzVector& momentumStart  = (TLorentzVector&)particle->Momentum(FirstHit);   
  TVector3 mcstartmom = particle->Momentum(FirstHit).Vect();
  Theta_XZ = std::atan2(momentumStart.Px(), momentumStart.Pz());
  Theta_YZ = std::atan2(momentumStart.Py(), momentumStart.Pz());
  Eta_XY   = std::atan2(momentumStart.Px(), momentumStart.Py());
  Eta_ZY   = std::atan2(momentumStart.Pz(), momentumStart.Py());
  Theta    = mcstartmom.Theta();
  Phi      = mcstartmom.Phi();
  
  return;
} // MCTruthInformation

// ******************************** Efficiency Filler ****************************************************
void ProtonIdentification::ProtonIdentification::CheatPIDAFiller ( struct ProtonIdentification::ProtonIdentification::CheatHists *HistPtr, 
								   std::vector<const anab::Calorimetry*> calos, double &PIDA ) {
  int BestPlaneHits, UsedHits;
  PIDA = BestPlaneHits = UsedHits = 0;
  for ( int Plane=0; Plane < (int)calos.size(); ++Plane ) { // Loop through planes
    if ( (int)calos[Plane]->dEdx().size() > BestPlaneHits ) // Work out which plane has the most hits
      BestPlaneHits = (int)calos[Plane]->dEdx().size();
    for ( int PlaneHit=0; PlaneHit < (int)calos[Plane]->dEdx().size(); ++PlaneHit ) { // loop through hits on each plane
      if ( calos[Plane]->ResidualRange()[PlaneHit] < 30 ) { // Only want PIDA for last 30 cm
	PIDA += calos[Plane]->dEdx()[PlaneHit] * pow(calos[Plane]->ResidualRange()[PlaneHit], PIDApower ); 
	++UsedHits;
      } // If ResRange < 30 cm
      HistPtr->dEdx          -> Fill( calos[Plane]->dEdx()[PlaneHit]          );
      HistPtr->ResRange      -> Fill( calos[Plane]->ResidualRange()[PlaneHit] );
      HistPtr->dEdx_ResRange -> Fill( calos[Plane]->ResidualRange()[PlaneHit], calos[Plane]->dEdx()[PlaneHit] );
    } // Loop over hits on each plane
  } // Loop over planes
  
  if ( UsedHits != 0 ) { // If had any hits, work out PIDA and calculate
    PIDA = PIDA / UsedHits;
    HistPtr->PIDA -> Fill( PIDA );
  }
  return;
} // HistoFiller

// ******************************** Make Efficiencies ****************************************************
void ProtonIdentification::ProtonIdentification::Make_Efficiencies ( struct ProtonIdentification::ProtonIdentification::TEfficiencies *AllHistPtr, bool matched, 
								 double Energy, double EnergyDeposited, double TPCLength, 
								 double Theta_XZ, double Theta_YZ, double Eta_XY, double Eta_ZY, double Theta, double Phi ) {
  AllHistPtr->Effic_Length      -> Fill ( matched, TPCLength       );
  AllHistPtr->Effic_Energy      -> Fill ( matched, Energy          );
  AllHistPtr->Effic_Energy_Depos-> Fill ( matched, EnergyDeposited );
  AllHistPtr->Effic_Theta       -> Fill ( matched, Theta           );
  AllHistPtr->Effic_Theta_XZ    -> Fill ( matched, Theta_XZ        );
  AllHistPtr->Effic_Theta_YZ    -> Fill ( matched, Theta_YZ        );
  AllHistPtr->Effic_Phi         -> Fill ( matched, Phi             );
  AllHistPtr->Effic_Eta_XY      -> Fill ( matched, Eta_XY          );
  AllHistPtr->Effic_Eta_ZY      -> Fill ( matched, Eta_ZY          );

  AllHistPtr->Effic_Phi_v_Theta -> Fill ( matched, Theta, Phi      );
}
// ******************************** Reset Variables *****************************************************
void ProtonIdentification::ProtonIdentification::ResetVars() {
}
// ******************************** Define Module *****************************************************
DEFINE_ART_MODULE(ProtonIdentification::ProtonIdentification)
