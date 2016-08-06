////////////////////////////////////////////////////////////////////////
// Class:       NeutronDecayN2Ana
// Module Type: analyzer
// File:        NeutronDecayN2Ana_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
//
//  Having a look at FD reconstruction quantities
//
//  Thomas Karl Warburton
//  k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTracker.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TTree.h"

//standard library includes
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits> // std::numeric_limits<>

struct IDEYLess {
  bool operator()(const sim::IDE& first, const sim::IDE& second) {
    return first.y < second.y;
  }
};

#define MaxPart 100

namespace NeutronDecayN2Ana {
  class NeutronDecayN2Ana;
}

class NeutronDecayN2Ana::NeutronDecayN2Ana : public art::EDAnalyzer {
public:
  explicit NeutronDecayN2Ana(fhicl::ParameterSet const & p);
  virtual ~NeutronDecayN2Ana();

  void analyze(art::Event const & e) override;

  void beginRun(art::Run& run);
  void beginJob() override;
  void endJob();
  void endRun();
  //void reconfigure(fhicl::ParameterSet const & p) override;
  
private:
  // ------ My functions ------
  void ResetVars();
  void FillVars( std::vector<int> &TrackIDVec, int &numParts, double EDep[MaxPart], double Start[MaxPart][4], double End[MaxPart][4],
	    int ThisID, unsigned int ThisTDC, sim::IDE ThisIDE );
 
  // Handles
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTracker> bktrk;

  // Parameter List
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fMCTruthT0ModuleLabel;

  double ActiveBounds[6]; // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )

  TTree* fDecayTree;
  int Run;
  int Event;
  double PrimMuonRange;
  double PrimMuonEDep;
  double PrimMuonShadowEDep;
  double DistEdge[3][2];
  double TotalEDep;
  // Muon
  int nMuon;
  double MuonEDep[MaxPart], MuonStart[MaxPart][4], MuonEnd[MaxPart][4];
  // Charged Pion
  int nPion;
  double PionEDep[MaxPart], PionStart[MaxPart][4], PionEnd[MaxPart][4];
  // Pi0
  int nPi0;
  double Pi0EDep[MaxPart], Pi0Start[MaxPart][4], Pi0End[MaxPart][4];
  // K+
  int nKPlus;
  double KPlusEDep[MaxPart], KPlusStart[MaxPart][4], KPlusEnd[MaxPart][4];
  // K-
  int nKMinus;
  double KMinusEDep[MaxPart], KMinusStart[MaxPart][4], KMinusEnd[MaxPart][4];
  // e+
  int nEPlus;
  double EPlusEDep[MaxPart], EPlusStart[MaxPart][4], EPlusEnd[MaxPart][4];
  // e-
  int nEMinus;
  double EMinusEDep[MaxPart], EMinusStart[MaxPart][4], EMinusEnd[MaxPart][4];
 };
// ******************************** Reset Variables *****************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::ResetVars() {
  Run = Event = 0;
  PrimMuonRange = PrimMuonEDep = PrimMuonShadowEDep = TotalEDep = 0;
  // DistEdge
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      DistEdge[i][j]=99999999.;
  // Each particle type
  nMuon = nPion = nPi0 = nKPlus = nKMinus = nEPlus = nEMinus = 0;
  for (int i=0; i<MaxPart; ++i) {
    MuonEDep[i] = PionEDep[i] = Pi0EDep[i] = KPlusEDep[i] = KMinusEDep[i] = EPlusEDep[i] = EMinusEDep[i] = 0;
    for (int j=0; j<4; ++j) {
      MuonStart[i][j] = PionStart[i][j] = Pi0Start[i][j] = KPlusStart[i][j] = KMinusStart[i][j] = EPlusStart[i][j] = EMinusStart[i][j] = 0;
      MuonEnd[i][j]   = PionEnd[i][j]   = Pi0End[i][j]   = KPlusEnd[i][j]   = KMinusEnd[i][j]   = EPlusEnd[i][j]   = EMinusEnd[i][j]   = 0;
    }
  }
}
// ********************************** Begin Run *******************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::beginRun(art::Run& run) {

}
// *********************************** Begin Job ********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::beginJob()
{
  // Build my Cryostat boundaries array...Taken from Tyler Alion in Geometry Core. Should still return the same values for uBoone.
  ActiveBounds[0] = ActiveBounds[2] = ActiveBounds[4] = DBL_MAX;
  ActiveBounds[1] = ActiveBounds[3] = ActiveBounds[5] = -DBL_MAX;

  // ----- FixMe: Assume single cryostats ------
  auto const* geom = lar::providerFrom<geo::Geometry>();
  for (geo::TPCGeo const& TPC: geom->IterateTPCs()) {
    // get center in world coordinates
    const double origin[3] = {0.};
    double center[3] = {0.};
    TPC.LocalToWorld(origin, center);
    //double tpcDim[3] = {TPC.ActiveHalfWidth(), TPC.ActiveHalfHeight(), 0.5*TPC.ActiveLength() };
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() }; // Gives same result as what Matt has
    if( center[0] - tpcDim[0] < ActiveBounds[0] ) ActiveBounds[0] = center[0] - tpcDim[0];
    if( center[0] + tpcDim[0] > ActiveBounds[1] ) ActiveBounds[1] = center[0] + tpcDim[0];
    if( center[1] - tpcDim[1] < ActiveBounds[2] ) ActiveBounds[2] = center[1] - tpcDim[1];
    if( center[1] + tpcDim[1] > ActiveBounds[3] ) ActiveBounds[3] = center[1] + tpcDim[1];
    if( center[2] - tpcDim[2] < ActiveBounds[4] ) ActiveBounds[4] = center[2] - tpcDim[2];
    if( center[2] + tpcDim[2] > ActiveBounds[5] ) ActiveBounds[5] = center[2] + tpcDim[2];
  } // for all TPC
  std::cout << "Active Boundaries: "
	    << "\n\tx: " << ActiveBounds[0] << " to " << ActiveBounds[1]
	    << "\n\ty: " << ActiveBounds[2] << " to " << ActiveBounds[3]
	    << "\n\tz: " << ActiveBounds[4] << " to " << ActiveBounds[5]
	    << std::endl;
  /*
  double minx_, maxx_, miny_, maxy_, minz_, maxz_;
  minx_ = miny_ = minz_ = DBL_MAX;
  maxx_ = maxy_ = maxz_ = -DBL_MAX;
  for (unsigned int c=0; c<geom->Ncryostats(); c++)
    {
      const geo::CryostatGeo& cryostat=geom->Cryostat(c);
      for (unsigned int t=0; t<cryostat.NTPC(); t++)
	{
	  geo::TPCID id;
	  id.Cryostat=c;
	  id.TPC=t;
	  id.isValid=true;
	  const geo::TPCGeo& tpc=cryostat.TPC(id);
	  std::cout << t << "\t" << (tpc.Length()/2) << ", " << (tpc.ActiveLength()/2) << std::endl;
	  if (tpc.MinX()<minx_) minx_=tpc.MinX();
	  if (tpc.MaxX()>maxx_) maxx_=tpc.MaxX();
	  if (tpc.MinY()<miny_) miny_=tpc.MinY();
	  if (tpc.MaxY()>maxy_) maxy_=tpc.MaxY();
	  if (tpc.MinZ()<minz_) minz_=tpc.MinZ();
	  if (tpc.MaxZ()>maxz_) maxz_=tpc.MaxZ();
	  
	}
    }
  std::cout << "Matt's positions. " << minx_ << ", " << maxx_ << "\t" <<  miny_ << ", " << maxy_ << "\t" <<  minz_ << ", " << maxz_ << std::endl;
  */
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fDecayTree = tfs->make<TTree>("ReconstructedTree","analysis tree");
  fDecayTree->Branch("Run"          ,&Run          ,"Run/I"            );
  fDecayTree->Branch("Event"        ,&Event        ,"Event/I"          );
  fDecayTree->Branch("PrimMuonRange",&PrimMuonRange,"PrimMuonRange/D"  );
  fDecayTree->Branch("PrimMuonEDep" ,&PrimMuonEDep ,"PrimMuonEDep/D"   );
  fDecayTree->Branch("PrimMuonShadowEDep" ,&PrimMuonShadowEDep ,"PrimMuonShadowEDep/D"   );
  fDecayTree->Branch("DistEdge"     ,&DistEdge     ,"DistEdge[3][2]/D" );
  fDecayTree->Branch("TotalEDep"    ,&TotalEDep    ,"TotalEDep/D"      );
  // Muon
  fDecayTree->Branch("nMuon"     ,&nMuon      ,"nMuon/I"               );
  fDecayTree->Branch("MuonEDep"  ,&MuonEDep   ,"MuonEDep[nMuon]/D"     );
  fDecayTree->Branch("MuonStart" ,&MuonStart  ,"MuonStart[nMuon][4]/D" );
  fDecayTree->Branch("MuonEnd"   ,&MuonEnd    ,"MuonEnd[nMuon][4]/D"   );
  // Pion
  fDecayTree->Branch("nPion"     ,&nPion      ,"nPion/I"               );
  fDecayTree->Branch("PionEDep"  ,&PionEDep   ,"PionEDep[nPion]/D"     );
  fDecayTree->Branch("PionStart" ,&PionStart  ,"PionStart[nPion][4]/D" );
  fDecayTree->Branch("PionEnd"   ,&PionEnd    ,"PionEnd[nPion][4]/D"   );
  // Pi0
  fDecayTree->Branch("nPi0"     ,&nPi0      ,"nPi0/I"              );
  fDecayTree->Branch("Pi0EDep"  ,&Pi0EDep   ,"Pi0EDep[nPi0]/D"     );
  fDecayTree->Branch("Pi0Start" ,&Pi0Start  ,"Pi0Start[nPi0][4]/D" );
  fDecayTree->Branch("Pi0End"   ,&Pi0End    ,"Pi0End[nPi0][4]/D"   );
  // KPlus
  fDecayTree->Branch("nKPlus"     ,&nKPlus      ,"nKPlus/I"                );
  fDecayTree->Branch("KPlusEDep"  ,&KPlusEDep   ,"KPlusEDep[nKPlus]/D"     );
  fDecayTree->Branch("KPlusStart" ,&KPlusStart  ,"KPlusStart[nKPlus][4]/D" );
  fDecayTree->Branch("KPlusEnd"   ,&KPlusEnd    ,"KPlusEnd[nKPlus][4]/D"   );
  // KMinus
  fDecayTree->Branch("nKMinus"     ,&nKMinus      ,"nKMinus/I"                 );
  fDecayTree->Branch("KMinusEDep"  ,&KMinusEDep   ,"KMinusEDep[nKMinus]/D"     );
  fDecayTree->Branch("KMinusStart" ,&KMinusStart  ,"KMinusStart[nKMinus][4]/D" );
  fDecayTree->Branch("KMinusEnd"   ,&KMinusEnd    ,"KMinusEnd[nKMinus][4]/D"   );
  // EPlus
  fDecayTree->Branch("nEPlus"     ,&nEPlus      ,"nEPlus/I"                );
  fDecayTree->Branch("EPlusEDep"  ,&EPlusEDep   ,"EPlusEDep[nEPlus]/D"     );
  fDecayTree->Branch("EPlusStart" ,&EPlusStart  ,"EPlusStart[nEPlus][4]/D" );
  fDecayTree->Branch("EPlusEnd"   ,&EPlusEnd    ,"EPlusEnd[nEPlus][4]/D"   );
  // EMinus
  fDecayTree->Branch("nEMinus"     ,&nEMinus      ,"nEMinus/I"                 );
  fDecayTree->Branch("EMinusEDep"  ,&EMinusEDep   ,"EMinusEDep[nEMinus]/D"     );
  fDecayTree->Branch("EMinusStart" ,&EMinusStart  ,"EMinusStart[nEMinus][4]/D" );
  fDecayTree->Branch("EMinusEnd"   ,&EMinusEnd    ,"EMinusEnd[nEMinus][4]/D"   );
}
// ************************************ End Job *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::endJob() {
  
}
// ************************************ End Run *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::endRun() {
}

// ********************************** pset param *******************************************************
NeutronDecayN2Ana::NeutronDecayN2Ana::NeutronDecayN2Ana(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fHitsModuleLabel         ( pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel        ( pset.get< std::string >("TrackModuleLabel"))
  , fMCTruthT0ModuleLabel    ( pset.get< std::string >("MCTruthT0ModuleLabel"))
{

}
// ******************************************************************************************************
NeutronDecayN2Ana::NeutronDecayN2Ana::~NeutronDecayN2Ana()
{
  // Clean up dynamic memory and other resources here.

}
// ************************************ Analyse *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::analyze(art::Event const & evt)
{
  ResetVars();
  Run   = evt.run();
  Event = evt.event();
  std::cout << "\n\n************* New Event / Module running - Run " << Run << ", Event " << Event << " *************\n\n" << std::endl;

  // Any providers I need.
  auto const* geo = lar::providerFrom<geo::Geometry>();

  // Implementation of required member function here. 
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
  
  // Make a map of MCParticles which I can access later.
  art::Handle<std::vector<simb::MCParticle> > truth;
  evt.getByLabel("largeant", truth);
  std::map<int, const simb::MCParticle*> truthmap;
  for (size_t i=0; i<truth->size(); i++)
    truthmap[truth->at(i).TrackId()]=&((*truth)[i]);

  // Get a vector of sim channels.
  art::Handle<std::vector<sim::SimChannel> > simchannels;
  evt.getByLabel("largeant", simchannels);
  
  // Make a vector for primary IDEs to work out total Edep and length of primary muon
  std::vector<sim::IDE> priides;

  // Make vectors to hold all of my particle TrackIDs
  std::vector<int> MuonVec, PionVec, Pi0Vec, KPlusVec, KMinusVec, EPlusVec, EMinusVec;
  std::vector<int> AllTrackIDs;
  
  // Want to loop through the simchannels, and the ides.
  int chanIt = 0;
  for (auto const& simchannel:*simchannels) {
    // ------ Only want to look at collection plane hits ------
    if (geo->SignalType(simchannel.Channel()) != geo::kCollection) continue;
    int tdcideIt = 0;
    // ------ Loop through all the IDEs for this channel ------
    for (auto const& tdcide:simchannel.TDCIDEMap()) {
      unsigned int tdc = tdcide.first;
      auto const& idevec=tdcide.second;
      //std::cout << "TDC IDE " << tdcideIt << " of " << simchannel.TDCIDEMap().size() << ", has tdc " << tdc << ", and idevec of size " << idevec.size() << std::endl;
      int ideIt = 0;
      // ------ Look at each individual IDE ------
      for (auto const& ide:idevec) {
	int ideTrackID = ide.trackID;
	float ideEnergy = ide.energy;
	double idePos[3];
	idePos[0] = ide.x;
	idePos[1] = ide.y;
	idePos[2] = ide.z;
	/*
	if (ideTrackID == 1) { 
	float ideNumEl  = ide.numElectrons;
	std::cout << "      IDE " << ideIt << " of " << idevec.size() << " has TrackId " << ideTrackID << ", energy " << ideEnergy << ", NumEl " << ideNumEl << ", tdc " << tdc << ". ";
	std::cout << "Position " << idePos[0] << ", " << idePos[1] << ", " << idePos[2] << std::endl;
	}
	//*/

	// ------ Work out if this hit is in a TPC ------
	bool InTPC = false;
	if ( idePos[0] > ActiveBounds[0] && idePos[0] < ActiveBounds[1] )
	  if ( idePos[1] > ActiveBounds[2] && idePos[1] < ActiveBounds[3] )
	    if ( idePos[2] > ActiveBounds[4] && idePos[2] < ActiveBounds[5] )
	      InTPC = true;
	if (! InTPC ) {
	  std::cout << "Outside the Active volume I found at the top!" << std::endl;
	  continue;
	}	  

	// ------ Add the energy deposition from this IDE to the sum of IDEs
	TotalEDep += ideEnergy;
	if ( ideTrackID == -1 ) PrimMuonShadowEDep += ideEnergy;
	
	// ------ I want to make a vector of the IDEs for the primary muon ------
	if (ideTrackID==1) priides.push_back(ide);
	  
	// ------ I want to work the closest IDE to an edge of the active volume ------
	if ( DistEdge[0][0] > ( ide.x - ActiveBounds[0]) ) DistEdge[0][0] =  ide.x - ActiveBounds[0];
	if ( DistEdge[0][1] > (-ide.x + ActiveBounds[1]) ) DistEdge[0][1] = -ide.x + ActiveBounds[1];
	if ( DistEdge[1][0] > ( ide.y - ActiveBounds[2]) ) DistEdge[1][0] =  ide.y - ActiveBounds[2];
	if ( DistEdge[1][1] > (-ide.y + ActiveBounds[3]) ) DistEdge[1][1] = -ide.y + ActiveBounds[3];
	if ( DistEdge[2][0] > ( ide.z - ActiveBounds[4]) ) DistEdge[2][0] =  ide.z - ActiveBounds[4];
	if ( DistEdge[2][1] > (-ide.z + ActiveBounds[5]) ) DistEdge[2][1] = -ide.z + ActiveBounds[5];
	/*
	  std::cout << "DistEdge[0][0] = " <<  ide.x - ActiveBounds[0] << " = " <<  ide.x << " - " << ActiveBounds[0] << ". The minimum is " << DistEdge[0][0] << std::endl;
	  std::cout << "DistEdge[0][1] = " << -ide.x + ActiveBounds[1] << " = " << -ide.x << " + " << ActiveBounds[1] << ". The minimum is " << DistEdge[0][1] << std::endl;
	  std::cout << "DistEdge[1][0] = " <<  ide.y - ActiveBounds[2] << " = " <<  ide.y << " - " << ActiveBounds[2] << ". The minimum is " << DistEdge[1][0] << std::endl;
	  std::cout << "DistEdge[1][1] = " << -ide.y + ActiveBounds[3] << " = " << -ide.y << " + " << ActiveBounds[3] << ". The minimum is " << DistEdge[1][1] << std::endl;
	  std::cout << "DistEdge[2][0] = " <<  ide.z - ActiveBounds[4] << " = " <<  ide.z << " - " << ActiveBounds[4] << ". The minimum is " << DistEdge[2][0] << std::endl;
	  std::cout << "DistEdge[2][1] = " << -ide.z + ActiveBounds[5] << " = " << -ide.z << " + " << ActiveBounds[5] << ". The minimum is " << DistEdge[2][1] << std::endl;
	//*/
	
	// ------ I can't get the truth information for negative track IDs ------
	if ( ideTrackID < 0 ) continue;

	// ------ Now to work out which particles in particular I want to save more information about... ------
	const simb::MCParticle& part=*(truthmap[ideTrackID]);
	int PdgCode=part.PdgCode();
	
	if      ( PdgCode == -13  || PdgCode == 13  ) FillVars( MuonVec  , nMuon  , MuonEDep  , MuonStart  , MuonEnd  , ideTrackID, tdc, ide ); // Muon
	else if ( PdgCode == -211 || PdgCode == 211 ) FillVars( PionVec  , nPion  , PionEDep  , PionStart  , PionEnd  , ideTrackID, tdc, ide ); // Pion
	else if ( PdgCode == 111  )                   FillVars( Pi0Vec   , nPi0   , Pi0EDep   , Pi0Start   , Pi0End   , ideTrackID, tdc, ide ); // Pi0
	else if ( PdgCode == 321  )                   FillVars( KPlusVec , nKPlus , KPlusEDep , KPlusStart , KPlusEnd , ideTrackID, tdc, ide ); // KPlus
	else if ( PdgCode == -321 )                   FillVars( KMinusVec, nKMinus, KMinusEDep, KMinusStart, KMinusEnd, ideTrackID, tdc, ide ); // KMinus
	else if ( PdgCode == -11  )                   FillVars( EPlusVec , nEPlus , EPlusEDep , EPlusStart , EPlusEnd , ideTrackID, tdc, ide ); // EPlus
	else if ( PdgCode == 11   )                   FillVars( EMinusVec, nEMinus, EMinusEDep, EMinusStart, EMinusEnd, ideTrackID, tdc, ide ); // EMinus
	
	std::vector<int>::iterator it=std::find(AllTrackIDs.begin(), AllTrackIDs.end(), ideTrackID);
	if ( it==AllTrackIDs.end() ) {
	  AllTrackIDs.push_back( ideTrackID );
	  std::cerr << "I also had trackID " << ideTrackID << " in this event, it was from a " << PdgCode << std::endl;
	}

	ideIt++;
      } // Each IDE ( ide:idevec )
      ++tdcideIt;
    } // IDE vector for SimChannel ( tdcide:simcahnnel.TPCIDEMap() )
    ++chanIt;
  } // Loop through simchannels

  // Work out some properties of the primary muon
  if (priides.size()) {
    // ------ Sort the IDEs in y ------
    std::sort(priides.begin(), priides.end(), IDEYLess());
    // ------ Work out the range ------
    double dx=priides.front().x-priides.back().x;
    double dy=priides.front().y-priides.back().y;
    double dz=priides.front().z-priides.back().z;
    PrimMuonRange=sqrt(dx*dx+dy*dy+dz*dz);
    // ------ Work out the energy deposited ------
    for (unsigned int muIt=0; muIt<priides.size(); ++muIt) PrimMuonEDep += priides[muIt].energy;
    PrimMuonEDep = PrimMuonEDep / priides.size();
  } 
  std::cout << "Primary muon track length = " << PrimMuonRange << " cm, and energy despoited = " << PrimMuonEDep << ", shadow EDep = " << PrimMuonShadowEDep << " and total Edep is " << TotalEDep << ".\n"
	    << "There were " << nMuon << " muons, " << nPion << " pions, " << nPi0 << " pi0s, " << nKPlus << " KPlus, " << nKMinus << " KMinus, " << nEPlus << " EPlus, " << nEMinus << " EMinus.\n"
	    << std::endl;
  // If nMuon
  if (nMuon) {
    double MyRange = pow( pow((MuonStart[0][0]-MuonEnd[0][0]),2) + pow((MuonStart[0][1]-MuonEnd[0][1]),2) + pow((MuonStart[0][2]-MuonEnd[0][2]),2), 0.5 );
    std::cout << "PrimMuon start " << MuonStart[0][0] << ", " << MuonStart[0][1] << ", " << MuonStart[0][2] << ", " << MuonStart[0][3]
	      << ". PrimMuon end " << MuonEnd[0][0] << ", " << MuonEnd[0][1] << ", " << MuonEnd[0][2] << ", " << MuonEnd[0][3]
	      << ". Distance of " << pow( pow((MuonStart[0][0]-MuonEnd[0][0]),2) + pow((MuonStart[0][1]-MuonEnd[0][1]),2) + pow((MuonStart[0][2]-MuonEnd[0][2]),2), 0.5 )
	      << std::endl;
    
    if ( PrimMuonRange - MyRange > 2 ) {
      std::cout << "\n NOT MATCHING!! \n" << std::endl;
      for (size_t qq=0; qq<priides.size(); ++qq)
	std::cout << "PrimMuon start " << priides[qq].x << ", " << priides[qq].y << ", " << priides[qq].z << std::endl;
    }
  }  
  // ------ Fill the Tree ------
  fDecayTree->Fill();
} // Analyse
// ******************************** Fill variables *****************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::FillVars( std::vector<int> &TrackIDVec, int &numParts, double EDep[MaxPart], double Start[MaxPart][4], double End[MaxPart][4],
						     int ThisID, unsigned int ThisTDC, sim::IDE ThisIDE ) {

  std::vector<int>::iterator it=std::find(TrackIDVec.begin(), TrackIDVec.end(), ThisID);
  int partNum;
  if ( it==TrackIDVec.end() ) {
    partNum = numParts;
    ++numParts;
    TrackIDVec.push_back( ThisID );
    std::cerr << "Pushing back a new ideTrackID " << ThisID << std::endl;
    Start[partNum][0] = End[partNum][0] = ThisIDE.x;
    Start[partNum][1] = End[partNum][1] = ThisIDE.y;
    Start[partNum][2] = End[partNum][2] = ThisIDE.z;
    Start[partNum][3] = End[partNum][3] = ThisTDC;
  } else {
    partNum = it - TrackIDVec.begin();
    if ( ThisTDC < Start[partNum][3] ) {
      Start[partNum][0] = ThisIDE.x;
      Start[partNum][1] = ThisIDE.y;
      Start[partNum][2] = ThisIDE.z;
      Start[partNum][3] = ThisTDC;
    } else if ( ThisTDC > End[partNum][3] ) {
      End[partNum][0] = ThisIDE.x;
      End[partNum][1] = ThisIDE.y;
      End[partNum][2] = ThisIDE.z;
      End[partNum][3] = ThisTDC;
    }
  }
  EDep[partNum] += ThisIDE.energy;
} // FillVars
// ******************************** Define Module *****************************************************
DEFINE_ART_MODULE(NeutronDecayN2Ana::NeutronDecayN2Ana)
