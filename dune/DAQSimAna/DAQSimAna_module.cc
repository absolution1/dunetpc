////////////////////////////////////////////////////////////////////////
// Class:       DAQSimAna
// Module Type: analyzer
// File:        DAQSimAna_module.cc
//
// Generated by Michael Baird using the old copy and paste...
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

// C++ includes

// ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

// Framework includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//const int nMaxDigs = 4492; // unused

enum PType{ kUnknown=0, kMarl, kAPA, kCPA, kAr39, kAr42, kNeutron, kKryp, kPlon, kRdon, kNPTypes };

class DAQSimAna : public art::EDAnalyzer {

public:

  explicit DAQSimAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  DAQSimAna(DAQSimAna const &) = delete;
  DAQSimAna(DAQSimAna &&) = delete;
  DAQSimAna & operator = (DAQSimAna const &) = delete;
  DAQSimAna & operator = (DAQSimAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & evt) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;

private:

  // --- Some of our own functions.
  void ResetVariables();
    void  FillMyMaps  ( std::map< int, simb::MCParticle> &MyMap,
                        art::FindManyP<simb::MCParticle> Assn,
                        art::ValidHandle< std::vector<simb::MCTruth> > Hand,
                        std::map<int, int>* indexMap=nullptr);
  PType WhichParType( int TrID );
  PType WhichParType( const art::ValidHandle<simb::MCTruth>& truthHand );
  bool  InMyMap     ( int TrID, std::map< int, simb::MCParticle> ParMap );
  void  CalcAdjHits ( std::vector< recob::Hit > MyVec, TH1I* MyHist, bool HeavDebug="false" );
  void SaveIDEs(art::Event const & evt);

  // --- Our fcl parameter labels for the modules that made the data products
  std::string fRawDigitLabel;
  std::string fHitLabel;
  bool fDoCalcAdjHits;

  std::string fGEANTLabel;
  std::string fMARLLabel; std::map< int, simb::MCParticle > MarlParts;
  std::string fAPALabel;  std::map< int, simb::MCParticle > APAParts;
  std::string fCPALabel;  std::map< int, simb::MCParticle > CPAParts;
  std::string fAr39Label; std::map< int, simb::MCParticle > Ar39Parts;
  std::string fAr42Label; std::map< int, simb::MCParticle > Ar42Parts;
  std::string fNeutLabel; std::map< int, simb::MCParticle > NeutParts;
  std::string fKrypLabel; std::map< int, simb::MCParticle > KrypParts;
  std::string fPlonLabel; std::map< int, simb::MCParticle > PlonParts;
  std::string fRdonLabel; std::map< int, simb::MCParticle > RdonParts;

  std::map<int, int> trkIDToMarleyIndex;

  // Mapping from track ID to particle type, for use in WhichParType()
  std::map<int, PType> trkIDToPType;

  // --- Other variables
  //int nADC; // no longer used

  // --- Our TTree, and its associated variables.
  TTree* fDAQSimTree;
  // General event info.
  int Run;
  int SubRun;
  int Event;
  // Raw digits
  //int NTotDigs; // unused

  // The reconstructed hits
  int   NTotHits;
  int   NColHits;
  int   NIndHits;

  std::vector<int>   HitView; ///< View i.e Coll, U, V
  std::vector<int>   HitSize; ///< Time width (ticks) Start - End time
  std::vector<int>   HitTPC; ///< The TPC which the hit occurs in
  std::vector<int>   HitChan; ///< The channel which the hit occurs on
  std::vector<float> HitTime; ///< The time of the hit (ticks)
  std::vector<float> HitRMS; ///< The RMS of the hit
  std::vector<float> HitSADC; ///< The summed ADC of the hit
  std::vector<float> HitInt; ///< The ADC integral of the hit
  std::vector<float> HitPeak; ///< The peak ADC value of the hit
  std::vector<int>   GenType; ///< The generator which generated the particle responsible for the hit
  std::vector<int>   MarleyIndex; ///< Which SN in the list of Marley interactions this hit is from (-1 if not from SN)

  int   TotGen_Marl;
  int   TotGen_APA;
  int   TotGen_CPA;
  int   TotGen_Ar39;
  int   TotGen_Ar42;
  int   TotGen_Neut;
  int   TotGen_Kryp;
  int   TotGen_Plon;
  int   TotGen_Rdon;

  int   NTotIDEs;
  std::vector<int>   IDEChannel;
  std::vector<int>   IDEStartTime;
  std::vector<int>   IDEEndTime;
  std::vector<float> IDEEnergy;
  std::vector<float> IDEElectrons;
  std::vector<int>   IDEParticle;

  // histograms to fill about Collection plane hits
  TH1I* hAdjHits_Marl;
  TH1I* hAdjHits_APA;
  TH1I* hAdjHits_CPA;
  TH1I* hAdjHits_Ar39;
  TH1I* hAdjHits_Ar42;
  TH1I* hAdjHits_Neut;
  TH1I* hAdjHits_Kryp;
  TH1I* hAdjHits_Plon;
  TH1I* hAdjHits_Rdon;
  TH1I* hAdjHits_Oth;

  // --- Declare our services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  
};

//......................................................
DAQSimAna::DAQSimAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}

//......................................................
void DAQSimAna::reconfigure(fhicl::ParameterSet const & p)
{
  fRawDigitLabel = p.get<std::string> ("RawDigitLabel");  
  fHitLabel      = p.get<std::string> ("HitLabel");

  fGEANTLabel = p.get<std::string> ("GEANT4Label");
  fMARLLabel = p.get<std::string> ("MARLEYLabel");
  fAPALabel  = p.get<std::string> ("APALabel");
  fCPALabel  = p.get<std::string> ("CPALabel");
  fAr39Label = p.get<std::string> ("Argon39Label");
  fAr42Label = p.get<std::string> ("Argon42Label");
  fNeutLabel = p.get<std::string> ("NeutronLabel");
  fKrypLabel = p.get<std::string> ("KryptonLabel");
  fPlonLabel = p.get<std::string> ("PoloniumLabel");
  fRdonLabel = p.get<std::string> ("RadonLabel");

  fDoCalcAdjHits = p.get<bool>("DoCalcAdjHits", false);
} // Reconfigure

//......................................................
void DAQSimAna::ResetVariables()
{
  // Clear my MCParticle maps.
  MarlParts.clear(); APAParts .clear(); CPAParts .clear(); Ar39Parts.clear();
  Ar42Parts.clear();
  NeutParts.clear(); KrypParts.clear(); PlonParts.clear(); RdonParts.clear();

  trkIDToPType.clear();

  // General event info.
  Run = SubRun = Event = -1;
  
  // Set Number of GenParts to 0
  TotGen_Marl = TotGen_APA  = TotGen_CPA  = TotGen_Ar39 = 0;
  TotGen_Ar42=0;
  TotGen_Neut = TotGen_Kryp = TotGen_Plon = TotGen_Rdon = 0;

  // reconstructed hits
  NTotHits = NColHits = NIndHits = 0; 

  HitView.clear();
  HitSize.clear();
  HitTPC.clear();
  HitChan.clear();
  HitTime.clear();
  HitRMS.clear();
  HitSADC.clear();
  HitInt.clear();
  HitPeak.clear();
  GenType.clear();
  MarleyIndex.clear();

  // IDEs
  NTotIDEs=0;
  IDEChannel.clear();
  IDEStartTime.clear();
  IDEEndTime.clear();
  IDEEnergy.clear();
  IDEElectrons.clear();
  IDEParticle.clear();

} // ResetVariables

//......................................................
void DAQSimAna::beginJob()
{
  // --- Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  // --- Our TTree
  fDAQSimTree = tfs->make<TTree>("DAQSimTree","DAQ simulation analysis tree");
  // General event information...
  fDAQSimTree -> Branch( "Run"   , &Run   , "Run/I"    );
  fDAQSimTree -> Branch( "SubRun", &SubRun, "SubRun/I" );
  fDAQSimTree -> Branch( "Event" , &Event , "Event/I"  );
  // Reconstructed hits...
  fDAQSimTree -> Branch( "NTotHits"  , &NTotHits  , "NTotHits/I" );
  fDAQSimTree -> Branch( "NColHits"  , &NColHits  , "NColHits/I" );
  fDAQSimTree -> Branch( "NIndHits"  , &NIndHits  , "NIndHits/I" );

  fDAQSimTree->Branch("HitView", &HitView);
  fDAQSimTree->Branch("HitSize", &HitSize);
  fDAQSimTree->Branch("HitTPC", &HitTPC);
  fDAQSimTree->Branch("HitChan", &HitChan);
  fDAQSimTree->Branch("HitTime", &HitTime);
  fDAQSimTree->Branch("HitRMS", &HitRMS);
  fDAQSimTree->Branch("HitSADC", &HitSADC);
  fDAQSimTree->Branch("HitInt", &HitInt);
  fDAQSimTree->Branch("HitPeak", &HitPeak);
  fDAQSimTree->Branch("GenType", &GenType);
  fDAQSimTree->Branch("MarleyIndex", &MarleyIndex);

  fDAQSimTree -> Branch( "TotGen_Marl", &TotGen_Marl, "TotGen_Marl/I" );
  fDAQSimTree -> Branch( "TotGen_APA" , &TotGen_APA , "TotGen_APA/I"  );
  fDAQSimTree -> Branch( "TotGen_CPA" , &TotGen_CPA , "TotGen_CPA/I"  );
  fDAQSimTree -> Branch( "TotGen_Ar39", &TotGen_Ar39, "TotGen_Ar39/I" );
  fDAQSimTree -> Branch( "TotGen_Ar42", &TotGen_Ar42, "TotGen_Ar42/I" );
  fDAQSimTree -> Branch( "TotGen_Neut", &TotGen_Neut, "TotGen_Neut/I" );
  fDAQSimTree -> Branch( "TotGen_Kryp", &TotGen_Kryp, "TotGen_Kryp/I" );
  fDAQSimTree -> Branch( "TotGen_Plon", &TotGen_Plon, "TotGen_Plon/I" );
  fDAQSimTree -> Branch( "TotGen_Rdon", &TotGen_Rdon, "TotGen_Rdon/I" );

  // IDEs
  fDAQSimTree -> Branch( "NTotIDEs"  , &NTotIDEs  , "NTotIDEs/I" );
  fDAQSimTree->Branch("IDEChannel", &IDEChannel);
  fDAQSimTree->Branch("IDEStartTime", &IDEStartTime);
  fDAQSimTree->Branch("IDEEndTime", &IDEEndTime);
  fDAQSimTree->Branch("IDEEnergy", &IDEEnergy);
  fDAQSimTree->Branch("IDEElectrons", &IDEElectrons);
  fDAQSimTree->Branch("IDEParticle", &IDEParticle);

  // --- Our Histograms...
  hAdjHits_Marl = tfs->make<TH1I>("hAdjHits_Marl", "Number of adjacent collection plane hits for MARLEY; Number of adjacent collection plane hits; Number of events"  , 21, -0.5, 20.5 );
  hAdjHits_APA  = tfs->make<TH1I>("hAdjHits_APA" , "Number of adjacent collection plane hits for APAs; Number of adjacent collection plane hits; Number of events"    , 21, -0.5, 20.5 );
  hAdjHits_CPA  = tfs->make<TH1I>("hAdjHits_CPA" , "Number of adjacent collection plane hits for CPAs; Number of adjacent collection plane hits; Number of events"    , 21, -0.5, 20.5 );
  hAdjHits_Ar39 = tfs->make<TH1I>("hAdjHits_Ar39", "Number of adjacent collection plane hits for Argon39; Number of adjacent collection plane hits; Number of events" , 21, -0.5, 20.5 );
  hAdjHits_Ar42 = tfs->make<TH1I>("hAdjHits_Ar42", "Number of adjacent collection plane hits for Argon42; Number of adjacent collection plane hits; Number of events" , 21, -0.5, 20.5 );
  hAdjHits_Neut = tfs->make<TH1I>("hAdjHits_Neut", "Number of adjacent collection plane hits for Neutrons; Number of adjacent collection plane hits; Number of events", 21, -0.5, 20.5 );
  hAdjHits_Kryp = tfs->make<TH1I>("hAdjHits_Kryp", "Number of adjacent collection plane hits for Krypton; Number of adjacent collection plane hits; Number of events" , 21, -0.5, 20.5 );
  hAdjHits_Plon = tfs->make<TH1I>("hAdjHits_Plon", "Number of adjacent collection plane hits for Polonium; Number of adjacent collection plane hits; Number of events", 21, -0.5, 20.5 );
  hAdjHits_Rdon = tfs->make<TH1I>("hAdjHits_Rdon", "Number of adjacent collection plane hits for Radon; Number of adjacent collection plane hits; Number of events"   , 21, -0.5, 20.5 );
  hAdjHits_Oth  = tfs->make<TH1I>("hAdjHits_Oth" , "Number of adjacent collection plane hits for Others; Number of adjacent collection plane hits; Number of events"  , 21, -0.5, 20.5 );
} // BeginJob

//......................................................
void DAQSimAna::SaveIDEs(art::Event const & evt)
{
    auto allParticles = evt.getValidHandle<std::vector<simb::MCParticle> >(fGEANTLabel);
    art::FindMany<simb::MCTruth> assn(allParticles,evt,fGEANTLabel);
    std::map<int, const simb::MCTruth*> idToTruth;
    for(size_t i=0; i<allParticles->size(); ++i){
        const simb::MCParticle& particle=allParticles->at(i);
        const std::vector<const simb::MCTruth*> truths=assn.at(i);
        if(truths.size()==1){
            idToTruth[particle.TrackId()]=truths[0];
        }
        else{
            mf::LogDebug("DAQSimAna") << "Particle " << particle.TrackId() << " has " << truths.size() << " truths";
            idToTruth[particle.TrackId()]=nullptr;
        }
    }

    // Get the SimChannels so we can see where the actual energy depositions were
    auto& simchs=*evt.getValidHandle<std::vector<sim::SimChannel>>("largeant");

    for(auto&& simch: simchs){
        // We only care about collection channels
        if(geo->SignalType(simch.Channel())!=geo::kCollection) continue;

        // The IDEs record energy depositions at every tick, but
        // mostly we have several depositions at contiguous times. So
        // we're going to save just one output IDE for each contiguous
        // block of hits on a channel. Each item in vector is a list
        // of (TDC, IDE*) for contiguous-in-time IDEs
        std::vector<std::vector<std::pair<int, const sim::IDE*> > > contigIDEs;
        int prevTDC=0;
        for (const auto& TDCinfo: simch.TDCIDEMap()) {
            // Do we need to start a new group of IDEs? Yes if this is
            // the first IDE in this channel. Yes if this IDE is not
            // contiguous with the previous one
            if(contigIDEs.empty() || TDCinfo.first-prevTDC>5){
                contigIDEs.push_back(std::vector<std::pair<int, const sim::IDE*> >());
            }
            std::vector<std::pair<int, const sim::IDE*> >& currentIDEs=contigIDEs.back();
            
            // Add all the current tick's IDEs to the list
            for (const sim::IDE& ide: TDCinfo.second) {
                currentIDEs.push_back(std::make_pair(TDCinfo.first, &ide));
            }
            prevTDC=TDCinfo.first;
        }

        for(auto const& contigs : contigIDEs){
            float energy=0;
            float electrons=0;
            int startTime=99999;
            int endTime=0;
            std::map<PType, float> ptypeToEnergy;
            for(auto const& timeide : contigs){
                const int tdc=timeide.first;
                startTime=std::min(tdc, startTime);
                endTime=std::max(tdc, endTime);
                const sim::IDE& ide=*timeide.second;
                const float thisEnergy=ide.energy;
                const PType thisPType=WhichParType(std::abs(ide.trackID));
                energy+=thisEnergy;
                electrons+=ide.numElectrons;
                ptypeToEnergy[thisPType]+=thisEnergy;
            }
            float bestEnergy=0;
            PType bestPType=kUnknown;
            for(auto const& it : ptypeToEnergy){
                if(it.second>bestEnergy){
                    bestEnergy=it.second;
                    bestPType=it.first;
                }
            }
            // Ignore anything past the end of the readout window
            if(startTime<4492){
                IDEChannel.push_back(simch.Channel());
                IDEStartTime.push_back(startTime);
                IDEEndTime.push_back(endTime);
                IDEEnergy.push_back(energy);
                IDEElectrons.push_back(electrons);
                IDEParticle.push_back(bestPType);
            }
        } // loop over our compressed IDEs
    } // loop over SimChannels
}

//......................................................
void DAQSimAna::analyze(art::Event const & evt)
{

  // --- We want to reset all of our TTree variables...
  ResetVariables();
 
  // --- Set all of my general event information...
  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  // --- Lift out the TPC raw digits:
  //auto rawdigits = evt.getValidHandle<std::vector<raw::RawDigit> >(fRawDigitLabel);

  // --- Lift out the reco hits:
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);

  // --- Lift out the MARLEY particles.
  auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fMARLLabel);
  std::cout << "MarlTrue.size()=" << MarlTrue->size() << std::endl;
  art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue,evt,fGEANTLabel);
  FillMyMaps( MarlParts, MarlAssn, MarlTrue, &trkIDToMarleyIndex );
  TotGen_Marl = MarlParts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of MarleyParts is " << MarlParts.size();

  // --- Lift out the APA particles.
  auto APATrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fAPALabel);
  art::FindManyP<simb::MCParticle> APAAssn(APATrue,evt,fGEANTLabel);
  FillMyMaps( APAParts, APAAssn, APATrue );
  TotGen_APA = APAParts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of APAParts is " << APAParts.size();

  // --- Lift out the CPA particles.
  auto CPATrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fCPALabel);
  art::FindManyP<simb::MCParticle> CPAAssn(CPATrue,evt,fGEANTLabel);
  FillMyMaps( CPAParts, CPAAssn, CPATrue );
  TotGen_CPA = CPAParts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of CPAParts is " << CPAParts.size();

  // --- Lift out the Ar39 particles.
  auto Ar39True = evt.getValidHandle<std::vector<simb::MCTruth> >(fAr39Label);
  art::FindManyP<simb::MCParticle> Ar39Assn(Ar39True,evt,fGEANTLabel);
  FillMyMaps( Ar39Parts, Ar39Assn, Ar39True );
  TotGen_Ar39 = Ar39Parts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of Ar39Parts is " << Ar39Parts.size();

  // --- Lift out the Ar42 particles.
  auto Ar42True = evt.getValidHandle<std::vector<simb::MCTruth> >(fAr42Label);
  art::FindManyP<simb::MCParticle> Ar42Assn(Ar42True,evt,fGEANTLabel);
  FillMyMaps( Ar42Parts, Ar42Assn, Ar42True );
  TotGen_Ar42 = Ar42Parts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of Ar42Parts is " << Ar42Parts.size();

  // --- Lift out the Neut particles.
  auto NeutTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fNeutLabel);
  art::FindManyP<simb::MCParticle> NeutAssn(NeutTrue,evt,fGEANTLabel);
  FillMyMaps( NeutParts, NeutAssn, NeutTrue );
  TotGen_Neut = NeutParts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of NeutParts is " << NeutParts.size();

  // --- Lift out the Kryp particles.
  auto KrypTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fKrypLabel);
  art::FindManyP<simb::MCParticle> KrypAssn(KrypTrue,evt,fGEANTLabel);
  FillMyMaps( KrypParts, KrypAssn, KrypTrue );
  TotGen_Kryp = KrypParts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of KrypParts is " << KrypParts.size();

  // --- Lift out the Plon particles.
  auto PlonTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fPlonLabel);
  art::FindManyP<simb::MCParticle> PlonAssn(PlonTrue,evt,fGEANTLabel);
  FillMyMaps( PlonParts, PlonAssn, PlonTrue );
  TotGen_Plon = PlonParts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of PlonParts is " << PlonParts.size();

  // --- Lift out the Rdon particles.
  auto RdonTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fRdonLabel);
  art::FindManyP<simb::MCParticle> RdonAssn(RdonTrue,evt,fGEANTLabel);
  FillMyMaps( RdonParts, RdonAssn, RdonTrue );
  TotGen_Rdon = RdonParts.size();
  mf::LogDebug("DAQSimAna") << "--- The size of RdonParts is " << RdonParts.size();

  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
      {kMarl,    MarlParts},
      {kAPA,     APAParts},
      {kCPA,     CPAParts},
      {kAr39,    Ar39Parts},
      {kAr42,    Ar42Parts},
      {kNeutron, NeutParts},
      {kKryp,    KrypParts},
      {kPlon,    PlonParts},
      {kRdon,    RdonParts}
  };

  for(auto const& it : PTypeToMap){
      const PType p=it.first;
      if(p>kNPTypes){
          std::cout << "PType is " << (int)p << std::endl;
      }
      auto const& m=it.second;
      for(auto const& it2 : m){
          trkIDToPType.insert(std::make_pair(it2.first, p));
      }
  }
  
  // --- Finally, get a list of all of my particles in one chunk.
  const sim::ParticleList& PartList = pi_serv->ParticleList();
  mf::LogDebug("DAQSimAna") << "There are a total of " << PartList.size() << " MCParticles in the event ";

  // Now that we've filled all the truth maps, we can fill a list of the true energy deposititions (IDEs)
  SaveIDEs(evt);

  std::vector< recob::Hit > ColHits_Marl;
  std::vector< recob::Hit > ColHits_CPA;
  std::vector< recob::Hit > ColHits_APA;
  std::vector< recob::Hit > ColHits_Ar39;
  std::vector< recob::Hit > ColHits_Ar42;
  std::vector< recob::Hit > ColHits_Neut;
  std::vector< recob::Hit > ColHits_Kryp;
  std::vector< recob::Hit > ColHits_Plon;
  std::vector< recob::Hit > ColHits_Rdon;
  std::vector< recob::Hit > ColHits_Oth;

  //*
  // --- Loop over the reconstructed hits to determine the "size" of each hit 
  NTotHits = reco_hits->size();

  for(int hit = 0; hit < NTotHits; ++hit) {
    // --- Let access this particular hit.
    recob::Hit const& ThisHit = reco_hits->at(hit);  
    
    // --- Lets figure out which particle contributed the most charge to this hit...
    int MainTrID    = -1;
    double TopEFrac = -DBL_MAX;

    // HitToTrackIDEs opens a specific window around the hit. I want a
    // wider one, because the filtering can delay the hit. So this bit
    // is a copy of HitToTrackIDEs from the backtracker, with some
    // modification
    const double start = ThisHit.PeakTime()-20;
    const double end   = ThisHit.PeakTime()+ThisHit.RMS()+20;
    std::vector<sim::TrackIDE> ThisHitIDE = bt_serv->ChannelToTrackIDEs(ThisHit.Channel(), start, end);
    // The old method
    // std::vector< sim::TrackIDE > ThisHitIDE = bt_serv->HitToTrackIDEs( ThisHit );
    for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL) {
        if ( ThisHitIDE[ideL].energyFrac > TopEFrac ) {
            TopEFrac = ThisHitIDE[ideL].energyFrac;
            MainTrID = ThisHitIDE[ideL].trackID;
        }
    }
    // --- Lets figure out how that particle was generated...
    PType ThisPType = WhichParType( MainTrID );
    int thisMarleyIndex=-1;
    if(ThisPType==kMarl){
        auto const it=trkIDToMarleyIndex.find(MainTrID);
        if(it==trkIDToMarleyIndex.end()){
            std::cout << "Track ID " << MainTrID << " is not in Marley index map" << std::endl;
        }
        else{
            thisMarleyIndex=it->second;
        }
    }
    // --- Write out some information about this hit....
    // std::cout << "Looking at hit on channel " << ThisHit.Channel() << " corresponding to TPC " << ThisHit.WireID().TPC << ", wire " << ThisHit.WireID().Wire << ", plane " << ThisHit.WireID().Plane << ".\n"
	//       << "\tIt was at time " << ThisHit.PeakTime() << ", with amplitude " << ThisHit.PeakAmplitude() << ", it was caused by " << ThisHitIDE.size() << " particles, the main one being"
	//       << " TrackID " << MainTrID << " which was generated by " << ThisPType
	//       << std::endl;
    //*/
    // --- Check which view this hit is on...
    if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) {
      ++NIndHits;
    } else { // If not induction then must be collection.
      ++NColHits;
    }

    // --- Now fill in all of the hit level variables.
    HitView.push_back(ThisHit.View());
    HitSize.push_back(ThisHit.EndTick() - ThisHit.StartTick());
    HitTPC .push_back(ThisHit.WireID().TPC);
    HitChan.push_back(ThisHit.Channel());
    HitTime.push_back(ThisHit.PeakTime());
    HitRMS .push_back(ThisHit.RMS());
    HitSADC.push_back(ThisHit.SummedADC());
    HitInt .push_back(ThisHit.Integral());
    HitPeak.push_back(ThisHit.PeakAmplitude());
    GenType.push_back(ThisPType);
    MarleyIndex.push_back(thisMarleyIndex);
  
    // --- I want to fill a vector of coll plane hits, for each of the different kinds of generator.
    if (ThisHit.View() == 2) {
      if (ThisPType == kUnknown)      ColHits_Oth .push_back( ThisHit );
      else if (ThisPType == kMarl)    ColHits_Marl.push_back( ThisHit );
      else if (ThisPType == kAPA)     ColHits_APA .push_back( ThisHit );
      else if (ThisPType == kCPA)     ColHits_CPA .push_back( ThisHit );
      else if (ThisPType == kAr39)    ColHits_Ar39.push_back( ThisHit );
      else if (ThisPType == kAr42)    ColHits_Ar42.push_back( ThisHit );
      else if (ThisPType == kNeutron) ColHits_Neut.push_back( ThisHit );
      else if (ThisPType == kKryp)    ColHits_Kryp.push_back( ThisHit );
      else if (ThisPType == kPlon)    ColHits_Plon.push_back( ThisHit );
      else if (ThisPType == kRdon)    ColHits_Rdon.push_back( ThisHit );
    }
  } // Loop over reco_hits.

  // ---- Write out the Marley hits....
  mf::LogDebug("DAQSimAna") << "\n\nAfter all of that I have a total of " << ColHits_Marl.size() << " MARLEY col plane hits.";
  for (size_t hh=0; hh<ColHits_Marl.size(); ++hh) {
    mf::LogDebug("DAQSimAna") << "\tHit " << hh << " was on chan " << ColHits_Marl[hh].Channel() << " at " << ColHits_Marl[hh].PeakTime();
  }
  if(fDoCalcAdjHits){
      // --- Now calculate all of the hits...
      CalcAdjHits( ColHits_Marl, hAdjHits_Marl, false );
      mf::LogDebug("DAQSimAna") << "\nAnd now for APA hits...";
      CalcAdjHits( ColHits_APA , hAdjHits_APA , false  );
      mf::LogDebug("DAQSimAna") << "\nAnd now for CPA hits...";
      CalcAdjHits( ColHits_CPA , hAdjHits_CPA , false  );
      mf::LogDebug("DAQSimAna") << "\nAnd now for Ar39 hits...";
      CalcAdjHits( ColHits_Ar39, hAdjHits_Ar39, false );
      mf::LogDebug("DAQSimAna") << "\nAnd now for Ar42 hits...";
      CalcAdjHits( ColHits_Ar42, hAdjHits_Ar42, false );
      mf::LogDebug("DAQSimAna") << "\nAnd now for Neuton hits...";
      CalcAdjHits( ColHits_Neut, hAdjHits_Neut, false );
      mf::LogDebug("DAQSimAna") << "\nAnd now for Krypton hits...";
      CalcAdjHits( ColHits_Kryp, hAdjHits_Kryp, false );
      mf::LogDebug("DAQSimAna") << "\nAnd now for Polonium hits...";
      CalcAdjHits( ColHits_Plon, hAdjHits_Plon, false );
      mf::LogDebug("DAQSimAna") << "\nAnd now for Radon hits...";
      CalcAdjHits( ColHits_Rdon, hAdjHits_Rdon, false );
      mf::LogDebug("DAQSimAna") << "\nAnd now for Other hits...";
      CalcAdjHits( ColHits_Oth , hAdjHits_Oth , false  );
  }

  fDAQSimTree -> Fill();

} // Analyze DAQSimAna.


//......................................................
void DAQSimAna::FillMyMaps( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand,
                            std::map<int, int>* indexMap)
{
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar;
      if(indexMap) indexMap->insert({ThisPar.TrackId(), L1});
    }
  }
  return;
}

//......................................................
PType DAQSimAna::WhichParType( int TrID )
{
    PType ThisPType = kUnknown;
    auto const& it=trkIDToPType.find(TrID);
    if(it!=trkIDToPType.end()){
        ThisPType=it->second;
    }
    if(ThisPType>kNPTypes){
        std::cout << "In WhichParType, ptype is " << (int)ThisPType << std::endl;
    }   
    return ThisPType;
}

//......................................................
PType DAQSimAna::WhichParType( const art::ValidHandle<simb::MCTruth>& truthHand)
{
  const std::string& label=truthHand.provenance()->moduleLabel();
   if(label==fMARLLabel){ return kMarl; }
   if(label==fAPALabel) { return kAPA;  }
   if(label==fCPALabel) { return kCPA;  }
   if(label==fAr39Label){ return kAr39; }
   if(label==fAr42Label){ return kAr42; }
   if(label==fNeutLabel){ return kNeutron; }
   if(label==fKrypLabel){ return kKryp; }
   if(label==fPlonLabel){ return kPlon; }
   if(label==fRdonLabel){ return kRdon; }
   return kUnknown;
}

//......................................................
bool DAQSimAna::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap )
{
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if ( ParIt != ParMap.end() ) {
    return true;
  } else
    return false;
}

//......................................................
void DAQSimAna::CalcAdjHits( std::vector< recob::Hit > MyVec, TH1I* MyHist, bool HeavDebug ) {
  const double TimeRange = 20;
  const int    ChanRange = 2;
  unsigned int FilledHits = 0;
  unsigned int NumOriHits = MyVec.size();
  while( NumOriHits != FilledHits ) {
    if (HeavDebug) mf::LogDebug("DAQSimAna") << "\nStart of my while loop";
    std::vector< recob::Hit > AdjHitVec;
    AdjHitVec.push_back ( MyVec[0] );
    MyVec.erase( MyVec.begin()+0 );
    int LastSize = 0;
    int NewSize  = AdjHitVec.size();
    while ( LastSize != NewSize ) {
      std::vector<int> AddNow;
      for (size_t aL=0; aL < AdjHitVec.size(); ++aL) {
        for (size_t nL=0; nL < MyVec.size(); ++nL) {
          if (HeavDebug) {
            mf::LogDebug("DAQSimAna") << "\t\tLooping though AdjVec " << aL << " and  MyVec " << nL
                                      << " AdjHitVec - " << AdjHitVec[aL].Channel() << " & " << AdjHitVec[aL].PeakTime()
                                      << " MVec - " << MyVec[nL].Channel() << " & " << MyVec[nL].PeakTime()
                                      << " Channel " << abs( (int)AdjHitVec[aL].Channel()  - (int)MyVec[nL].Channel()  )  << " bool " << (bool)(abs( (int)AdjHitVec[aL].Channel() - (int)MyVec[nL].Channel()  ) <= ChanRange)
                                      << " Time " << abs( AdjHitVec[aL].PeakTime() - MyVec[nL].PeakTime() ) << " bool " << (bool)(abs( (double)AdjHitVec[aL].PeakTime() - (double)MyVec[nL].PeakTime() ) <= TimeRange);
		  
          }
          if ( abs( (int)AdjHitVec[aL].Channel()  - (int)MyVec[nL].Channel()  ) <= ChanRange &&
               abs( (double)AdjHitVec[aL].PeakTime() - (double)MyVec[nL].PeakTime() ) <= TimeRange ) {
            if (HeavDebug) mf::LogDebug("DAQSimAna") << "\t\t\tFound a new thing!!!";
            // --- Check that this element isn't already in AddNow.
            bool AlreadyPres = false;
            for (size_t zz=0; zz<AddNow.size(); ++zz) {
              if (AddNow[zz] == (int)nL) AlreadyPres = true;
            }
            if (!AlreadyPres)
              AddNow.push_back( nL );
          } // If this hit is within the window around one of my other hits.
        } // Loop through my vector of colleciton plane hits.
      } // Loop through AdjHitVec
      // --- Now loop through AddNow and remove from Marley whilst adding to AdjHitVec
      for (size_t aa=0; aa<AddNow.size(); ++aa) {
        if (HeavDebug) {
          mf::LogDebug("DAQSimAna") << "\tRemoving element " << AddNow.size()-1-aa << " from MyVec ===> "
                                    << MyVec[ AddNow[AddNow.size()-1-aa] ].Channel() << " & " << MyVec[ AddNow[AddNow.size()-1-aa] ].PeakTime();

        }
        AdjHitVec.push_back ( MyVec[ AddNow[AddNow.size()-1-aa] ] );
        MyVec.erase( MyVec.begin() + AddNow[AddNow.size()-1-aa] );
      }
      LastSize = NewSize;
      NewSize  = AdjHitVec.size();
      if (HeavDebug) {
        mf::LogDebug("DAQSimAna") << "\t---After that pass, AddNow was size " << AddNow.size() << " ==> LastSize is " << LastSize << ", and NewSize is " << NewSize
                  << "\nLets see what is in AdjHitVec....";
        for (size_t aL=0; aL < AdjHitVec.size(); ++aL) {
          mf::LogDebug("DAQSimAna") << "\tElement " << aL << " is ===> " << AdjHitVec[aL].Channel() << " & " << AdjHitVec[aL].PeakTime();
        }
      }

    } // while ( LastSize != NewSize )
    int NumAdjColHits = AdjHitVec.size();
    if (HeavDebug) mf::LogDebug("DAQSimAna") << "After that loop, I had " << NumAdjColHits << " adjacent collection plane hits.";
    MyHist -> Fill( NumAdjColHits );
    FilledHits += NumAdjColHits;
  }
  return;
}

//......................................................
DEFINE_ART_MODULE(DAQSimAna)
