#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SupernovaTruth.h"

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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

const int nMaxHits    = 500000;
//const int nMaxSimIDEs = 1000000; // unused
//const int nMaxDigs = 4492; // unused

enum PType{ kUnknown, kMarl, kAPA, kCPA, kAr39, kNeut, kKryp, kPlon, kRdon , kAr42};

class SNAna : public art::EDAnalyzer {

public:

  explicit SNAna(fhicl::ParameterSet const & p);

  SNAna(SNAna const &) = delete;
  SNAna(SNAna &&) = delete;
  SNAna & operator = (SNAna const &) = delete;
  SNAna & operator = (SNAna &&) = delete;

  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void endJob() override;

private:

  void ResetVariables();
  void FillMyMaps  ( std::map< int, simb::MCParticle> &MyMap, 
                     art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand );
  PType WhichParType( int TrID );
  bool  InMyMap     ( int TrID, std::map< int, simb::MCParticle> ParMap );

  int firstCatch;
  int secondCatch;
  int thirdCatch;

  std::string fRawDigitLabel;
  std::string fHitLabel;

  std::string fGEANTLabel;
  std::string fMARLLabel; std::map< int, simb::MCParticle > MarlParts;
  std::string fAPALabel;  std::map< int, simb::MCParticle > APAParts;
  std::string fCPALabel;  std::map< int, simb::MCParticle > CPAParts;
  std::string fAr39Label; std::map< int, simb::MCParticle > Ar39Parts;
  std::string fNeutLabel; std::map< int, simb::MCParticle > NeutParts;
  std::string fKrypLabel; std::map< int, simb::MCParticle > KrypParts;
  std::string fPlonLabel; std::map< int, simb::MCParticle > PlonParts;
  std::string fRdonLabel; std::map< int, simb::MCParticle > RdonParts;
  std::string fAr42Label; std::map< int, simb::MCParticle > Ar42Parts;
  std::map<int, const simb::MCParticle*> truthmap;

  TTree* fSNAnaTree;

  int Run;
  int SubRun;
  int Event;

  int   NTotHits;
  int   NColHits;
  int   NIndHits;
  int   nHitsNotBackTracked;
  int   HitView[nMaxHits];
  int   HitSize[nMaxHits];
  int   HitTPC [nMaxHits];
  int   HitChan[nMaxHits];
  float HitTime[nMaxHits];
  float HitRMS [nMaxHits];
  float HitSADC[nMaxHits];
  float HitInt [nMaxHits];
  float HitPeak[nMaxHits];
  int   GenType[nMaxHits];
  float Hit_X[nMaxHits];                  
  float Hit_Y[nMaxHits];                  
  float Hit_Z[nMaxHits];                  
  float Hit_Energy[nMaxHits];             
  float Hit_NumElectrons[nMaxHits];
  int   NCorrespondingIDEs[nMaxHits]; 
  int   TotGen_Marl;
  int   TotGen_APA;
  int   TotGen_CPA;
  int   TotGen_Ar39;
  int   TotGen_Neut;
  int   TotGen_Kryp;
  int   TotGen_Plon;
  int   TotGen_Rdon;
  int   TotGen_Ar42;

  int    VertexChan;
  int    Nu_Type;
  int    Nu_Lep_Type;
  int    Mode;
  int    CCNC;
  int    HitNucleon;
  int    Target;
  std::vector<int>    MarlSample;
  std::vector<double> MarlTime;
  std::vector<double> MarlWeight;
  double ENu;
  double ENu_Lep;
  double VertX;
  double VertY;
  double VertZ;
  double VertexT;
  double Px;
  double Py;
  double Pz;

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
};

SNAna::SNAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}


void SNAna::reconfigure(fhicl::ParameterSet const & p)
{
  fRawDigitLabel = p.get<std::string> ("RawDigitLabel");  
  fHitLabel      = p.get<std::string> ("HitLabel");

  fGEANTLabel = p.get<std::string> ("GEANT4Label");
  fMARLLabel = p.get<std::string> ("MARLEYLabel");
  fAPALabel  = p.get<std::string> ("APALabel");
  fCPALabel  = p.get<std::string> ("CPALabel");
  fAr39Label = p.get<std::string> ("Argon39Label");
  fNeutLabel = p.get<std::string> ("NeutronLabel");
  fKrypLabel = p.get<std::string> ("KryptonLabel");
  fPlonLabel = p.get<std::string> ("PoloniumLabel");
  fRdonLabel = p.get<std::string> ("RadonLabel");
  fAr42Label = p.get<std::string> ("Argon42Label");
} 


void SNAna::ResetVariables()
{
  MarlParts.clear(); APAParts .clear(); CPAParts .clear(); Ar39Parts.clear();
  NeutParts.clear(); KrypParts.clear(); PlonParts.clear(); RdonParts.clear();
  Ar42Parts.clear();

  Run = SubRun = Event = -1;
  
  TotGen_Marl = TotGen_APA  = TotGen_CPA  = TotGen_Ar39 = 0;
  TotGen_Neut = TotGen_Kryp = TotGen_Plon = TotGen_Rdon = 0;
  TotGen_Ar42 = 0;

  NTotHits = NColHits = NIndHits = 0; 
  nHitsNotBackTracked = 0;
  for (int hh=0; hh<nMaxHits; ++hh) {
    HitView[hh] = HitSize[hh] = HitChan[hh] = GenType[hh] = 0;
    HitTime[hh] = HitRMS [hh] = HitSADC[hh] = 0;
    HitInt [hh] = HitPeak[hh] = HitTPC [hh] = 0;
    NCorrespondingIDEs[hh] = 0;
    Hit_X[hh] = 0;                  
    Hit_Y[hh] = 0;                  
    Hit_Z[hh] = 0;                  
    Hit_Energy[hh] = 0;             
    Hit_NumElectrons[hh] = 0;
  } 

  MarlSample.clear();
  MarlTime.clear();
  MarlWeight.clear();
  ENu = 0;
  ENu_Lep = 0;
  VertX = 0;
  VertY = 0;
  VertZ = 0;
  VertexT = 0;
  Px    = 0;
  Py    = 0;
  Pz    = 0;
  VertexChan = 0;
} 


void SNAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fSNAnaTree = tfs->make<TTree>("SNSimTree","SN simulation analysis tree");

  fSNAnaTree -> Branch( "Run"   , &Run   , "Run/I"    );
  fSNAnaTree -> Branch( "SubRun", &SubRun, "SubRun/I" );
  fSNAnaTree -> Branch( "Event" , &Event , "Event/I"  );

  fSNAnaTree -> Branch( "NTotHits"  , &NTotHits  , "NTotHits/I" );
  fSNAnaTree -> Branch( "NColHits"  , &NColHits  , "NColHits/I" );
  fSNAnaTree -> Branch( "NIndHits"  , &NIndHits  , "NIndHits/I" );
  fSNAnaTree -> Branch( "nHitsNotBackTracked", &nHitsNotBackTracked, "nHitsNotBackTracked/I" );
  fSNAnaTree -> Branch( "HitView"   , &HitView   , "HitView[NColHits]/I" );
  fSNAnaTree -> Branch( "HitSize"   , &HitSize   , "HitSize[NColHits]/I" );
  fSNAnaTree -> Branch( "HitTPC"    , &HitTPC    , "HitTPC[NColHits]/I"  );
  fSNAnaTree -> Branch( "HitChan"   , &HitChan   , "HitChan[NColHits]/I" );
  fSNAnaTree -> Branch( "HitTime"   , &HitTime   , "HitTime[NColHits]/F" );
  fSNAnaTree -> Branch( "HitRMS"    , &HitRMS    , "HitRMS[NColHits]/F"  );
  fSNAnaTree -> Branch( "HitSADC"   , &HitSADC   , "HitSADC[NColHits]/F" );
  fSNAnaTree -> Branch( "HitInt"    , &HitInt    , "HitInt[NColHits]/F"  );
  fSNAnaTree -> Branch( "HitPeak"   , &HitPeak   , "HitPeak[NColHits]/F" );
  fSNAnaTree -> Branch( "GenType"   , &GenType   , "GenType[NColHits]/I" ); 
                                                                            
  fSNAnaTree -> Branch( "NCorrespondingIDEs", &NCorrespondingIDEs, "NCorrespondingIDEs[NColHits]/I" );
  fSNAnaTree -> Branch( "Hit_X"     , &Hit_X   , "Hit_X[NColHits]/F" );
  fSNAnaTree -> Branch( "Hit_Y"     , &Hit_Y   , "Hit_Y[NColHits]/F" );
  fSNAnaTree -> Branch( "Hit_Z"     , &Hit_Z     , "Hit_Z[NColHits]/F" );
  fSNAnaTree -> Branch( "Hit_Energy", &Hit_Energy, "Hit_Energy[NColHits]/F" );
  fSNAnaTree -> Branch( "Hit_NumElectrons", &Hit_NumElectrons, "Hit_NumElectrons[NColHits]/F" );

  fSNAnaTree -> Branch( "TotGen_Marl", &TotGen_Marl, "TotGen_Marl/I" );
  fSNAnaTree -> Branch( "TotGen_APA" , &TotGen_APA , "TotGen_APA/I"  );
  fSNAnaTree -> Branch( "TotGen_CPA" , &TotGen_CPA , "TotGen_CPA/I"  );
  fSNAnaTree -> Branch( "TotGen_Ar39", &TotGen_Ar39, "TotGen_Ar39/I" );
  fSNAnaTree -> Branch( "TotGen_Neut", &TotGen_Neut, "TotGen_Neut/I" );
  fSNAnaTree -> Branch( "TotGen_Kryp", &TotGen_Kryp, "TotGen_Kryp/I" );
  fSNAnaTree -> Branch( "TotGen_Plon", &TotGen_Plon, "TotGen_Plon/I" );
  fSNAnaTree -> Branch( "TotGen_Rdon", &TotGen_Rdon, "TotGen_Rdon/I" );
  fSNAnaTree -> Branch( "TotGen_Ar42", &TotGen_Ar42, "TotGen_Ar42/I" );

  fSNAnaTree->Branch("Nu_Type", &Nu_Type, "Nu_Type/I");
  fSNAnaTree->Branch("Nu_Lep_Type", &Nu_Lep_Type, "Nu_Lep_Type/I");
  fSNAnaTree->Branch("Mode", &Mode, "Mode/I");
  fSNAnaTree->Branch("CCNC", &CCNC, "CCNC/I");
  fSNAnaTree->Branch("HitNucleon", &HitNucleon, "HitNucleon/I");
  fSNAnaTree->Branch("Target", &Target, "Target/I");
  fSNAnaTree->Branch("MarlSample", &MarlSample);
  fSNAnaTree->Branch("MarlTime", &MarlTime);
  fSNAnaTree->Branch("MarlWeight",  &MarlWeight);
  fSNAnaTree->Branch("ENu", &ENu, "ENu/D");
  fSNAnaTree->Branch("Px",  &Px,  "Px/D");
  fSNAnaTree->Branch("Py",  &Py,  "Py/D");
  fSNAnaTree->Branch("Pz",  &Pz,  "Pz/D");
  fSNAnaTree->Branch("ENu_Lep", &ENu_Lep, "ENu_Lep/D");
  fSNAnaTree->Branch("VertX", &VertX, "VertX/D");
  fSNAnaTree->Branch("VertY", &VertY, "VertY/D");
  fSNAnaTree->Branch("VertZ", &VertZ, "VertZ/D");
  fSNAnaTree->Branch("VertexT", &VertexT, "VertexT/D");
  fSNAnaTree->Branch("VertexChan", &VertexChan, "VertexChan/I");
}


void SNAna::analyze(art::Event const & evt)
{
  ResetVariables();
 
  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  const detinfo::DetectorProperties* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  //GET THE RECO HITS.
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);

  //LIFT OUT THE MARLEY PARTICLES.
  auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fMARLLabel);
  art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue,evt,fGEANTLabel);
  FillMyMaps( MarlParts, MarlAssn, MarlTrue );
  TotGen_Marl = MarlParts.size();

  //SUPERNOVA TRUTH.
  art::FindManyP<sim::SupernovaTruth> SNTruth(MarlTrue, evt, fMARLLabel);

  double Px_(0), Py_(0), Pz_(0), Pnorm(1);
  for(unsigned int i = 0; i < MarlTrue->size(); i++)
  {
    Nu_Type     = MarlTrue->at(i).GetNeutrino().Nu().PdgCode();
    Nu_Lep_Type = MarlTrue->at(i).GetNeutrino().Lepton().PdgCode(); 
    Mode        = MarlTrue->at(i).GetNeutrino().Mode();
    CCNC        = MarlTrue->at(i).GetNeutrino().CCNC();
    Target      = MarlTrue->at(i).GetNeutrino().Target();
    HitNucleon  = MarlTrue->at(i).GetNeutrino().HitNuc();
    VertX   = MarlTrue->at(i).GetNeutrino().Lepton().Vx(); 
    VertY   = MarlTrue->at(i).GetNeutrino().Lepton().Vy(); 
    VertZ   = MarlTrue->at(i).GetNeutrino().Lepton().Vz(); 
    ENu_Lep = MarlTrue->at(i).GetNeutrino().Lepton().E();
    ENu     = MarlTrue->at(i).GetNeutrino().Nu().E();
    Px_     = MarlTrue->at(i).GetNeutrino().Lepton().Px();
    Py_     = MarlTrue->at(i).GetNeutrino().Lepton().Py();
    Pz_     = MarlTrue->at(i).GetNeutrino().Lepton().Pz();
    for (unsigned int j = 0; j < SNTruth.at(i).size(); j++) 
    {
      const sim::SupernovaTruth ThisTr = (*SNTruth.at(i).at(j));
      MarlTime.push_back(ThisTr.SupernovaTime);
      MarlWeight.push_back(ThisTr.Weight);
      MarlSample.push_back(ThisTr.SamplingMode);
    }
  }
  Pnorm = std::sqrt(Px_*Px_+Py_*Py_+Pz_*Pz_);
  Px = Px_/Pnorm;
  Py = Py_/Pnorm;
  Pz = Pz_/Pnorm;

  bool caught = false;
  double Vertex[3] = {VertX, VertY, VertZ};
  geo::WireID WireID;
  geo::PlaneID Plane(geo->FindTPCAtPosition(Vertex),geo::kZ);
  try
  {
    WireID = geo->NearestWireID(Vertex, Plane);
  }
  catch(...)
  {
    caught = true;
  }
  if(caught==true)
  {
    VertexChan = -1; 
  }
  else
  {
    VertexChan = geo->PlaneWireToChannel(WireID);
  }

  //CM/MICROSECOND.
  double drift_velocity = detp->DriftVelocity(detp->Efield(),detp->Temperature());
  //CM/TICK
  drift_velocity = drift_velocity*0.5;
  VertexT = VertX/drift_velocity;

  auto APATrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fAPALabel);
  art::FindManyP<simb::MCParticle> APAAssn(APATrue,evt,fGEANTLabel);
  FillMyMaps( APAParts, APAAssn, APATrue );
  TotGen_APA = APAParts.size();

  auto CPATrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fCPALabel);
  art::FindManyP<simb::MCParticle> CPAAssn(CPATrue,evt,fGEANTLabel);
  FillMyMaps( CPAParts, CPAAssn, CPATrue );
  TotGen_CPA = CPAParts.size();

  auto Ar39True = evt.getValidHandle<std::vector<simb::MCTruth> >(fAr39Label);
  art::FindManyP<simb::MCParticle> Ar39Assn(Ar39True,evt,fGEANTLabel);
  FillMyMaps( Ar39Parts, Ar39Assn, Ar39True );
  TotGen_Ar39 = Ar39Parts.size();

  auto NeutTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fNeutLabel);
  art::FindManyP<simb::MCParticle> NeutAssn(NeutTrue,evt,fGEANTLabel);
  FillMyMaps( NeutParts, NeutAssn, NeutTrue );
  TotGen_Neut = NeutParts.size();

  auto KrypTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fKrypLabel);
  art::FindManyP<simb::MCParticle> KrypAssn(KrypTrue,evt,fGEANTLabel);
  FillMyMaps( KrypParts, KrypAssn, KrypTrue );
  TotGen_Kryp = KrypParts.size();

  auto PlonTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fPlonLabel);
  art::FindManyP<simb::MCParticle> PlonAssn(PlonTrue,evt,fGEANTLabel);
  FillMyMaps( PlonParts, PlonAssn, PlonTrue );
  TotGen_Plon = PlonParts.size();

  auto RdonTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fRdonLabel);
  art::FindManyP<simb::MCParticle> RdonAssn(RdonTrue,evt,fGEANTLabel);
  FillMyMaps( RdonParts, RdonAssn, RdonTrue );
  TotGen_Rdon = RdonParts.size();

  auto Ar42True = evt.getValidHandle<std::vector<simb::MCTruth> >(fAr42Label);
  art::FindManyP<simb::MCParticle> Ar42Assn(Ar42True,evt,fGEANTLabel);
  FillMyMaps( Ar42Parts, Ar42Assn, Ar42True );
  TotGen_Ar42 = Ar42Parts.size();

  std::cout << "THE EVENTS NUMBER IS: " << Event << std::endl;

  std::vector< recob::Hit > ColHits_Marl;
  std::vector< recob::Hit > ColHits_CPA;
  std::vector< recob::Hit > ColHits_APA;
  std::vector< recob::Hit > ColHits_Ar39;
  std::vector< recob::Hit > ColHits_Neut;
  std::vector< recob::Hit > ColHits_Kryp;
  std::vector< recob::Hit > ColHits_Plon;
  std::vector< recob::Hit > ColHits_Rdon;
  std::vector< recob::Hit > ColHits_Oth;
  std::vector< recob::Hit > ColHits_Ar42;

  NTotHits = reco_hits->size();
  int colHitCount(0);
  int LoopHits = std::min( NTotHits, nMaxHits );
  std::cout << "---- There are " << NTotHits << " hits in the event, but array is of size " 
            << nMaxHits << ", so looping over first " << LoopHits << " hits." << std::endl;

  for(int hit = 0; hit < LoopHits; ++hit) 
  {
    recob::Hit const& ThisHit = reco_hits->at(hit);  
    if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) {
      ++NIndHits;
    } else { 
      ++NColHits;
    }
  }

  for(int hit = 0; hit < LoopHits; ++hit) 
  {
    recob::Hit const& ThisHit = reco_hits->at(hit);  

    if (ThisHit.View() == 2) {
    std::vector< sim::TrackIDE > ThisHitIDE; 
    //GETTING HOLD OF THE SIM::IDEs.
    std::vector<const sim::IDE*> ThisSimIDE;
    try
    {
      ThisHitIDE = bt_serv->HitToTrackIDEs( ThisHit );
    }
    catch(...)
    {
      std::cout << "FIRST CATCH" << std::endl;
      firstCatch++;
      try
      {
        ThisSimIDE = bt_serv->HitToSimIDEs_Ps(ThisHit);
      }
      catch(...)
      {
        std::cout << "SECOND CATCH" << std::endl; 
        secondCatch++;
        continue;
      }
      continue;
    }
    try
    {
      ThisSimIDE = bt_serv->HitToSimIDEs_Ps(ThisHit);
    }
    catch(...)
    {
      std::cout << "THIRD CATCH" << std::endl; 
      thirdCatch++;
      continue;
    }
    
    HitView[colHitCount] = ThisHit.View();
    HitSize[colHitCount] = ThisHit.EndTick() - ThisHit.StartTick();
    HitTPC [colHitCount] = ThisHit.WireID().TPC;
    HitChan[colHitCount] = ThisHit.Channel();
    HitTime[colHitCount] = ThisHit.PeakTime();
    HitRMS [colHitCount] = ThisHit.RMS();
    HitSADC[colHitCount] = ThisHit.SummedADC();
    HitInt [colHitCount] = ThisHit.Integral();
    HitPeak[colHitCount] = ThisHit.PeakAmplitude();
    NCorrespondingIDEs[colHitCount] = ThisHitIDE.size();
    if(ThisHitIDE.size()==0)
    {
      nHitsNotBackTracked++;
    }

    //WHICH PARTICLE CONTRIBUTED MOST TO THE HIT.
    int    MainTrID = -1;
    double TopEFrac = -DBL_MAX;
    for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL) 
    {
      if ( ThisHitIDE[ideL].energyFrac > TopEFrac ) 
      {
        TopEFrac = ThisHitIDE[ideL].energyFrac;
        MainTrID = ThisHitIDE[ideL].trackID;
      }
    }

    PType ThisPType      = WhichParType( MainTrID );
    GenType[colHitCount] = ThisPType;

    if(MainTrID == -1)
    {
      Hit_X[colHitCount]            = -1;
      Hit_Y[colHitCount]            = -1;
      Hit_Z[colHitCount]            = -1;
      Hit_Energy[colHitCount]       = -1;
      Hit_NumElectrons[colHitCount] = -1;
    }
    else
    {
      for(unsigned int i = 0; i < ThisSimIDE.size(); i++)
      {
        if(ThisSimIDE.at(i)->trackID==MainTrID)
        {
          Hit_X[colHitCount]            = ThisSimIDE.at(i)->x;
          Hit_Y[colHitCount]            = ThisSimIDE.at(i)->y;
          Hit_Z[colHitCount]            = ThisSimIDE.at(i)->z;
          Hit_Energy[colHitCount]       = ThisSimIDE.at(i)->energy;
          Hit_NumElectrons[colHitCount] = ThisSimIDE.at(i)->numElectrons;
          break;
        }
      }
    }
    
    if (ThisPType == 0){ColHits_Oth .push_back( ThisHit );
    }
    else if (ThisPType == 1) {ColHits_Marl.push_back( ThisHit );
    }
    else if (ThisPType == 2) {ColHits_APA .push_back( ThisHit );
    }
    else if (ThisPType == 3) {ColHits_CPA .push_back( ThisHit );
    }
    else if (ThisPType == 4) {ColHits_Ar39.push_back( ThisHit );
    }
    else if (ThisPType == 5) {ColHits_Neut.push_back( ThisHit );
    }
    else if (ThisPType == 6) {ColHits_Kryp.push_back( ThisHit );
    }
    else if (ThisPType == 7) {ColHits_Plon.push_back( ThisHit );
    }
    else if (ThisPType == 8) {ColHits_Rdon.push_back( ThisHit );
    }
    else if (ThisPType == 9) {ColHits_Ar42.push_back( ThisHit );
    }

    colHitCount++;
    }
  } 

  std::cerr << "\n\nAfter all of that I have a total of " << ColHits_Marl.size() << " MARLEY col plane hits." << std::endl;
  for (size_t hh=0; hh<ColHits_Marl.size(); ++hh) {
    std::cerr << "\tHit " << hh << " was on chan " << ColHits_Marl[hh].Channel() << " at " << ColHits_Marl[hh].PeakTime() << std::endl;
  }
  fSNAnaTree -> Fill();
} 


void SNAna::FillMyMaps( std::map< int, simb::MCParticle> &MyMap, 
                            art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand )
{
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar;
    }
  }
  return;
}


PType SNAna::WhichParType( int TrID )
{
  if ( InMyMap( TrID, Ar39Parts ) ) {
    return kAr39;
  } else  if ( InMyMap( TrID, MarlParts ) ) {
    return kMarl;
  } else if ( InMyMap( TrID, APAParts  ) ) {
    return kAPA;
  } else if ( InMyMap( TrID, CPAParts  ) ) {
    return kCPA;
  } else if ( InMyMap( TrID, NeutParts ) ) {
    return kNeut;
  } else if ( InMyMap( TrID, KrypParts ) ) {
    return kKryp;
  } else if ( InMyMap( TrID, PlonParts ) ) {
    return kPlon;
  } else if ( InMyMap( TrID, RdonParts ) ) {
    return kRdon;
  } else if (InMyMap( TrID, Ar42Parts)){
    return kAr42;
  }
  return kUnknown;
}


bool SNAna::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap )
{
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if ( ParIt != ParMap.end() ) {
    return true;
  } else
    return false;
}


void SNAna::endJob()
{
  std::cout << firstCatch << " " << secondCatch << " " << thirdCatch << std::endl; 
}

DEFINE_ART_MODULE(SNAna)
