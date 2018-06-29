
//ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"

//LArSoft includes
#include "dune/OpticalDetector/OpFlashSort.h"

#include "larcore/Geometry/Geometry.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SupernovaTruth.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"

//ART includes
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
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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
  std::string fCalDataModuleLabel;
  std::string fOpFlashModuleLabel;
  std::string fOpHitModuleLabel;

  std::string fGEANTLabel;
  std::string fMARLLabel ; std::map< int, simb::MCParticle > MarlParts;
  std::string fAPALabel  ; std::map< int, simb::MCParticle > APAParts ;
  std::string fCPALabel  ; std::map< int, simb::MCParticle > CPAParts ;
  std::string fAr39Label ; std::map< int, simb::MCParticle > Ar39Parts;
  std::string fNeutLabel ; std::map< int, simb::MCParticle > NeutParts;
  std::string fKrypLabel ; std::map< int, simb::MCParticle > KrypParts;
  std::string fPlonLabel ; std::map< int, simb::MCParticle > PlonParts;
  std::string fRdonLabel ; std::map< int, simb::MCParticle > RdonParts;
  std::string fAr42Label ; std::map< int, simb::MCParticle > Ar42Parts;
  std::map<int, const simb::MCParticle*> truthmap;

  // Mapping from track ID to particle type, for use in WhichParType()
  std::map<int, PType> trkIDToPType;

  bool fSaveNeighbourADCs;
  
  TTree* fSNAnaTree;

  int Run;
  int SubRun;
  int Event;
  
  int NTotHit   ;
  int NColHit   ;
  int NIndHit   ;
  int NHitNoBT  ;
  int NFlash    ;
  int NFlashNoBT;

  std::vector<int>   Hit_View              ;
  std::vector<int>   Hit_Size              ;
  std::vector<int>   Hit_TPC               ;
  std::vector<int>   Hit_Chan              ;
  std::vector<float> Hit_Time              ;
  std::vector<float> Hit_RMS               ;
  std::vector<float> Hit_SADC              ;
  std::vector<float> Hit_Int               ;
  std::vector<float> Hit_Peak              ;
  std::vector<int>   Hit_True_GenType      ;
  std::vector<int>   Hit_True_MainTrID     ;
  std::vector<float> Hit_True_EvEnergy     ;
  std::vector<float> Hit_True_X            ;                  
  std::vector<float> Hit_True_Y            ;                  
  std::vector<float> Hit_True_Z            ;                  
  std::vector<float> Hit_True_Energy       ;             
  std::vector<float> Hit_True_nElec        ;
  std::vector<int>   Hit_True_nIDEs        ; 
  std::vector<int>   Hit_AdjM5SADC         ;
  std::vector<int>   Hit_AdjM2SADC         ;
  std::vector<int>   Hit_AdjM1SADC         ;
  std::vector<int>   Hit_AdjP1SADC         ;
  std::vector<int>   Hit_AdjP2SADC         ;
  std::vector<int>   Hit_AdjP5SADC         ;
  std::vector<int>   Hit_AdjM5Chan         ;
  std::vector<int>   Hit_AdjM2Chan         ;
  std::vector<int>   Hit_AdjM1Chan         ;
  std::vector<int>   Hit_AdjP1Chan         ;
  std::vector<int>   Hit_AdjP2Chan         ;
  std::vector<int>   Hit_AdjP5Chan         ;
  std::vector<int>   PDS_Hit_FlashID       ;
  std::vector<float> PDS_Hit_YCenter       ;
  std::vector<float> PDS_Hit_ZCenter       ;
  std::vector<float> PDS_Hit_YWidth        ;
  std::vector<float> PDS_Hit_ZWidth        ;
  std::vector<float> PDS_Hit_Time          ;
  std::vector<float> PDS_Hit_TimeWidth     ;
  std::vector<float> PDS_Hit_TotalPE       ;
  std::vector<float> PDS_Hit_True_Purity   ;
  std::vector<float> PDS_Hit_Distance      ;
  std::vector<int>   True_VertexChan       ;
  std::vector<int>   True_Nu_Type          ;
  std::vector<int>   True_Nu_Lep_Type      ;
  std::vector<int>   True_Mode             ;
  std::vector<int>   True_CCNC             ;
  std::vector<int>   True_HitNucleon       ;
  std::vector<int>   True_Target           ;
  std::vector<int>   True_MarlSample       ;
  std::vector<float> True_MarlTime         ;
  std::vector<float> True_MarlWeight       ;
  std::vector<float> True_ENu              ;
  std::vector<float> True_ENu_Lep          ;
  std::vector<float> True_VertX            ;
  std::vector<float> True_VertY            ;
  std::vector<float> True_VertZ            ;
  std::vector<float> True_VertexT          ;
  std::vector<float> True_Px               ;
  std::vector<float> True_Py               ;
  std::vector<float> True_Pz               ;
  std::vector<float> True_Dirx             ;
  std::vector<float> True_Diry             ;
  std::vector<float> True_Dirz             ;
  std::vector<float> True_Time             ;

  int   TotGen_Marl;
  int   TotGen_APA ;
  int   TotGen_CPA ;
  int   TotGen_Ar39;
  int   TotGen_Neut;
  int   TotGen_Kryp;
  int   TotGen_Plon;
  int   TotGen_Rdon;
  int   TotGen_Ar42;


  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;

};

SNAna::SNAna(fhicl::ParameterSet const & p):EDAnalyzer(p)
{
  this->reconfigure(p);
}


void SNAna::reconfigure(fhicl::ParameterSet const & p)
{
  fRawDigitLabel      = p.get<std::string>("RawDigitLabel"     );  
  fHitLabel           = p.get<std::string>("HitLabel"          );
  fCalDataModuleLabel = p.get<std::string>("CalDataModuleLabel");
  fOpFlashModuleLabel = p.get<std::string>("OpFlashModuleLabel");
  fOpHitModuleLabel   = p.get<std::string>("OpHitModuleLabel"  );

  fGEANTLabel = p.get<std::string> ("GEANT4Label"  );
  fMARLLabel  = p.get<std::string> ("MARLEYLabel"  );
  fAPALabel   = p.get<std::string> ("APALabel"     );
  fCPALabel   = p.get<std::string> ("CPALabel"     );
  fAr39Label  = p.get<std::string> ("Argon39Label" );
  fNeutLabel  = p.get<std::string> ("NeutronLabel" );
  fKrypLabel  = p.get<std::string> ("KryptonLabel" );
  fPlonLabel  = p.get<std::string> ("PoloniumLabel");
  fRdonLabel  = p.get<std::string> ("RadonLabel"   );
  fAr42Label  = p.get<std::string> ("Argon42Label" );

  fSaveNeighbourADCs = p.get<bool> ("SaveNeighbourADCs");
  
} 


void SNAna::ResetVariables()
{
  trkIDToPType.clear();

  MarlParts.clear(); APAParts .clear(); CPAParts .clear(); Ar39Parts.clear();
  NeutParts.clear(); KrypParts.clear(); PlonParts.clear(); RdonParts.clear();
  Ar42Parts.clear();

  Run = SubRun = Event = -1;
  
  TotGen_Marl = TotGen_APA  = TotGen_CPA  = TotGen_Ar39 = 0;
  TotGen_Neut = TotGen_Kryp = TotGen_Plon = TotGen_Rdon = 0;
  TotGen_Ar42 = 0;

  NTotHit    = 0;
  NColHit    = 0;
  NIndHit    = 0;
  NHitNoBT   = 0;
  NFlash     = 0;
  NFlashNoBT = 0;

  Hit_View              .clear();
  Hit_Size              .clear();
  Hit_TPC               .clear();
  Hit_Chan              .clear();
  Hit_Time              .clear();
  Hit_RMS               .clear();
  Hit_SADC              .clear();
  Hit_Int               .clear();
  Hit_Peak              .clear();
  Hit_True_GenType      .clear();
  Hit_True_MainTrID     .clear();
  Hit_True_EvEnergy     .clear();
  Hit_True_X            .clear();
  Hit_True_Y            .clear();
  Hit_True_Z            .clear();
  Hit_True_Energy       .clear();
  Hit_True_nElec        .clear();
  Hit_True_nIDEs        .clear();
  Hit_AdjM5SADC         .clear();
  Hit_AdjM2SADC         .clear();
  Hit_AdjM1SADC         .clear();
  Hit_AdjP1SADC         .clear();
  Hit_AdjP2SADC         .clear();
  Hit_AdjP5SADC         .clear();
  Hit_AdjM5Chan         .clear();                  
  Hit_AdjM2Chan         .clear();                  
  Hit_AdjM1Chan         .clear();                  
  Hit_AdjP1Chan         .clear();             
  Hit_AdjP2Chan         .clear();
  Hit_AdjP5Chan         .clear(); 
  PDS_Hit_FlashID       .clear();
  PDS_Hit_YCenter       .clear();
  PDS_Hit_ZCenter       .clear();
  PDS_Hit_YWidth        .clear();
  PDS_Hit_ZWidth        .clear();
  PDS_Hit_Time          .clear();
  PDS_Hit_TimeWidth     .clear();
  PDS_Hit_TotalPE       .clear();
  PDS_Hit_True_Purity   .clear();
  PDS_Hit_Distance      .clear();
  True_VertexChan       .clear();
  True_Nu_Type          .clear();
  True_Nu_Lep_Type      .clear();
  True_Mode             .clear();
  True_CCNC             .clear();
  True_HitNucleon       .clear();
  True_Target           .clear();
  True_MarlSample       .clear();
  True_MarlTime         .clear();
  True_MarlWeight       .clear();
  True_ENu              .clear();
  True_ENu_Lep          .clear();
  True_VertX            .clear();
  True_VertY            .clear();
  True_VertZ            .clear();
  True_VertexT          .clear();
  True_Px               .clear();
  True_Py               .clear();
  True_Pz               .clear();
  True_Dirx             .clear();
  True_Diry             .clear();
  True_Dirz             .clear();
  True_Time             .clear();

  
}


void SNAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;


  fSNAnaTree = tfs->make<TTree>("SNSimTree","SN simulation analysis tree");
  //fSNAnaTree = tfs->make<TTree>("SNSimTree","SN simulation analysis tree");

  fSNAnaTree->Branch("Run"       , &Run       , "Run/I"       );
  fSNAnaTree->Branch("SubRun"    , &SubRun    , "SubRun/I"    );
  fSNAnaTree->Branch("Event"     , &Event     , "Event/I"     );
  fSNAnaTree->Branch("NTotHit"   , &NTotHit   , "NTotHits/I"  );
  fSNAnaTree->Branch("NColHit"   , &NColHit   , "NColHits/I"  );
  fSNAnaTree->Branch("NIndHit"   , &NIndHit   , "NIndHits/I"  );
  fSNAnaTree->Branch("NHitNoBT"  , &NHitNoBT  , "NHitNoBT/I"  );
  fSNAnaTree->Branch("NFlash"    , &NFlash    , "NFlash/I"    );
  fSNAnaTree->Branch("NFlashNoBT", &NFlashNoBT, "NFlashNoBT/I");

  fSNAnaTree->Branch("Hit_View"              , &Hit_View              );
  fSNAnaTree->Branch("Hit_Size"              , &Hit_Size              );
  fSNAnaTree->Branch("Hit_TPC"               , &Hit_TPC               );
  fSNAnaTree->Branch("Hit_Chan"              , &Hit_Chan              );
  fSNAnaTree->Branch("Hit_Time"              , &Hit_Time              );
  fSNAnaTree->Branch("Hit_RMS"               , &Hit_RMS               );
  fSNAnaTree->Branch("Hit_SADC"              , &Hit_SADC              );
  fSNAnaTree->Branch("Hit_Int"               , &Hit_Int               );
  fSNAnaTree->Branch("Hit_Peak"              , &Hit_Peak              );
  fSNAnaTree->Branch("Hit_True_GenType"      , &Hit_True_GenType      );
  fSNAnaTree->Branch("Hit_True_MainTrID"     , &Hit_True_MainTrID     );
  fSNAnaTree->Branch("Hit_True_EvEnergy"     , &Hit_True_EvEnergy     );
  fSNAnaTree->Branch("Hit_True_X"            , &Hit_True_X            );
  fSNAnaTree->Branch("Hit_True_Y"            , &Hit_True_Y            );
  fSNAnaTree->Branch("Hit_True_Z"            , &Hit_True_Z            );
  fSNAnaTree->Branch("Hit_True_Energy"       , &Hit_True_Energy       );
  fSNAnaTree->Branch("Hit_True_nElec"        , &Hit_True_nElec        );
  fSNAnaTree->Branch("Hit_True_nIDEs"        , &Hit_True_nIDEs        );
  fSNAnaTree->Branch("Hit_AdjM5SADC"         , &Hit_AdjM5SADC         );
  fSNAnaTree->Branch("Hit_AdjM2SADC"         , &Hit_AdjM2SADC         );
  fSNAnaTree->Branch("Hit_AdjM1SADC"         , &Hit_AdjM1SADC         );
  fSNAnaTree->Branch("Hit_AdjP1SADC"         , &Hit_AdjP1SADC         );
  fSNAnaTree->Branch("Hit_AdjP2SADC"         , &Hit_AdjP2SADC         );
  fSNAnaTree->Branch("Hit_AdjP5SADC"         , &Hit_AdjP5SADC         );
  fSNAnaTree->Branch("Hit_AdjM5Chan"         , &Hit_AdjM5Chan         );
  fSNAnaTree->Branch("Hit_AdjM2Chan"         , &Hit_AdjM2Chan         );
  fSNAnaTree->Branch("Hit_AdjM1Chan"         , &Hit_AdjM1Chan         );
  fSNAnaTree->Branch("Hit_AdjP1Chan"         , &Hit_AdjP1Chan         );
  fSNAnaTree->Branch("Hit_AdjP2Chan"         , &Hit_AdjP2Chan         );
  fSNAnaTree->Branch("Hit_AdjP5Chan"         , &Hit_AdjP5Chan         );
  fSNAnaTree->Branch("PDS_Hit_FlashID"       , &PDS_Hit_FlashID       );
  fSNAnaTree->Branch("PDS_Hit_YCenter"       , &PDS_Hit_YCenter       );
  fSNAnaTree->Branch("PDS_Hit_ZCenter"       , &PDS_Hit_ZCenter       );
  fSNAnaTree->Branch("PDS_Hit_YWidth"        , &PDS_Hit_YWidth        );
  fSNAnaTree->Branch("PDS_Hit_ZWidth"        , &PDS_Hit_ZWidth        );
  fSNAnaTree->Branch("PDS_Hit_Time"          , &PDS_Hit_Time          );
  fSNAnaTree->Branch("PDS_Hit_TimeWidth"     , &PDS_Hit_TimeWidth     );
  fSNAnaTree->Branch("PDS_Hit_TotalPE"       , &PDS_Hit_TotalPE       );
  fSNAnaTree->Branch("PDS_Hit_True_Purity"   , &PDS_Hit_True_Purity   );
  fSNAnaTree->Branch("PDS_Hit_Distance"      , &PDS_Hit_Distance      );
  
  fSNAnaTree->Branch("TotGen_Marl", &TotGen_Marl, "TotGen_Marl/I");
  fSNAnaTree->Branch("TotGen_APA" , &TotGen_APA , "TotGen_APA/I" );
  fSNAnaTree->Branch("TotGen_CPA" , &TotGen_CPA , "TotGen_CPA/I" );
  fSNAnaTree->Branch("TotGen_Ar39", &TotGen_Ar39, "TotGen_Ar39/I");
  fSNAnaTree->Branch("TotGen_Neut", &TotGen_Neut, "TotGen_Neut/I");
  fSNAnaTree->Branch("TotGen_Kryp", &TotGen_Kryp, "TotGen_Kryp/I");
  fSNAnaTree->Branch("TotGen_Plon", &TotGen_Plon, "TotGen_Plon/I");
  fSNAnaTree->Branch("TotGen_Rdon", &TotGen_Rdon, "TotGen_Rdon/I");
  fSNAnaTree->Branch("TotGen_Ar42", &TotGen_Ar42, "TotGen_Ar42/I");

  fSNAnaTree->Branch("True_VertexChan" , &True_VertexChan );
  fSNAnaTree->Branch("True_Nu_Type"    , &True_Nu_Type    );
  fSNAnaTree->Branch("True_Nu_Lep_Type", &True_Nu_Lep_Type);
  fSNAnaTree->Branch("True_Mode"       , &True_Mode       );
  fSNAnaTree->Branch("True_CCNC"       , &True_CCNC       );
  fSNAnaTree->Branch("True_HitNucleon" , &True_HitNucleon );
  fSNAnaTree->Branch("True_Target"     , &True_Target     );
  fSNAnaTree->Branch("True_MarlSample" , &True_MarlSample );
  fSNAnaTree->Branch("True_MarlTime"   , &True_MarlTime   );
  fSNAnaTree->Branch("True_MarlWeight" , &True_MarlWeight );
  fSNAnaTree->Branch("True_ENu"        , &True_ENu        );
  fSNAnaTree->Branch("True_ENu_Lep"    , &True_ENu_Lep    );
  fSNAnaTree->Branch("True_VertX"      , &True_VertX      );
  fSNAnaTree->Branch("True_VertY"      , &True_VertY      );
  fSNAnaTree->Branch("True_VertZ"      , &True_VertZ      );
  fSNAnaTree->Branch("True_VertexT"    , &True_VertexT    );
  fSNAnaTree->Branch("True_Px"         , &True_Px         );
  fSNAnaTree->Branch("True_Py"         , &True_Py         );
  fSNAnaTree->Branch("True_Pz"         , &True_Pz         );
  fSNAnaTree->Branch("True_Dirx"       , &True_Dirx       );
  fSNAnaTree->Branch("True_Diry"       , &True_Diry       );
  fSNAnaTree->Branch("True_Dirz"       , &True_Dirz       );
  fSNAnaTree->Branch("True_Time"       , &True_Time       );


  
} 


void SNAna::analyze(art::Event const & evt)
{
  ResetVariables();
 
  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  //GET THE RECO HITS.
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);

  //LIFT OUT THE MARLEY PARTICLES.
  auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fMARLLabel);
  art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue,evt,fGEANTLabel);
  FillMyMaps( MarlParts, MarlAssn, MarlTrue );
  TotGen_Marl = MarlParts.size();

  //SUPERNOVA TRUTH.
  art::FindManyP<sim::SupernovaTruth> SNTruth(MarlTrue, evt, fMARLLabel);
  std::set<int> signal_trackids;

  double Px_(0), Py_(0), Pz_(0), Pnorm(1);
  for (size_t i=0; i<MarlAssn.size(); ++i) {
    auto parts = MarlAssn.at(i);
    for (auto part = parts.begin(); part != parts.end(); part++) {
      signal_trackids.emplace((*part)->TrackId());
    }
  }

  for(size_t i = 0; i < MarlTrue->size(); i++)
  {
    True_Nu_Type    .push_back(MarlTrue->at(i).GetNeutrino().Nu().PdgCode());
    True_Nu_Lep_Type.push_back(MarlTrue->at(i).GetNeutrino().Lepton().PdgCode()); 
    True_Mode       .push_back(MarlTrue->at(i).GetNeutrino().Mode());
    True_CCNC       .push_back(MarlTrue->at(i).GetNeutrino().CCNC());
    True_Target     .push_back(MarlTrue->at(i).GetNeutrino().Target());
    True_HitNucleon .push_back(MarlTrue->at(i).GetNeutrino().HitNuc());
    True_VertX      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vx()); 
    True_VertY      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vy()); 
    True_VertZ      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vz()); 
    True_ENu_Lep    .push_back(MarlTrue->at(i).GetNeutrino().Lepton().E());
    True_ENu        .push_back(MarlTrue->at(i).GetNeutrino().Nu().E());
    True_Px         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Px());
    True_Py         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Py());
    True_Pz         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Pz());
    Pnorm = std::sqrt(Px_*Px_+Py_*Py_+Pz_*Pz_);
    double Px = Px_/Pnorm;
    double Py = Py_/Pnorm;
    double Pz = Pz_/Pnorm;
    True_Dirx       .push_back(Px);
    True_Diry       .push_back(Py);
    True_Dirz       .push_back(Pz);
    True_Time       .push_back(MarlTrue->at(i).GetNeutrino().Lepton().T());
    for (size_t j = 0; j < SNTruth.at(i).size(); j++) 
    {
      const sim::SupernovaTruth ThisTr = (*SNTruth.at(i).at(j));
      True_MarlTime  .push_back(ThisTr.SupernovaTime);
      True_MarlWeight.push_back(ThisTr.Weight);
      True_MarlSample.push_back(ThisTr.SamplingMode);
    }
  }

  for(size_t i=0; i<True_VertX.size(); ++i) {
    bool caught = false;
    double Vertex[3] = {True_VertX[i],
                        True_VertY[i],
                        True_VertZ[i]};
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
      True_VertexChan.push_back(-1); 
    }
    else
    {
      True_VertexChan.push_back(geo->PlaneWireToChannel(WireID));
    }

    //CM/MICROSECOND.
    double drift_velocity = detp->DriftVelocity(detp->Efield(),detp->Temperature());
    //CM/TICK
    drift_velocity = drift_velocity*0.5;
    True_VertexT.push_back(True_VertX.back()/drift_velocity);
  }

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

  std::vector<simb::MCParticle> allTruthParts;
  for(auto& it: APAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Ar39Parts)
    allTruthParts.push_back(it.second);
  for(auto& it: NeutParts)
    allTruthParts.push_back(it.second);
  for(auto& it: KrypParts)
    allTruthParts.push_back(it.second);
  for(auto& it: PlonParts)
    allTruthParts.push_back(it.second);
  for(auto& it: RdonParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Ar42Parts)
    allTruthParts.push_back(it.second);
  
  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
      {kMarl, MarlParts},
      {kAPA,  APAParts },
      {kCPA,  CPAParts },
      {kAr39, Ar39Parts},
      {kAr42, Ar42Parts},
      {kNeut, NeutParts},
      {kKryp, KrypParts},
      {kPlon, PlonParts},
      {kRdon, RdonParts}
  };

  for(auto const& it : PTypeToMap){
      const PType p=it.first;
      auto const& m=it.second;
      for(auto const& it2 : m){
          trkIDToPType.insert(std::make_pair(it2.first, p));
      }
  }
  
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

  NTotHit = reco_hits->size();
  int colHitCount(0);
  int LoopHits = NTotHit;
  std::cout << "---- There are " << NTotHit << " hits in the event." << std::endl;

  for(int hit = 0; hit < LoopHits; ++hit) 
  {
    recob::Hit const& ThisHit = reco_hits->at(hit);  
    if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) {
      ++NIndHit;
    } else { 
      ++NColHit;
    }
  }
  art::Handle< std::vector<recob::Wire> > wireVecHandle;
  art::Handle< std::vector<raw::RawDigit> > rawDigitsVecHandle;
  evt.getByLabel(fCalDataModuleLabel, wireVecHandle);
  evt.getByLabel(fRawDigitLabel, rawDigitsVecHandle);

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
    
    Hit_View.push_back(ThisHit.View());
    Hit_Size.push_back(ThisHit.EndTick() - ThisHit.StartTick());
    Hit_TPC .push_back(ThisHit.WireID().TPC);
    int channel = ThisHit.Channel();
    Hit_Chan.push_back(channel);
    
    Hit_AdjM5SADC.push_back(0);
    Hit_AdjM2SADC.push_back(0);
    Hit_AdjM1SADC.push_back(0);
    Hit_AdjP1SADC.push_back(0);
    Hit_AdjP2SADC.push_back(0);
    Hit_AdjP5SADC.push_back(0);
    Hit_AdjM5Chan.push_back(0);
    Hit_AdjM2Chan.push_back(0);
    Hit_AdjM1Chan.push_back(0);
    Hit_AdjP1Chan.push_back(0);
    Hit_AdjP2Chan.push_back(0);
    Hit_AdjP5Chan.push_back(0);
   
  
    std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
    for(size_t i=0; i<rawDigitsVecHandle->size(); ++i)
    {
      int rawWireChannel=(*rawDigitsVecHandle)[i].Channel();
      raw::RawDigit::ADCvector_t ADCs((*rawDigitsVecHandle)[i].Samples());
      std::vector<geo::WireID> adjacentwire = geo->ChannelToWire(rawWireChannel);

      if (adjacentwire.size() < 1) continue;
      if (adjacentwire[0].Plane == geo::kU ||
          adjacentwire[0].Plane == geo::kV)
        continue;
      
      switch(rawWireChannel-channel)
      {
      case -5:
        Hit_AdjM5SADC[colHitCount] = 0;
        Hit_AdjM5Chan[colHitCount] = rawWireChannel;
        raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                        (*rawDigitsVecHandle)[i].Compression());
        for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
          Hit_AdjM5SADC[colHitCount]+=ADCs[i];
        break;
      case -2:
        Hit_AdjM2SADC[colHitCount] = 0;
        Hit_AdjM2Chan[colHitCount] = rawWireChannel;
        raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                        (*rawDigitsVecHandle)[i].Compression());
        for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
          Hit_AdjM2SADC[colHitCount]+=ADCs[i];
        break;
      case -1:
        Hit_AdjM1SADC[colHitCount] = 0;
        Hit_AdjM1Chan[colHitCount] = rawWireChannel;
        raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                        (*rawDigitsVecHandle)[i].Compression());
        for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
          Hit_AdjM1SADC[colHitCount]+=ADCs[i];
        break;
      case  1:
        Hit_AdjP1SADC[colHitCount] = 0;
        Hit_AdjP1Chan[colHitCount] = rawWireChannel;
        raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                        (*rawDigitsVecHandle)[i].Compression());
        for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
          Hit_AdjP1SADC[colHitCount]+=ADCs[i];
        break;
      case  2:
        Hit_AdjP2SADC[colHitCount] = 0;
        Hit_AdjP2Chan[colHitCount] = rawWireChannel;
        raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                        (*rawDigitsVecHandle)[i].Compression());
        for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
          Hit_AdjP2SADC[colHitCount]+=ADCs[i];
        break;
      case  5:
        Hit_AdjP5SADC[colHitCount] = 0;
        Hit_AdjP5Chan[colHitCount] = rawWireChannel;
        raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                        (*rawDigitsVecHandle)[i].Compression());
        for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
          Hit_AdjP5SADC[colHitCount]+=ADCs[i];
        break;
      default:
        break;
      }
    }

    Hit_Time.push_back(ThisHit.PeakTime());
    Hit_RMS .push_back(ThisHit.RMS());
    Hit_SADC.push_back(ThisHit.SummedADC());
    Hit_Int .push_back(ThisHit.Integral());
    Hit_Peak.push_back(ThisHit.PeakAmplitude());
    Hit_True_nIDEs.push_back(ThisHitIDE.size());

    if(ThisHitIDE.size()==0)
      NHitNoBT++;

    //WHICH PARTICLE CONTRIBUTED MOST TO THE HIT.
    double TopEFrac = -DBL_MAX;

    Hit_True_EvEnergy.push_back(0);
    Hit_True_MainTrID.push_back(0);
    for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL) {
      for (size_t ipart=0; ipart<allTruthParts.size(); ++ipart) {

        if (allTruthParts[ipart].TrackId() == ThisHitIDE[ideL].trackID) {
          Hit_True_EvEnergy.at(colHitCount) += allTruthParts[ipart].E();
        }
      }
      if (ThisHitIDE[ideL].energyFrac > TopEFrac) {
        TopEFrac = ThisHitIDE[ideL].energyFrac;
        Hit_True_MainTrID.at(colHitCount) = ThisHitIDE[ideL].trackID;

      }
    }

    PType ThisPType      = WhichParType(Hit_True_MainTrID.at(colHitCount));
    Hit_True_GenType.push_back(ThisPType);

    if(Hit_True_MainTrID[colHitCount] == -1)
    {
      Hit_True_X     .push_back(-1);
      Hit_True_Y     .push_back(-1);
      Hit_True_Z     .push_back(-1);
      Hit_True_Energy.push_back(-1);
      Hit_True_nElec .push_back(-1);
    }
    else
    {
      for(unsigned int i = 0; i < ThisSimIDE.size(); i++)
      {
        if(ThisSimIDE.at(i)->trackID==Hit_True_MainTrID[colHitCount])
        {
          Hit_True_X     .push_back(ThisSimIDE.at(i)->x           );
          Hit_True_Y     .push_back(ThisSimIDE.at(i)->y           );
          Hit_True_Z     .push_back(ThisSimIDE.at(i)->z           );
          Hit_True_Energy.push_back(ThisSimIDE.at(i)->energy      );
          Hit_True_nElec .push_back(ThisSimIDE.at(i)->numElectrons);
          break;
        }
      }
    }
    
    if      (ThisPType == 0) { ColHits_Oth .push_back( ThisHit ); }
    else if (ThisPType == 1) { ColHits_Marl.push_back( ThisHit ); }
    else if (ThisPType == 2) { ColHits_APA .push_back( ThisHit ); }
    else if (ThisPType == 3) { ColHits_CPA .push_back( ThisHit ); }
    else if (ThisPType == 4) { ColHits_Ar39.push_back( ThisHit ); }
    else if (ThisPType == 5) { ColHits_Neut.push_back( ThisHit ); }
    else if (ThisPType == 6) { ColHits_Kryp.push_back( ThisHit ); }
    else if (ThisPType == 7) { ColHits_Plon.push_back( ThisHit ); }
    else if (ThisPType == 8) { ColHits_Rdon.push_back( ThisHit ); }
    else if (ThisPType == 9) { ColHits_Ar42.push_back( ThisHit ); }

    colHitCount++;
    }
  } 

  std::cerr << "\n\nAfter all of that I have a total of " << ColHits_Marl.size() << " MARLEY col plane hits." << std::endl;
  for (size_t hh=0; hh<ColHits_Marl.size(); ++hh) {
    std::cerr << "-- Hit " << hh << " chan " << ColHits_Marl[hh].Channel() << " at " << ColHits_Marl[hh].PeakTime() << std::endl;
  }

  // Get flashes from event
  art::Handle< std::vector< recob::OpFlash > > FlashHandle;
  std::vector<art::Ptr<recob::OpFlash> > flashlist;
  if (evt.getByLabel(fOpFlashModuleLabel, FlashHandle)) {
    art::fill_ptr_vector(flashlist, FlashHandle);
    std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
  }
  else {
    mf::LogError("FlashMatchAna") << "Cannot load any flashes. Failing";
    abort();
  }
  // std::cout << "Rebuilding event" << std::endl;
  // pbt_serv->Rebuild(evt);
  // std::cout << "Rebuilt event" << std::endl;
  // std::vector< art::Ptr< sim::OpDetBacktrackerRecord >> opdbtr = pbt_serv->OpDetBTRs();
  // std::cout << "OpDetBacktrackerRecord size " << opdbtr.size() << std::endl;
  // Get assosciations between flashes and hits
  art::FindManyP< recob::OpHit > Assns(flashlist, evt, fOpFlashModuleLabel);
  for(size_t i = 0; i < flashlist.size(); ++i)
  {
    // Get OpFlash and associated hits
    recob::OpFlash TheFlash = *flashlist[i];
    std::vector< art::Ptr<recob::OpHit> > matchedHits = Assns.at(i);
    //if(matchedHits.size() == 0) 
    // std::cout << "matchedHits " << matchedHits.size() << std::endl;
    // for (size_t hit=0; hit<matchedHits.size(); ++hit) {
    //   std::cout << "------   " << hit   << std::endl;
    //   std::cout << "OpChannel();   " << matchedHits[hit]->OpChannel()   << std::endl;
    //   std::cout << "PeakTimeAbs(); " << matchedHits[hit]->PeakTimeAbs() << std::endl;
    //   std::cout << "PeakTime();    " << matchedHits[hit]->PeakTime()    << std::endl; 
    //   std::cout << "Frame();       " << matchedHits[hit]->Frame()       << std::endl; 
    //   std::cout << "Width();       " << matchedHits[hit]->Width()       << std::endl;
    //   std::cout << "Area();        " << matchedHits[hit]->Area()        << std::endl;
    //   std::cout << "Amplitude();   " << matchedHits[hit]->Amplitude()   << std::endl;
    //   std::cout << "PE();          " << matchedHits[hit]->PE()          << std::endl;
    //   std::cout << "FastToTotal(); " << matchedHits[hit]->FastToTotal() << std::endl; 
    //   std::unordered_set<const sim::SDP*> spds = pbt_serv->OpHitToEveSimSDPs_Ps(matchedHits[hit]);
    //   std::cout << "spds.size() " << spds.size() << std::endl;
    // }

    double purity = pbt_serv->OpHitCollectionPurity(signal_trackids, matchedHits);
    if (purity!=0) {
      std::cout << "purity " << purity << std::endl;
    }
    // std::set<int> setseventid = pbt_serv->GetSetOfEveIds(matchedHits);
    // std::cout << "Size setseventid " << setseventid.size() << std::endl;
    // Put flash info into variables
    PDS_Hit_FlashID       .push_back(i);
    PDS_Hit_YCenter       .push_back(TheFlash.YCenter());
    PDS_Hit_ZCenter       .push_back(TheFlash.ZCenter());
    PDS_Hit_YWidth        .push_back(TheFlash.YWidth());
    PDS_Hit_ZWidth        .push_back(TheFlash.ZWidth());
    PDS_Hit_Time          .push_back(TheFlash.Time());
    PDS_Hit_TimeWidth     .push_back(TheFlash.TimeWidth());
    PDS_Hit_TotalPE       .push_back(TheFlash.TotalPE());
    PDS_Hit_True_Purity   .push_back(purity);
    PDS_Hit_Distance      .push_back(TMath::Sqrt(TMath::Power(True_VertY.back()-TheFlash.YCenter(),2)+
                                                 TMath::Power(True_VertZ.back()-TheFlash.ZCenter(),2)));
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
  PType ThisPType = kUnknown;
  auto const& it=trkIDToPType.find(TrID);
  if(it!=trkIDToPType.end()){
    ThisPType=it->second;
  }
  return ThisPType;

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
