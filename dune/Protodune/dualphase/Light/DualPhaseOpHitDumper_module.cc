#include <functional>   // std::greater

//ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TGeoMatrix.h" // TGeoHMatrix

// LArSoft libraries

//LArSoft includes
#include "larcore/Geometry/Geometry.h"

#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
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
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

enum PType{kUnknown=0, kSign, kAr39, kNeut, kKryp, kRdon , kAr42};
std::map<PType, std::string> PTypeString{
  {kUnknown,"Unknown"},
  {kSign   ,"Sign"   },
  {kAr39   ,"Ar39"   },
  {kNeut   ,"Neut"   },
  {kKryp   ,"Kryp"   },
  {kRdon   ,"Rdon"   },
  {kAr42   ,"Ar42"   }};


class DualPhaseOpHitDumper : public art::EDAnalyzer {

public:
  explicit DualPhaseOpHitDumper(fhicl::ParameterSet const & p);

  DualPhaseOpHitDumper(DualPhaseOpHitDumper const &) = delete;
  DualPhaseOpHitDumper(DualPhaseOpHitDumper &&) = delete;
  DualPhaseOpHitDumper & operator = (DualPhaseOpHitDumper const &) = delete;
  DualPhaseOpHitDumper & operator = (DualPhaseOpHitDumper &&) = delete;

  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void endJob() override;

private:

  void ResetVariables();
  void FillMyMaps  ( std::map< int, simb::MCParticle> &MyMap,
                     art::FindManyP<simb::MCParticle> Assn,
                     art::Handle< std::vector<simb::MCTruth> > Hand,
                     std::map<int, int>* indexMap=nullptr);
  PType WhichParType( int TrID );
  bool  InMyMap     ( int TrID, std::map< int, simb::MCParticle> ParMap );
  void FillTruth(const art::FindManyP<simb::MCParticle> Assn,
                 const art::Handle<std::vector<simb::MCTruth>>& Hand,
                 const PType type) ;
  void SaveNeighbourADC(int channel,
                        art::Handle< std::vector<raw::RawDigit> >rawDigitsVecHandle,
                        std::set<int> badChannels,
                        recob::Hit const& hits);

  

  int firstCatch;
  int secondCatch;
  int thirdCatch;

  std::string fOpHitModuleLabel;

  std::string fGEANTLabel;
  std::string fSIGNLabel ; std::map< int, simb::MCParticle > SignParts;
  std::string fAr39Label ; std::map< int, simb::MCParticle > Ar39Parts;
  std::string fNeutLabel ; std::map< int, simb::MCParticle > NeutParts;
  std::string fKrypLabel ; std::map< int, simb::MCParticle > KrypParts;
  std::string fRdonLabel ; std::map< int, simb::MCParticle > RdonParts;
  std::string fAr42Label ; std::map< int, simb::MCParticle > Ar42Parts;
  std::map<int, const simb::MCParticle*> truthmap;
  // Which SignEY interaction (if any) caused this true track ID?
  std::map<int, int> trkIDToMarleyIndex;

  // Mapping from track ID to particle type, for use in WhichParType()
  std::map<int, PType> trkIDToPType;
  
  bool fSaveTruth;
  bool fSavePDS;

  TTree* fAnaTree;

  int Run;
  int SubRun;
  int Event;

  std::vector<int>                  PDS_OpHit_OpChannel      ;
  std::vector<double>               PDS_OpHit_X              ;
  std::vector<double>               PDS_OpHit_Y              ;
  std::vector<double>               PDS_OpHit_Z              ;
  std::vector<double>               PDS_OpHit_PeakTimeAbs    ;
  std::vector<double>               PDS_OpHit_PeakTime       ;
  std::vector<unsigned short>       PDS_OpHit_Frame          ;
  std::vector<double>               PDS_OpHit_Width          ;
  std::vector<double>               PDS_OpHit_Area           ;
  std::vector<double>               PDS_OpHit_Amplitude      ;
  std::vector<double>               PDS_OpHit_PE             ;
  std::vector<double>               PDS_OpHit_FastToTotal    ;
  std::vector<int>                  PDS_OpHit_True_GenType   ;
  std::vector<int>                  PDS_OpHit_True_Index     ;
  std::vector<double>               PDS_OpHit_True_Energy    ;
  std::vector<int>                  PDS_OpHit_True_TrackID   ;
  std::vector<int>                  PDS_OpHit_True_GenTypeAll;
  std::vector<double>               PDS_OpHit_True_EnergyAll ;
  std::vector<int>                  PDS_OpHit_True_TrackIDAll;
  std::vector<int>                  PDS_OpHit_True_IndexAll  ;

  std::vector<int>                  True_Bck_Mode            ;
  std::vector<double>               True_Bck_VertX           ;
  std::vector<double>               True_Bck_VertY           ;
  std::vector<double>               True_Bck_VertZ           ;
  std::vector<double>               True_Bck_Time            ;
  std::vector<double>               True_Bck_Energy          ;
  std::vector<int>                  True_Bck_PDG             ;
  std::vector<int>                  True_Bck_ID              ;

  int   TotGen_Sign;
  int   TotGen_Ar39;
  int   TotGen_Neut;
  int   TotGen_Kryp;
  int   TotGen_Rdon;
  int   TotGen_Ar42;


  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;

};

DualPhaseOpHitDumper::DualPhaseOpHitDumper(fhicl::ParameterSet const & p):EDAnalyzer(p)
{
  this->reconfigure(p);
}


void DualPhaseOpHitDumper::reconfigure(fhicl::ParameterSet const & p)
{
  fOpHitModuleLabel   = p.get<std::string>("OpHitModuleLabel"  );



  fGEANTLabel = p.get<std::string> ("GEANT4Label"  );
  fSIGNLabel  = p.get<std::string> ("SIGNALLabel"  );
  fAr39Label  = p.get<std::string> ("Argon39Label" );
  fNeutLabel  = p.get<std::string> ("NeutronLabel" );
  fKrypLabel  = p.get<std::string> ("KryptonLabel" );
  fRdonLabel  = p.get<std::string> ("RadonLabel"   );
  fAr42Label  = p.get<std::string> ("Argon42Label" );

  fSaveTruth = p.get<bool>("SaveTruth",0);
  fSavePDS   = p.get<bool>("SavePDS",1);

  std::cout << "Reconfigured " << this->processName() << " with "
            << " SaveNeighbourADCs: " << std::boolalpha
            << " SaveTruth: " << std::boolalpha << fSaveTruth
            << " SaveIDEs: " << std::boolalpha << std::endl;
}


void DualPhaseOpHitDumper::ResetVariables()
{
  trkIDToPType.clear();

  SignParts.clear();  Ar39Parts.clear();
  NeutParts.clear(); KrypParts.clear(); RdonParts.clear();
  Ar42Parts.clear();

  Run = SubRun = Event = -1;

  TotGen_Sign = TotGen_Ar39 = 0;
  TotGen_Neut = TotGen_Kryp = TotGen_Rdon = 0;
  TotGen_Ar42 = 0;


  PDS_OpHit_OpChannel      .clear();
  PDS_OpHit_X              .clear();
  PDS_OpHit_Y              .clear();
  PDS_OpHit_Z              .clear();
  PDS_OpHit_PeakTimeAbs    .clear();
  PDS_OpHit_PeakTime       .clear();
  PDS_OpHit_Frame          .clear();
  PDS_OpHit_Width          .clear();
  PDS_OpHit_Area           .clear();
  PDS_OpHit_Amplitude      .clear();
  PDS_OpHit_PE             .clear();
  PDS_OpHit_FastToTotal    .clear();
  PDS_OpHit_True_GenType   .clear();
  PDS_OpHit_True_Index     .clear();
  PDS_OpHit_True_Energy    .clear();
  PDS_OpHit_True_TrackID   .clear();
  PDS_OpHit_True_GenTypeAll.clear();
  PDS_OpHit_True_EnergyAll .clear();
  PDS_OpHit_True_TrackIDAll.clear();
  PDS_OpHit_True_IndexAll  .clear();

  True_Bck_Mode            .clear();
  True_Bck_VertX           .clear();
  True_Bck_VertY           .clear();
  True_Bck_VertZ           .clear();
  True_Bck_Time            .clear();
  True_Bck_Energy          .clear();
  True_Bck_PDG             .clear();
  True_Bck_ID              .clear();


}


void DualPhaseOpHitDumper::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;


  fAnaTree = tfs->make<TTree>("OpHitDumpTree","OpHitDumper analysis tree");

  fAnaTree->Branch("Run"       , &Run       , "Run/I"       );
  fAnaTree->Branch("SubRun"    , &SubRun    , "SubRun/I"    );
  fAnaTree->Branch("Event"     , &Event     , "Event/I"     );


  if (fSavePDS) {
    fAnaTree->Branch("PDS_OpHit_OpChannel"      , &PDS_OpHit_OpChannel      );
    fAnaTree->Branch("PDS_OpHit_X"              , &PDS_OpHit_X              );
    fAnaTree->Branch("PDS_OpHit_Y"              , &PDS_OpHit_Y              );
    fAnaTree->Branch("PDS_OpHit_Z"              , &PDS_OpHit_Z              );
    fAnaTree->Branch("PDS_OpHit_PeakTimeAbs"    , &PDS_OpHit_PeakTimeAbs    );
    fAnaTree->Branch("PDS_OpHit_PeakTime"       , &PDS_OpHit_PeakTime       );
    fAnaTree->Branch("PDS_OpHit_Frame"          , &PDS_OpHit_Frame          );
    fAnaTree->Branch("PDS_OpHit_Width"          , &PDS_OpHit_Width          );
    fAnaTree->Branch("PDS_OpHit_Area"           , &PDS_OpHit_Area           );
    fAnaTree->Branch("PDS_OpHit_Amplitude"      , &PDS_OpHit_Amplitude      );
    fAnaTree->Branch("PDS_OpHit_PE"             , &PDS_OpHit_PE             );
    fAnaTree->Branch("PDS_OpHit_FastToTotal"    , &PDS_OpHit_FastToTotal    );
    fAnaTree->Branch("PDS_OpHit_True_GenType"   , &PDS_OpHit_True_GenType   );
    fAnaTree->Branch("PDS_OpHit_True_Energy"    , &PDS_OpHit_True_Energy    );
    fAnaTree->Branch("PDS_OpHit_True_TrackID"   , &PDS_OpHit_True_TrackID   );
    fAnaTree->Branch("PDS_OpHit_True_GenTypeAll", &PDS_OpHit_True_GenTypeAll);
    fAnaTree->Branch("PDS_OpHit_True_EnergyAll" , &PDS_OpHit_True_EnergyAll );
    fAnaTree->Branch("PDS_OpHit_True_TrackIDAll", &PDS_OpHit_True_TrackIDAll);
    fAnaTree->Branch("PDS_OpHit_True_IndexAll"  , &PDS_OpHit_True_IndexAll  );
  }

  fAnaTree->Branch("True_Bck_Mode"            , &True_Bck_Mode            );
  fAnaTree->Branch("True_Bck_VertX"           , &True_Bck_VertX           );
  fAnaTree->Branch("True_Bck_VertY"           , &True_Bck_VertY           );
  fAnaTree->Branch("True_Bck_VertZ"           , &True_Bck_VertZ           );
  fAnaTree->Branch("True_Bck_Time"            , &True_Bck_Time            );
  fAnaTree->Branch("True_Bck_Energy"          , &True_Bck_Energy          );
  fAnaTree->Branch("True_Bck_PDG"             , &True_Bck_PDG             );
  fAnaTree->Branch("True_Bck_ID"              , &True_Bck_ID              );



  fAnaTree->Branch("TotGen_Sign", &TotGen_Sign, "TotGen_Sign/I");
  fAnaTree->Branch("TotGen_Ar39", &TotGen_Ar39, "TotGen_Ar39/I");
  fAnaTree->Branch("TotGen_Neut", &TotGen_Neut, "TotGen_Neut/I");
  fAnaTree->Branch("TotGen_Kryp", &TotGen_Kryp, "TotGen_Kryp/I");
  fAnaTree->Branch("TotGen_Rdon", &TotGen_Rdon, "TotGen_Rdon/I");
  fAnaTree->Branch("TotGen_Ar42", &TotGen_Ar42, "TotGen_Ar42/I");


}


void DualPhaseOpHitDumper::analyze(art::Event const & evt)
{
  ResetVariables();

  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
  auto const* geo = lar::providerFrom<geo::Geometry>();

  //LIFT OUT THE SignEY PARTICLES.
  //auto SignTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fSIGNLabel);

  auto SignTrue = evt.getHandle< std::vector<simb::MCTruth> >(fSIGNLabel);
   if (SignTrue) {
    art::FindManyP<simb::MCParticle> SignAssn(SignTrue,evt,fGEANTLabel);

    FillMyMaps( SignParts, SignAssn, SignTrue );
    TotGen_Sign = SignParts.size();
    std::cout << "Found " << TotGen_Sign << " signal tracks" << std::endl;
    if(fSaveTruth) FillTruth(SignAssn , SignTrue , kSign );
  }

  auto Ar39True = evt.getHandle< std::vector<simb::MCTruth> >(fAr39Label);
  if (Ar39True) {
    art::FindManyP<simb::MCParticle> Ar39Assn(Ar39True,evt,fGEANTLabel);
    FillMyMaps( Ar39Parts, Ar39Assn, Ar39True );
    TotGen_Ar39 = Ar39Parts.size();
    if(fSaveTruth) FillTruth(Ar39Assn, Ar39True, kAr39);
  }

  auto NeutTrue = evt.getHandle< std::vector<simb::MCTruth> >(fNeutLabel);
  if (NeutTrue) {
    art::FindManyP<simb::MCParticle> NeutAssn(NeutTrue,evt,fGEANTLabel);
    FillMyMaps( NeutParts, NeutAssn, NeutTrue );
    TotGen_Neut = NeutParts.size();
    if(fSaveTruth) FillTruth(NeutAssn, NeutTrue, kNeut);
  }

  auto KrypTrue = evt.getHandle< std::vector<simb::MCTruth> >(fKrypLabel);
  if (KrypTrue) {
    art::FindManyP<simb::MCParticle> KrypAssn(KrypTrue,evt,fGEANTLabel);
    FillMyMaps( KrypParts, KrypAssn, KrypTrue );
    TotGen_Kryp = KrypParts.size();
    if(fSaveTruth) FillTruth(KrypAssn, KrypTrue, kKryp);
  }

  auto RdonTrue = evt.getHandle< std::vector<simb::MCTruth> >(fRdonLabel);
  if (RdonTrue) {
    art::FindManyP<simb::MCParticle> RdonAssn(RdonTrue,evt,fGEANTLabel);
    FillMyMaps( RdonParts, RdonAssn, RdonTrue );
    TotGen_Rdon = RdonParts.size();
    if(fSaveTruth) FillTruth(RdonAssn, RdonTrue, kRdon);
  }

  auto Ar42True = evt.getHandle< std::vector<simb::MCTruth> >(fAr42Label);
  if (Ar42True) {
    art::FindManyP<simb::MCParticle> Ar42Assn(Ar42True,evt,fGEANTLabel);
    FillMyMaps( Ar42Parts, Ar42Assn, Ar42True );
    TotGen_Ar42 = Ar42Parts.size();
    if(fSaveTruth) FillTruth(Ar42Assn, Ar42True, kAr42);
  }

  std::vector<simb::MCParticle> allTruthParts;
  for(auto& it: SignParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Ar39Parts)
    allTruthParts.push_back(it.second);
  for(auto& it: NeutParts)
    allTruthParts.push_back(it.second);
  for(auto& it: KrypParts)
    allTruthParts.push_back(it.second);
  for(auto& it: RdonParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Ar42Parts)
    allTruthParts.push_back(it.second);

  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
    {kSign, SignParts},
    {kAr39, Ar39Parts},
    {kAr42, Ar42Parts},
    {kNeut, NeutParts},
    {kKryp, KrypParts},
    {kRdon, RdonParts}
  };

  std::map<PType,std::set<int>> PTypeToTrackID;

  for(auto const& it : PTypeToMap){
    const PType p=it.first;
    auto const& m=it.second;
    for(auto const& it2 : m){
      trkIDToPType.insert(std::make_pair(it2.first, p)); //std::cout << "Inserting trackID " << it2.first <<" to PType "  << p << " inn trkIDToPType map. " << std::endl;
      PTypeToTrackID[p].insert(it2.first);
    }
  }

  std::cout << "THE EVENTS NUMBER IS: " << Event << std::endl;

  if (fSavePDS) {

    std::vector<art::Ptr<recob::OpHit> > ophitlist;
    std::map<PType, std::vector<art::Ptr<recob::OpHit> > > map_of_ophit;

    auto OpHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fOpHitModuleLabel);
    if (OpHitHandle) {
      art::fill_ptr_vector(ophitlist, OpHitHandle);

      std::cout << "There are " << ophitlist.size() << " optical hits in the event." << std::endl;

      for(size_t i = 0; i < ophitlist.size(); ++i)
      {
        
        std::vector<sim::TrackSDP> vec_tracksdp = pbt_serv->OpHitToTrackSDPs(ophitlist.at(i));
        PType gen = kUnknown;
//        std::cout << "Analyizing hit " << i << " " << vec_tracksdp.size()<<  std::endl;
        

        std::sort(vec_tracksdp.begin(), vec_tracksdp.end(),
                  [](const sim::TrackSDP& a, const sim::TrackSDP& b) -> bool { return a.energyFrac > b.energyFrac; });

        for (size_t iSDP=0; iSDP<vec_tracksdp.size(); ++iSDP) {
          PDS_OpHit_True_TrackIDAll.push_back(vec_tracksdp[iSDP].trackID);
          PDS_OpHit_True_GenTypeAll.push_back(WhichParType(vec_tracksdp[iSDP].trackID));
          PDS_OpHit_True_EnergyAll .push_back(vec_tracksdp[iSDP].energy);
          PDS_OpHit_True_IndexAll  .push_back((int)i);
        }

        if (vec_tracksdp.size()>0){

          //std::cout << "track assing to this hit " << vec_tracksdp[0].trackID << " " <<WhichParType(vec_tracksdp[0].trackID) << " " << pi_serv->TrackIdToParticle(vec_tracksdp[0].trackID).TrackId()<<  std::endl;
          //int MainTrID = vec_tracksdp[0].trackID;
          PDS_OpHit_True_TrackID.push_back(pi_serv->TrackIdToParticle(vec_tracksdp[0].trackID).TrackId());
          gen = WhichParType(pi_serv->TrackIdToParticle(vec_tracksdp[0].trackID).TrackId());
          PDS_OpHit_True_GenType.push_back(gen);
          PDS_OpHit_True_Energy .push_back(vec_tracksdp[0].energy);
        } else {

          //std::cout << "NO track assing to this hit!!!! " << std::endl;
          PDS_OpHit_True_TrackID.push_back(-1);
          PDS_OpHit_True_GenType.push_back(kUnknown);
          PDS_OpHit_True_Energy .push_back(-1);
        }

        map_of_ophit[gen].push_back(ophitlist.at(i));

        double xyz_optdet[3]={0,0,0};
        double xyz_world [3]={0,0,0};

        geo->OpDetGeoFromOpChannel(ophitlist[i]->OpChannel()).LocalToWorld(xyz_optdet,xyz_world);
        PDS_OpHit_OpChannel   .push_back(ophitlist[i]->OpChannel());
        PDS_OpHit_X           .push_back(xyz_world[0]);
        PDS_OpHit_Y           .push_back(xyz_world[1]);
        PDS_OpHit_Z           .push_back(xyz_world[2]);
        PDS_OpHit_PeakTimeAbs .push_back(ophitlist[i]->PeakTimeAbs());
        PDS_OpHit_PeakTime    .push_back(ophitlist[i]->PeakTime());
        PDS_OpHit_Frame       .push_back(ophitlist[i]->Frame());
        PDS_OpHit_Width       .push_back(ophitlist[i]->Width());
        PDS_OpHit_Area        .push_back(ophitlist[i]->Area());
        PDS_OpHit_Amplitude   .push_back(ophitlist[i]->Amplitude());
        PDS_OpHit_PE          .push_back(ophitlist[i]->PE());
        PDS_OpHit_FastToTotal .push_back(ophitlist[i]->FastToTotal());
      }
      std::cerr << "Total of:" << std::endl;
      std::cerr << " - Othe : " << map_of_ophit[kUnknown].size() << " opt hits" << std::endl;
      std::cerr << " - Sign : " << TotGen_Sign << " true parts\t| " << map_of_ophit[kSign].size() << " opt hits" << std::endl;
      std::cerr << " - Ar39 : " << TotGen_Ar39 << " true parts\t| " << map_of_ophit[kAr39].size() << " opt hits" << std::endl;
      std::cerr << " - Neut : " << TotGen_Neut << " true parts\t| " << map_of_ophit[kNeut].size() << " opt hits" << std::endl;
      std::cerr << " - Kryp : " << TotGen_Kryp << " true parts\t| " << map_of_ophit[kKryp].size() << " opt hits" << std::endl;
      std::cerr << " - Rdon : " << TotGen_Rdon << " true parts\t| " << map_of_ophit[kRdon].size() << " opt hits" << std::endl;
      std::cerr << " - Ar42 : " << TotGen_Ar42 << " true parts\t| " << map_of_ophit[kAr42].size() << " opt hits" << std::endl;
    }
    else {
      mf::LogWarning("DualPhaseOpHitDumperModule") << "Requested to save optical hits, but cannot load any ophits";
    }
  }

  


  std::cout << "True_Bck_Mode.size() " << True_Bck_Mode.size() << std::endl;
  fAnaTree->Fill();
}

void DualPhaseOpHitDumper::FillTruth(const art::FindManyP<simb::MCParticle> Assn,
                      const art::Handle<std::vector<simb::MCTruth>>& Hand,
                      const PType type) {

  for(size_t i1=0; i1<Hand->size(); ++i1) {
    for ( size_t i2=0; i2 < Assn.at(i1).size(); ++i2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(i1).at(i2));
      True_Bck_Mode  .push_back(type);
      True_Bck_VertX .push_back(ThisPar.Vx());
      True_Bck_VertY .push_back(ThisPar.Vy());
      True_Bck_VertZ .push_back(ThisPar.Vz());
      True_Bck_Time  .push_back(ThisPar.T());
      True_Bck_Energy.push_back(ThisPar.E());
      True_Bck_PDG   .push_back(ThisPar.PdgCode());
      True_Bck_ID    .push_back(ThisPar.TrackId());
    }
  }
}


//   auto SignTrue = evt.getHandle< std::vector<simb::MCTruth> >(fSIGNLabel);
//   if (SignTrue) {
//    art::FindManyP<simb::MCParticle> SignAssn(SignTrue,evt,fGEANTLabel);


void DualPhaseOpHitDumper::FillMyMaps(std::map< int, simb::MCParticle> &MyMap,
                       art::FindManyP<simb::MCParticle> Assn,
                       art::Handle< std::vector<simb::MCTruth> > Hand,
                       std::map<int, int>* indexMap) {
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    //std::cout << "Eve particle " << L1 << " " << Hand->at(L1).NParticles() << std::endl;
    //for(int i=0;i<Hand->at(L1).NParticles();i++) std::cout << Hand->at(L1).GetParticle(i).TrackId() << " " << Hand->at(L1).GetParticle(i).PdgCode() << std::endl;
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar; //std::cout << ThisPar.TrackId() << "trackID added to MyMap, PDG: " << ThisPar.PdgCode() <<std::endl;
      if(indexMap) indexMap->insert({ThisPar.TrackId(), L1});
    }
  }
  return;
}


PType DualPhaseOpHitDumper::WhichParType( int TrID )
{
  PType ThisPType = kUnknown;
  auto const& it=trkIDToPType.find(TrID);
  if(it!=trkIDToPType.end()){
    ThisPType=it->second;
  }
  return ThisPType;

}


bool DualPhaseOpHitDumper::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap )
{
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if ( ParIt != ParMap.end() ) {
    return true;
  } else
    return false;
}

//......................................................



void DualPhaseOpHitDumper::endJob()
{
  std::cout << firstCatch << " " << secondCatch << " " << thirdCatch << std::endl;
}

DEFINE_ART_MODULE(DualPhaseOpHitDumper)
