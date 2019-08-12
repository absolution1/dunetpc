////////////////////////////////////////////////////////////////////////
// Class:       PDSPmatch
// Plugin Type: analyzer (art v3_00_00)
// File:        PDSPmatch_module.cc
//
// Generated at Thu Feb 14 21:47:27 2019 by Bryan Joseph Ramson using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include <numeric>
#include <sstream>
#include <string.h>
#include <bitset>
#include <vector>

#include "larana/OpticalDetector/OpFlashAlg.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
//#include "art/Framework/Services/Optional/TFileService.h"


#include "dunetpc/dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"
//#include "dunetpc/dune/Geometry/ProtoDUNESPCRTSorter.h"

#include "dune-raw-data/Overlays/CRTFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace pdsp {
  class PDSPmatch;
}


class pdsp::PDSPmatch : public art::EDAnalyzer {
public:
  explicit PDSPmatch(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPmatch(PDSPmatch const&) = delete;
  PDSPmatch(PDSPmatch&&) = delete;
  PDSPmatch& operator=(PDSPmatch const&) = delete;
  PDSPmatch& operator=(PDSPmatch&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  const art::InputTag fTimeLabel; //Label of timing products
  const art::InputTag fCTBLabel; //Label of pdsctbdata products
  const art::InputTag fCRTLabel; //Label of crt products
  const art::InputTag fOpHitLabel; //Label of ophits after rawdecoding
  const art::InputTag fPandoLabel; //Label for Pandora Tracks
  const art::InputTag fPFParListLabel; // Label for PF Particle containers
  const art::InputTag fPMLabel; //Label for PM Tracks
  
  const int64_t fCRTCTBOffset; //CRT offset to global trigger in 20mus ticks
  const uint64_t fCRTWindow; //coincidence window for CRT triggers in 20mus ticks
  const bool fMC; //MC flag
  
  TTree *fTree;
  int run;
  int subrun;
  int event;
  int64_t timing_time;
  int64_t fCTB_time;
  uint32_t fCTBChan;
  
  //  TH1F *CTBCRTMap[32];
  //TH1F *CTBPairPE[16];
  //TH1F *CosTrigger[16];
  //TH1F *CRTPairPE[16];
  //TH2F *PandaTrackProUp[16];
  //TH2F *PandaTrackProDown[16];
  
  //TH1F *PMTrack[16];
  //TH1F *CRTPairPE[16];
  //TH1F *trig

  std::vector<uint32_t> fCRTChan;
  std::vector<int64_t> fPDS_time, fPando_time, fPM_time, fCRT_time;
  
  std::vector<double> fOpChan, fPE;//, ft0, fx0, fy0, fz0;
  
  std::vector<double> fTrkStartx_Pando, fTrkStarty_Pando, fTrkStartz_Pando, fTrkEndx_Pando, fTrkEndy_Pando, fTrkEndz_Pando;
  std::vector<double> fTrkTheta_Pando, fTrkPhi_Pando, fTrkUSPos1x_Pando, fTrkUSPos1y_Pando, fTrkUSPos2x_Pando, fTrkUSPos2y_Pando, fTrkDSPosx_Pando, fTrkDSPosy_Pando; 
  std::vector<double> fTrigCos_Pando;
  std::vector<uint32_t> fTrkProUS_Pando, fTrkProDS_Pando;
  std::vector<double> fTrkUSPixelProx_Pando, fTrkUSPixelProy_Pando; 
  std::vector<double> fTrkDSPixelProx_Pando, fTrkDSPixelProy_Pando; 
  
  std::vector<double> fTrkStartx_PM, fTrkStarty_PM, fTrkStartz_PM, fTrkEndx_PM, fTrkEndy_PM, fTrkEndz_PM;
  std::vector<double> fTrkTheta_PM, fTrkPhi_PM, fTrkUSPos1x_PM, fTrkUSPos1y_PM, fTrkUSPos2x_PM, fTrkUSPos2y_PM, fTrkDSPosx_PM, fTrkDSPosy_PM; 
  std::vector<double> fTrigCos_PM;
  std::vector<uint32_t> fTrkProUS_PM, fTrkProDS_PM;
  std::vector<double> fTrkUSPixelProx_PM, fTrkUSPixelProy_PM; 
  std::vector<double> fTrkDSPixelProx_PM, fTrkDSPixelProy_PM; 
 
  
  std::vector<TVector3> PandoTracks, PMTracks;
  // Declare member data here.
  
  TVector3 CTBPixelCtr[32];
};


pdsp::PDSPmatch::PDSPmatch(fhicl::ParameterSet const& p)
  : 
  EDAnalyzer(p), 
  fTimeLabel(p.get<art::InputTag>("TimingLabel")),
  fCTBLabel(p.get<art::InputTag>("CTBLabel")),
  fCRTLabel(p.get<art::InputTag>("CRTLabel")),
  fOpHitLabel(p.get<art::InputTag>("OpHitLabel")),
  fPandoLabel(p.get<art::InputTag>("PandoraLabel")),
  fPFParListLabel(p.get<art::InputTag>("PFParListLabel")),
  fPMLabel(p.get<art::InputTag>("PMLabel")),
  
  fCRTCTBOffset(p.get<int64_t>("CRTCTBOffset")),
  fCRTWindow(p.get<uint64_t>("CRTWindow")),
  fMC(p.get<bool>("MC"))
  //: EDAnalyzer{p}  // ,
  //More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  consumes<std::vector<raw::RDTimeStamp>>(fTimeLabel);
  consumes<std::vector<raw::ctb::pdspctb>>(fCTBLabel);
  consumes<std::vector<CRT::Trigger>>(fCRTLabel);
  consumes<std::vector<recob::OpHit>>(fOpHitLabel);
  consumes<std::vector<recob::Track>>(fPandoLabel);
  consumes<std::vector<recob::PFParticle>>(fPFParListLabel);
  consumes<std::vector<recob::Track>>(fPMLabel);
}

void pdsp::PDSPmatch::beginJob(){
  art::ServiceHandle<art::TFileService> tFileService;
  //for(int o=0;o<32;o++) CTBCRTMap[o] = tFileService->make<TH1F>(Form("CTB_%d",o),Form("CTB_%d",o),34,-1,32);
  /*for(int w=0;w<16;w++){
    CTBPairPE[w] = tFileService->make<TH1F>(Form("PE_%d",w),Form("PE_%d",w),288,-1,288);
    PandaTrackProUp[w] = tFileService->make<TH2F>(Form("PandaTrackUS_%d",w),Form("PandaTrackUS_%d",w),120,-600,600,120,-600,600);
    PandaTrackProDown[w] = tFileService->make<TH2F>(Form("PandaTrackDS_%d",w),Form("PandaTrackDS_%d",w),120,-600,600,120,-600,600);
    CosTrigger[w] = tFileService->make<TH1F>(Form("CosTrigger_%d",w),Form("CosTrigger_%d",w),25,0,1.0);
    }*/
  CTBPixelCtr[0].SetXYZ(-111.494,580.249,-972.204);
  CTBPixelCtr[1].SetXYZ(58.5057,580.249,-972.204);
  CTBPixelCtr[2].SetXYZ(307.506,580.249,-267.204);
  CTBPixelCtr[3].SetXYZ(477.506,580.249,-267.204);
  CTBPixelCtr[4].SetXYZ(477.506,410.249,-267.204);
  CTBPixelCtr[5].SetXYZ(477.506,235.249,-267.204);
  CTBPixelCtr[6].SetXYZ(477.506,65.2494,-267.204);
  CTBPixelCtr[7].SetXYZ(307.506,65.2494,-267.204);
  CTBPixelCtr[8].SetXYZ(58.5057,65.2494,-972.204);
  CTBPixelCtr[9].SetXYZ(-111.494,65.2494,-972.204);
  CTBPixelCtr[10].SetXYZ(-111.494,235.249,-972.204);
  CTBPixelCtr[11].SetXYZ(-111.494,410.249,-972.204);
  CTBPixelCtr[12].SetXYZ(58.5057,410.249,-972.204);
  CTBPixelCtr[13].SetXYZ(307.506,410.249,-267.204);
  CTBPixelCtr[14].SetXYZ(307.506,235.249,-267.204);
  CTBPixelCtr[15].SetXYZ(58.5057,235.249,-972.204);
  CTBPixelCtr[16].SetXYZ(293.506,580.249,1085.55);
  CTBPixelCtr[17].SetXYZ(123.506,580.249,1085.55);
  CTBPixelCtr[18].SetXYZ(-51.4993,580.249,1085.55);
  CTBPixelCtr[19].SetXYZ(-221.494,580.249,1085.55);
  CTBPixelCtr[20].SetXYZ(-221.494,410.249,1085.55);
  CTBPixelCtr[21].SetXYZ(-221.494,235.249,1085.55);
  CTBPixelCtr[22].SetXYZ(-221.494,65.2494,1085.55);
  CTBPixelCtr[23].SetXYZ(-51.4943,65.2494,1085.55);
  CTBPixelCtr[24].SetXYZ(123.506,65.2494,1085.55);
  CTBPixelCtr[25].SetXYZ(293.506,65.2494,1085.55);
  CTBPixelCtr[26].SetXYZ(293.506,235.249,1085.55);
  CTBPixelCtr[27].SetXYZ(293.506,410.249,1085.55);
  CTBPixelCtr[28].SetXYZ(123.506,410.249,1085.55);
  CTBPixelCtr[29].SetXYZ(-51.4993,410.249,1085.55);
  CTBPixelCtr[30].SetXYZ(-51.4943,235.249,1085.55);
  CTBPixelCtr[31].SetXYZ(123.506,235.249,1085.55);
  
  for(int h=0;h<32;h++){
    if(h==0 || h==1 || h==11 || h==12 || h ==10 || h==15 || h==9 || h==8) CTBPixelCtr[h].SetY(CTBPixelCtr[h].Y()-46.0); //US Jura/Beam Left
    if(h==2 || h==3 || h==13 || h==4 || h ==14 || h==5 || h==7 || h==6) CTBPixelCtr[h].SetY(CTBPixelCtr[h].Y()-47.5); //US Saleve/Beam Right
    if(h==16 || h==17 || h==27 || h==28 || h ==26 || h==31 || h==25 || h==24) CTBPixelCtr[h].SetY(CTBPixelCtr[h].Y()-151.0); //DS Saleve/Beam Right
    if(h==18 || h==19 || h==29 || h==20 || h ==21 || h==30 || h==22 || h==23) CTBPixelCtr[h].SetY(CTBPixelCtr[h].Y()-150.0); //DS Jura/Beam Left
  }								    
  art::ServiceHandle < geo::Geometry > geom;
  /*  for (unsigned int i=0; i<32; i++){
    const auto & trigGeo = geom -> AuxDet(i);
    //    float xpos=0.0,ypos=0.0,zpos=0.0;
    //if(i<4){
    const auto & hit1Geo = trigGeo.SensitiveVolume(31);
    const auto & hit2Geo = trigGeo.SensitiveVolume(32);
    //const auto & hit3Geo = trigGeo.SensitiveVolume(47);
    //const auto & hit4Geo = trigGeo.SensitiveVolume(48);
    const auto hit1Center = hit1Geo.GetCenter();
    const auto hit2Center = hit2Geo.GetCenter();
    //const auto hit3Center = hit3Geo.GetCenter();
    //const auto hit4Center = hit4Geo.GetCenter();
    
    //for (unsigned int j=0; j<64; j++){
    
    std::cout<<"Module "<<i<<" First Position X: "<<(hit1Center.X()+hit2Center.X())/2.0<<" First Position Y: "<<(hit1Center.Y()+hit2Center.Y())/2.0<<" First Position Z: "<<(hit1Center.Z()+hit2Center.Z())/2.0<<std::endl;
    //std::cout<<"Module "<<i<<" Second Position X: " <<(hit3Center.X()+hit4Center.X())/2.0<< " Second Position Y: "<<(hit3Center.Y()+hit4Center.Y())/2.0 << (hit3Center.Z()+hit4Center.Z())/2.0 << std::endl;
    }*/

  fTree = tFileService->make<TTree>("ProtoDUNE_Evt_Match","ProtoDUNE Matched Event Tree");
  
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("timing_time",&timing_time,"timing_time/L");
  
  fTree->Branch("fCRT_time",&fCRT_time);
  fTree->Branch("fCRTChan", &fCRTChan);
  fTree->Branch("fCTB_time",&fCTB_time,"fCTB_time/L");
  fTree->Branch("fCTBChan", &fCTBChan,"fCTBChan/i");
  
  fTree->Branch("fPDS_time",&fPDS_time);
  fTree->Branch("fOpChan", &fOpChan);
  fTree->Branch("fPE",&fPE);
 
  fTree->Branch("fPando_time",&fPando_time);
  fTree->Branch("fTrkStartx_Pando",&fTrkStartx_Pando);
  fTree->Branch("fTrkStarty_Pando",&fTrkStarty_Pando);
  fTree->Branch("fTrkStartz_Pando",&fTrkStartz_Pando);
  fTree->Branch("fTrkEndx_Pando",&fTrkEndx_Pando);
  fTree->Branch("fTrkEndy_Pando",&fTrkEndy_Pando);
  fTree->Branch("fTrkEndz_Pando",&fTrkEndz_Pando);
  /*fTree->Branch("fTrkProUS_Pando",&fTrkProUS_Pando);
  fTree->Branch("fTrkProDS_Pando",&fTrkProDS_Pando);
  fTree->Branch("fTrkTheta_Pando",&fTrkTheta_Pando);
  fTree->Branch("fTrkPhi_Pando",&fTrkPhi_Pando);
  fTree->Branch("fTrkUSPos1x_Pando",&fTrkUSPos1x_Pando);
  fTree->Branch("fTrkUSPos1y_Pando",&fTrkUSPos1y_Pando);
  fTree->Branch("fTrkUSPos2x_Pando",&fTrkUSPos2x_Pando);
  fTree->Branch("fTrkUSPos2y_Pando",&fTrkUSPos2y_Pando);
  fTree->Branch("fTrkDSPosx_Pando",&fTrkDSPosx_Pando);
  fTree->Branch("fTrkDSPosy_Pando",&fTrkDSPosy_Pando);
  fTree->Branch("fTrigCos_Pando",&fTrigCos_Pando);
  fTree->Branch("fTrkUSPixelProx_Pando",&fTrkUSPixelProx_Pando);
  fTree->Branch("fTrkUSPixelProy_Pando",&fTrkUSPixelProy_Pando);
  fTree->Branch("fTrkDSPixelProx_Pando",&fTrkDSPixelProx_Pando);
  fTree->Branch("fTrkDSPixelProy_Pando",&fTrkDSPixelProy_Pando);*/
  
  fTree->Branch("fPM_time",&fPM_time);
  fTree->Branch("fTrkStartx_PM",&fTrkStartx_PM);
  fTree->Branch("fTrkStarty_PM",&fTrkStarty_PM);
  fTree->Branch("fTrkStartz_PM",&fTrkStartz_PM);
  fTree->Branch("fTrkEndx_PM",&fTrkEndx_PM);
  fTree->Branch("fTrkEndy_PM",&fTrkEndy_PM);
  fTree->Branch("fTrkEndz_PM",&fTrkEndz_PM);  
  /*fTree->Branch("fTrkProUS_PM",&fTrkProUS_PM);
  fTree->Branch("fTrkProDS_PM",&fTrkProDS_PM);
  fTree->Branch("fTrkTheta_PM",&fTrkTheta_PM);
  fTree->Branch("fTrkPhi_PM",&fTrkPhi_PM);
  fTree->Branch("fTrkUSPos1x_PM",&fTrkUSPos1x_PM);
  fTree->Branch("fTrkUSPos1y_PM",&fTrkUSPos1y_PM);
  fTree->Branch("fTrkUSPos2x_PM",&fTrkUSPos2x_PM);
  fTree->Branch("fTrkUSPos2y_PM",&fTrkUSPos2y_PM);
  fTree->Branch("fTrkDSPosx_PM",&fTrkDSPosx_PM);
  fTree->Branch("fTrkDSPosy_PM",&fTrkDSPosy_PM);
  fTree->Branch("fTrigCos_PM",&fTrigCos_PM);
  fTree->Branch("fTrkUSPixelProx_PM",&fTrkUSPixelProx_PM);
  fTree->Branch("fTrkUSPixelProy_PM",&fTrkUSPixelProy_PM);
  fTree->Branch("fTrkDSPixelProx_PM",&fTrkDSPixelProx_PM);
  fTree->Branch("fTrkDSPixelProy_PM",&fTrkDSPixelProy_PM);*/

}

void pdsp::PDSPmatch::analyze(art::Event const& e){
  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  
  fCRT_time.clear();
  fCRTChan.clear();
  
  //  fCTB_time.clear();
  //fCTBChan.clear();
  
  fPDS_time.clear();
  fOpChan.clear();
  fPE.clear();
  
  fPando_time.clear();
  fTrkStartx_Pando.clear();
  fTrkStarty_Pando.clear();
  fTrkStartz_Pando.clear();
  fTrkEndx_Pando.clear();
  fTrkEndy_Pando.clear();
  fTrkEndz_Pando.clear();
  //fTrkProUS_Pando.clear();
  //fTrkProDS_Pando.clear();
  fTrkTheta_Pando.clear();
  fTrkPhi_Pando.clear();
  //fTrkUSPos1x_Pando.clear();
  //fTrkUSPos1y_Pando.clear();
  //fTrkUSPos2x_Pando.clear();
  //fTrkUSPos2y_Pando.clear();
  //fTrkDSPosx_Pando.clear();
  //fTrkDSPosy_Pando.clear();
  //fTrigCos_Pando.clear();
  //fTrkUSPixelProx_Pando.clear();
  //fTrkUSPixelProy_Pando.clear();
  //fTrkDSPixelProx_Pando.clear();
  //fTrkDSPixelProy_Pando.clear();

  
  fPM_time.clear();
  fTrkStartx_PM.clear();
  fTrkStarty_PM.clear();
  fTrkStartz_PM.clear();
  fTrkEndx_PM.clear();
  fTrkEndy_PM.clear();
  fTrkEndz_PM.clear();
  //fTrkProUS_PM.clear();
  //fTrkProDS_PM.clear();
  fTrkTheta_PM.clear();
  fTrkPhi_PM.clear();
  //fTrkUSPos1x_PM.clear();
  //fTrkUSPos1y_PM.clear();
  //fTrkUSPos2x_PM.clear();
  //fTrkUSPos2y_PM.clear();
  //fTrkDSPosx_PM.clear();
  //fTrkDSPosy_PM.clear();
  //fTrigCos_PM.clear();
  //fTrkUSPixelProx_PM.clear();
  //fTrkUSPixelProy_PM.clear();
  //fTrkDSPixelProx_PM.clear();
  //fTrkDSPixelProy_PM.clear();

  int trigger = -1;
  const auto timeStamps = e.getValidHandle<std::vector<raw::RDTimeStamp>>(fTimeLabel);
  if(timeStamps.isValid() && timeStamps->size() == 1){
    const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
    trigger = timeStamp.GetFlags();
    timing_time = timeStamp.GetTimeStamp();
    //std::cout << "Timing: " << std::setw(20) << timing_time << std::endl;
  }
  else {
    mf::LogWarning("Empty Event TimeStamp") << "Invalid Event TimeStamp. Skipping. \n";
    return;
  }
    
  if (trigger !=13) return;
  
  const auto crtHandle = e.getValidHandle<std::vector<CRT::Trigger>>(fCRTLabel);
  if(crtHandle->empty()){
    mf::LogWarning("Empty CRT Fragment") << "Empty CRT fragments for this event. Skipping. \n";
    return;
  }
    
  const auto OpHitHandle = e.getValidHandle<std::vector<recob::OpHit>>(fOpHitLabel);
  if(OpHitHandle->empty()){
    mf::LogWarning("Empty OpHit Object") << "Empty OpHit Object. Error in retrieval. Skipping. \n";
    return;
  }
    
  const auto ctbHandle = e.getValidHandle<std::vector<raw::ctb::pdspctb>>(fCTBLabel);
  if(ctbHandle->empty()){
    mf::LogWarning("Empty CTB Fragment") << "Empty CTB fragment for this event. Skipping. \n";
    return;
  }
  
  if(ctbHandle->size() > 1){
    mf::LogWarning("Multiple CTB Triggers") << "Found " << ctbHandle->size() << " CTB data products. Skipping. \n";
    return;
  }
  
  const auto& ctbStatus = ctbHandle->front();
  const auto statuses = ctbStatus.GetChStatusAfterHLTs();
  
  //bool ctbus = false;
  //bool ctbds = false;
  //int crtus = 0;
  //int crtds = 0;
  bool pmt = false;
  bool pat = false;
  //bool ctbtrig[32] = { false };
  //bool crttrig[32] = { false };

  for(const auto& status: statuses){
    //ctbus = false;
    //ctbds = false;
    //crtus = 0;
    //crtds = 0;
    pmt = false;
    pat = false;
    //std::memset(ctbtrig,false,32);
    //std::memset(crttrig,false,32);

    
    fCRT_time.clear();
    fCRTChan.clear();
    
    //fCTB_time.clear();
    //fCTBChan.clear();
    
    fPDS_time.clear();
    fOpChan.clear();
    fPE.clear();
    
    fPando_time.clear();
    fTrkStartx_Pando.clear();
    fTrkStarty_Pando.clear();
    fTrkStartz_Pando.clear();
    fTrkEndx_Pando.clear();
    fTrkEndy_Pando.clear();
    fTrkEndz_Pando.clear();
    //fTrkProUS_Pando.clear();
    //fTrkProDS_Pando.clear();
    //fTrkTheta_Pando.clear();
    //fTrkPhi_Pando.clear();
    //fTrkUSPos1x_Pando.clear();
    //fTrkUSPos1y_Pando.clear();
    //fTrkUSPos2x_Pando.clear();
    //fTrkUSPos2y_Pando.clear();
    //fTrkDSPosx_Pando.clear();
    //fTrkDSPosy_Pando.clear();
    //fTrigCos_Pando.clear();
    //fTrkUSPixelProx_Pando.clear();
    //fTrkUSPixelProy_Pando.clear();
    //fTrkDSPixelProx_Pando.clear();
    //fTrkDSPixelProy_Pando.clear();


    fPM_time.clear();
    fTrkStartx_PM.clear();
    fTrkStarty_PM.clear();
    fTrkStartz_PM.clear();
    fTrkEndx_PM.clear();
    fTrkEndy_PM.clear();
    fTrkEndz_PM.clear();
    //fTrkProUS_PM.clear();
    //fTrkProDS_PM.clear();
    //fTrkTheta_PM.clear();
    //fTrkPhi_PM.clear();
    //fTrkUSPos1x_PM.clear();
    //fTrkUSPos1y_PM.clear();
    //fTrkUSPos2x_PM.clear();
    //fTrkUSPos2y_PM.clear();
    //fTrkDSPosx_PM.clear();
    //fTrkDSPosy_PM.clear();
    //fTrigCos_PM.clear();
    //fTrkUSPixelProx_PM.clear();
    //fTrkUSPixelProy_PM.clear();
    //fTrkDSPixelProx_PM.clear();
    //fTrkDSPixelProy_PM.clear();

   
    if(status.crt != 0) {
      fCTB_time = status.timestamp;
      fCTBChan = status.crt;
      if(!crtHandle->empty()){ 
	std::vector<CRT::Trigger> inWindow(crtHandle->size());
	const int64_t ctbetime = status.timestamp;
	const auto newEnd = std::copy_if(crtHandle->begin(), crtHandle->end(), inWindow.begin(),[this,ctbetime](const auto& trigger){
	    return (uint64_t)(labs(ctbetime - ((int64_t)(trigger.Timestamp()) - fCRTCTBOffset))) < fCRTWindow; });
	inWindow.resize(std::distance(inWindow.begin(), newEnd));   
	for(size_t k=0;k<inWindow.size();++k){
	  //std::cout << "CRT: " << inWindow[k].Timestamp() << " fCRTChan: " << inWindow[k].Channel() << std::endl;
	  fCRT_time.push_back(inWindow[k].Timestamp());
	  fCRTChan.push_back(inWindow[k].Channel());
	  //if(inWindow[k].Channel()<16) crtus++;
	  //else crtds++;
	}
      }
      //const std::string binary = std::bitset<32>(status.crt).to_string();
      //std::cout << "CTB: " << status.timestamp << " fCTBChan: " << status.crt << " " << binary.length() <<  std::endl;
      //std::cout << "CTB: " << status.timestamp << " fCTBChan: " << binary << std::endl;
	 
      //for(std::string::size_type l = binary.length();static_cast<int>(l)>-1;--l) {
      //if(binary[l]=='1') {
	  //std::cout << "CTB: " << status.timestamp << " fCTBChan: " << (binary.length()-1)-l << std::endl;
     	  //fCTBChan.push_back((binary.length()-1)-l);
	  //if((binary.length()-1)-l < 16) ctbus = true;
	  //else ctbds = true;
	    //}
	    //}
      //      if(us == false || ds == false) continue;
     
      
      for(const auto& OpHit: *OpHitHandle){
	// if((timing_time-OpHit.PeakTime == 112){
	if(OpHit.PE() < 10.0 || OpHit.PE() > 1000.00) continue;
	//std::cout << "SSP Peaktime: " << static_cast<int64_t>(OpHit.PeakTime())/3 << " SSP Peaktime Abs: " << static_cast<int64_t>(OpHit.PeakTimeAbs())/3 << " OpChannel: " << OpHit.OpChannel() << " PE: " << OpHit.PE() <<  std::endl;
	fPDS_time.push_back((OpHit.PeakTime()/3));
	fOpChan.push_back(OpHit.OpChannel());
	fPE.push_back(OpHit.PE());
	//}
      }
      
      art::Handle<std::vector<recob::Track>> PandoTrkHandle;
      std::vector<art::Ptr<recob::Track>> PandoTrk;
      if(e.getByLabel(fPandoLabel,PandoTrkHandle)) art::fill_ptr_vector(PandoTrk,PandoTrkHandle);
    
      art::Handle<std::vector<recob::PFParticle>> PFParListHandle;
      e.getByLabel(fPFParListLabel,PFParListHandle);
      art::FindManyP<recob::PFParticle> PFPar(PandoTrkHandle,e,fPandoLabel);
      art::FindManyP<anab::T0> PFT0(PFParListHandle,e,fPFParListLabel);
    
   
      for(size_t p = 0;p<PandoTrk.size();++p){
	auto & Trk = PandoTrk[p];
	if(!((Trk->Vertex().Z() < 40) || (Trk->End().Z() < 40))) continue;
	if(!((Trk->Vertex().Z() > 660) || (Trk->End().Z() > 660))) continue;
	fTrkStartx_Pando.push_back(Trk->Vertex().X());
	fTrkStarty_Pando.push_back(Trk->Vertex().Y());
	fTrkStartz_Pando.push_back(Trk->Vertex().Z());
	fTrkEndx_Pando.push_back(Trk->End().X());
	fTrkEndy_Pando.push_back(Trk->End().Y());
	fTrkEndz_Pando.push_back(Trk->End().Z());
	double t0temp = 0;
	auto &PFPS = PFPar.at(Trk.key());
	if(!PFPS.empty()){
	  auto &T0S = PFT0.at(PFPS[0].key());
	  if(!T0S.empty()){
	    t0temp = T0S[0]->Time();
	  }
	}
	fPando_time.push_back(t0temp);
	//TVector3 ttrack;
	//ttrack.SetXYZ((Trk->End().X()-Trk->Vertex().X()),(Trk->End().Y()-Trk->Vertex().Y()),(Trk->End().Z()-Trk->Vertex().Z()));
	//PandoTracks.push_back(ttrack);
	//float paramus1 = (-972.204-(fTrkStartz_Pando[p]))/PandoTracks[p].Z();
	//float paramus2 = (-267.204-(fTrkStartz_Pando[p]))/PandoTracks[p].Z();
	//float paramds = (1085.55-(fTrkStartz_Pando[p]))/PandoTracks[p].Z();
	//uint32_t TrkUS=99;
	//uint32_t TrkDS=99;
	//float ftrig = -9.9;
	//float prousx = -999.0;
	//float prousy = -999.0;
	//float prodsx = -999.0;
	//float prodsy = -999.0;

	//fTrkTheta_Pando.push_back(PandoTracks[p].Theta());
	//fTrkPhi_Pando.push_back(PandoTracks[p].Phi());
	//fTrkUSPos1x_Pando.push_back(PandoTracks[p].X()*paramus1+fTrkStartx_Pando[p]);
	//fTrkUSPos1y_Pando.push_back(PandoTracks[p].Y()*paramus1+fTrkStarty_Pando[p]);
	//fTrkUSPos2x_Pando.push_back(PandoTracks[p].X()*paramus2+fTrkStartx_Pando[p]);
	//fTrkUSPos2y_Pando.push_back(PandoTracks[p].Y()*paramus2+fTrkStarty_Pando[p]);
	//fTrkDSPosx_Pando.push_back(PandoTracks[p].X()*paramds+fTrkStartx_Pando[p]);
	//fTrkDSPosy_Pando.push_back(PandoTracks[p].Y()*paramds+fTrkStarty_Pando[p]);
	/*for(size_t v=0;v<32;v++){
	  if(((CTBPixelCtr[v].X()+80.0)>fTrkUSPos1x_Pando[p]) && (fTrkUSPos1x_Pando[p]>(CTBPixelCtr[v].X()-80.0)) && ((CTBPixelCtr[v].Y()+80.0)>fTrkUSPos1y_Pando[p]) && (fTrkUSPos1y_Pando[p]>(CTBPixelCtr[v].Y()-80.0))){
	    TrkUS=v;
	    prousx = fTrkUSPos1x_Pando[p]-CTBPixelCtr[v].X();
	    prousy = fTrkUSPos1y_Pando[p]-CTBPixelCtr[v].Y();
	  }
	  if(((CTBPixelCtr[v].X()+80.0)>fTrkUSPos2x_Pando[p]) && (fTrkUSPos2x_Pando[p]>(CTBPixelCtr[v].X()-80.0)) && ((CTBPixelCtr[v].Y()+80.0)>fTrkUSPos2y_Pando[p]) && (fTrkUSPos2y_Pando[p]>(CTBPixelCtr[v].Y()-80.0))){
	    TrkUS=v;
	    prousx = fTrkUSPos2x_Pando[p]-CTBPixelCtr[v].X();
	    prousy = fTrkUSPos2y_Pando[p]-CTBPixelCtr[v].Y();
	  }
	  if(((CTBPixelCtr[v].X()+80.0)>fTrkDSPosx_Pando[p]) && (fTrkDSPosx_Pando[p]>(CTBPixelCtr[v].X()-80.0)) && ((CTBPixelCtr[v].Y()+80.0)>fTrkDSPosy_Pando[p]) && (fTrkDSPosy_Pando[p]>(CTBPixelCtr[v].Y()-80.0))) {
	    TrkDS=v;
	    prodsx = fTrkDSPosx_Pando[p]-CTBPixelCtr[v].X();
	    prodsy = fTrkDSPosy_Pando[p]-CTBPixelCtr[v].Y();
	  }  
	}
	fTrkProUS_Pando.push_back(TrkUS);
	fTrkProDS_Pando.push_back(TrkDS);
        fTrkUSPixelProx_Pando.push_back(prousx);
        fTrkUSPixelProy_Pando.push_back(prousy);
	fTrkDSPixelProx_Pando.push_back(prodsx);
        fTrkDSPixelProy_Pando.push_back(prodsy);
	
	if(TrkUS < 99 && TrkDS < 99) {
	  TVector3 CompVec = CTBPixelCtr[TrkDS]-CTBPixelCtr[TrkUS];
	  ftrig = CompVec.Unit()*ttrack.Unit();
	}
	fTrigCos_Pando.push_back(ftrig);
	//std::cout << "Pandora: T0: " << t0temp << " X0: " <<  Trk->Vertex().X() << " Y0: " <<  Trk->Vertex().Y() << " Z0: " <<  Trk->Vertex().Z() << " XF: " <<  Trk->End().X() << " YF: " << Trk->End().Y() << " ZF: " <<  Trk->End().Z() << std::endl;
	*/
      }
       if(PandoTrk.size() > 0) pat = true;

      art::Handle<std::vector<recob::Track>> PMTrkHandle;
      std::vector<art::Ptr<recob::Track>> PMTrks;
      if(e.getByLabel(fPMLabel,PMTrkHandle)) art::fill_ptr_vector(PMTrks,PMTrkHandle);
      
      art::FindManyP<anab::T0> PPMT0(PMTrkHandle,e,fPMLabel);
      
      for(size_t r = 0;r<PMTrks.size();++r){
	auto & Trk = PMTrks[r];
	if(!((Trk->Vertex().Z() < 40) || (Trk->End().Z() < 40))) continue;
	if(!((Trk->Vertex().Z() > 660) || (Trk->End().Z() > 660))) continue;
	fTrkStartx_PM.push_back(Trk->Vertex().X());
	fTrkStarty_PM.push_back(Trk->Vertex().Y());
	fTrkStartz_PM.push_back(Trk->Vertex().Z());
	fTrkEndx_PM.push_back(Trk->End().X());
	fTrkEndy_PM.push_back(Trk->End().Y());
	fTrkEndz_PM.push_back(Trk->End().Z());
	double t0temp = 0;
	auto &T0S = PPMT0.at(Trk.key());
	if(!T0S.empty()){
	  t0temp = T0S[0]->Time();
	}
	fPM_time.push_back(t0temp);
	//std::cout << "PM: T0: " << t0temp << " X0: " <<  Trk->Vertex().X() << " Y0: " <<  Trk->Vertex().Y() << " Z0: " <<  Trk->Vertex().Z() << " XF: " <<  Trk->End().X() << " YF: " << Trk->End().Y() << " ZF: " <<  Trk->End().Z() << std::endl;
	/*TVector3 ttrack;
	ttrack.SetXYZ((Trk->End().X()-Trk->Vertex().X()),(Trk->End().Y()-Trk->Vertex().Y()),(Trk->End().Z()-Trk->Vertex().Z()));
	PMTracks.push_back(ttrack);
	float paramus1 = (-972.204-(fTrkStartz_PM[r]))/PMTracks[r].Z();
	float paramus2 = (-267.204-(fTrkStartz_PM[r]))/PMTracks[r].Z();
	float paramds = (1085.55-(fTrkStartz_PM[r]))/PMTracks[r].Z();
	uint32_t TrkUS=99;
	uint32_t TrkDS=99;
	float ftrig = -9.9;
	float prousx = -999.0;
	float prousy = -999.0;
	float prodsx = -999.0;
	float prodsy = -999.0;

	
	fTrkTheta_PM.push_back(PMTracks[r].Theta());
	fTrkPhi_PM.push_back(PMTracks[r].Phi());
	fTrkUSPos1x_PM.push_back(PMTracks[r].X()*paramus1+fTrkStartx_PM[r]);
	fTrkUSPos1y_PM.push_back(PMTracks[r].Y()*paramus1+fTrkStarty_PM[r]);
	fTrkUSPos2x_PM.push_back(PMTracks[r].X()*paramus2+fTrkStartx_PM[r]);
	fTrkUSPos2y_PM.push_back(PMTracks[r].Y()*paramus2+fTrkStarty_PM[r]);
	fTrkDSPosx_PM.push_back(PMTracks[r].X()*paramds+fTrkStartx_PM[r]);
	fTrkDSPosy_PM.push_back(PMTracks[r].Y()*paramds+fTrkStarty_PM[r]);
	for(size_t v=0;v<32;v++){
	  if(((CTBPixelCtr[v].X()+80.0)>fTrkUSPos1x_PM[r]) && (fTrkUSPos1x_PM[r]>(CTBPixelCtr[v].X()-80.0)) && ((CTBPixelCtr[v].Y()+80.0)>fTrkUSPos1y_PM[r]) && (fTrkUSPos1y_PM[r]>(CTBPixelCtr[v].Y()-80.0))){
	    TrkUS=v;
	    prousx = fTrkUSPos1x_PM[r]-CTBPixelCtr[v].X();
	    prousy = fTrkUSPos1y_PM[r]-CTBPixelCtr[v].Y();
	  }
	  if(((CTBPixelCtr[v].X()+80.0)>fTrkUSPos2x_PM[r]) && (fTrkUSPos2x_PM[r]>(CTBPixelCtr[v].X()-80.0)) && ((CTBPixelCtr[v].Y()+80.0)>fTrkUSPos2y_PM[r]) && (fTrkUSPos2y_PM[r]>(CTBPixelCtr[v].Y()-80.0))){
	    TrkUS=v;
	    prousx = fTrkUSPos2x_PM[r]-CTBPixelCtr[v].X();
	    prousy = fTrkUSPos2y_PM[r]-CTBPixelCtr[v].Y();
	  }
	  if(((CTBPixelCtr[v].X()+80.0)>fTrkDSPosx_PM[r]) && (fTrkDSPosx_PM[r]>(CTBPixelCtr[v].X()-80.0)) && ((CTBPixelCtr[v].Y()+80.0)>fTrkDSPosy_PM[r]) && (fTrkDSPosy_PM[r]>(CTBPixelCtr[v].Y()-80.0))){
	    TrkDS=v;  
	    prousx = fTrkDSPosx_PM[r]-CTBPixelCtr[v].X();
	    prousy = fTrkDSPosy_PM[r]-CTBPixelCtr[v].Y();
	  }
	}
	fTrkProUS_PM.push_back(TrkUS);
	fTrkProDS_PM.push_back(TrkDS);
	fTrkUSPixelProx_PM.push_back(prousx);
        fTrkUSPixelProy_PM.push_back(prousy);
	fTrkDSPixelProx_PM.push_back(prodsx);
        fTrkDSPixelProy_PM.push_back(prodsy);

	if(TrkUS < 99 && TrkDS < 99) {
	  TVector3 CompVec = CTBPixelCtr[TrkDS]-CTBPixelCtr[TrkUS];
	  ftrig = CompVec.Unit()*ttrack.Unit();
	}
	fTrigCos_PM.push_back(ftrig);
	*/
      }
      if(PMTrks.size()>0) pmt = true;
    
  
  
    //fCRTChan.push_back(inWindow[k].Channel());
    //fCTBChan.push_back(status.crt);
      /* if(((ctbus && ctbds) || ((crtus > 2) && (crtds > 2)))){  
	if(pmt || pat){
	  for(size_t c=0;c<fCTBChan.size();c++) ctbtrig[fCTBChan[c]] = true;
	  for(size_t d=0;d<fCRTChan.size();d++) crttrig[fCRTChan[d]] = true;

	  if((ctbtrig[9] && ctbtrig[22]) || (crttrig[6] && crttrig[5] && crttrig[22] && crttrig[21])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[0]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[0]->Fill(fTrkUSPos1x_Pando[f],fTrkUSPos1y_Pando[f]);
	      PandaTrackProDown[0]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[22]-CTBPixelCtr[9];
	      CosTrigger[0]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[0]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[8] && ctbtrig[23]) || (crttrig[7] && crttrig[5] && crttrig[23] && crttrig[21])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[1]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[1]->Fill(fTrkUSPos1x_Pando[f],fTrkUSPos1y_Pando[f]);
	      PandaTrackProDown[1]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[23]-CTBPixelCtr[8];
	      CosTrigger[1]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[1]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[7] && ctbtrig[24]) || (crttrig[8] && crttrig[10] && crttrig[24] && crttrig[26])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[2]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[2]->Fill(fTrkUSPos2x_Pando[f],fTrkUSPos2y_Pando[f]);
	      PandaTrackProDown[2]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[24]-CTBPixelCtr[7];
	      CosTrigger[2]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[2]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[6] && ctbtrig[25]) || (crttrig[9] && crttrig[10] && crttrig[25] && crttrig[26])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[3]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[3]->Fill(fTrkUSPos2x_Pando[f],fTrkUSPos2y_Pando[f]);
	      PandaTrackProDown[3]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[25]-CTBPixelCtr[6];
	      CosTrigger[3]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[3]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[10] && ctbtrig[21]) || (crttrig[6] && crttrig[4] && crttrig[25] && crttrig[27])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[4]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[4]->Fill(fTrkUSPos1x_Pando[f],fTrkUSPos1y_Pando[f]);
	      PandaTrackProDown[4]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[21]-CTBPixelCtr[10];
	      CosTrigger[4]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[4]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[15] && ctbtrig[30]) || (crttrig[7] && crttrig[4] && crttrig[24] && crttrig[27])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[5]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[5]->Fill(fTrkUSPos1x_Pando[f],fTrkUSPos1y_Pando[f]);
	      PandaTrackProDown[5]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[30]-CTBPixelCtr[15];
	      CosTrigger[5]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[5]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[14] && ctbtrig[31]) || (crttrig[8] && crttrig[11] && crttrig[20] && crttrig[23])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[6]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[6]->Fill(fTrkUSPos2x_Pando[f],fTrkUSPos2y_Pando[f]);
	      PandaTrackProDown[6]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[31]-CTBPixelCtr[14];
	      CosTrigger[6]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[6]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[5] && ctbtrig[26]) || (crttrig[9] && crttrig[11] && crttrig[20] && crttrig[22])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[7]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[7]->Fill(fTrkUSPos2x_Pando[f],fTrkUSPos2y_Pando[f]);
	      PandaTrackProDown[7]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[26]-CTBPixelCtr[5];
	      CosTrigger[7]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[7]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[11] && ctbtrig[20]) || (crttrig[1] && crttrig[3] && crttrig[17] && crttrig[19])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[8]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[8]->Fill(fTrkUSPos1x_Pando[f],fTrkUSPos1y_Pando[f]);
	      PandaTrackProDown[8]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[20]-CTBPixelCtr[11];
	      CosTrigger[8]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[8]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[12] && ctbtrig[29]) || (crttrig[0] && crttrig[5] && crttrig[28] && crttrig[31])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[9]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[9]->Fill(fTrkUSPos1x_Pando[f],fTrkUSPos1y_Pando[f]);
	      PandaTrackProDown[9]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[29]-CTBPixelCtr[12];
	      CosTrigger[9]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[9]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[13] && ctbtrig[28]) || (crttrig[12] && crttrig[15] && crttrig[16] && crttrig[19])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[10]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[10]->Fill(fTrkUSPos2x_Pando[f],fTrkUSPos2y_Pando[f]);
	      PandaTrackProDown[10]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[28]-CTBPixelCtr[13];
	      CosTrigger[10]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[10]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[4] && ctbtrig[27]) || (crttrig[12] && crttrig[14] && crttrig[17] && crttrig[19])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[11]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[11]->Fill(fTrkUSPos2x_Pando[f],fTrkUSPos2y_Pando[f]);
	      PandaTrackProDown[11]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[27]-CTBPixelCtr[4];
	      CosTrigger[11]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[11]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[0] && ctbtrig[19]) || (crttrig[1] && crttrig[2] && crttrig[29] && crttrig[30])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[12]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[12]->Fill(fTrkUSPos1x_Pando[f],fTrkUSPos1y_Pando[f]);
	      PandaTrackProDown[12]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[19]-CTBPixelCtr[0];
	      CosTrigger[12]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[12]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[1] && ctbtrig[18]) || (crttrig[0] && crttrig[2] && crttrig[29] && crttrig[31])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[13]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[13]->Fill(fTrkUSPos1x_Pando[f],fTrkUSPos1y_Pando[f]);
	      PandaTrackProDown[13]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[18]-CTBPixelCtr[1];
	      CosTrigger[13]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[13]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[2] && ctbtrig[17]) || (crttrig[13] && crttrig[14] && crttrig[16] && crttrig[18])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[14]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[14]->Fill(fTrkUSPos2x_Pando[f],fTrkUSPos2y_Pando[f]);
	      PandaTrackProDown[14]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[17]-CTBPixelCtr[2];
	      CosTrigger[14]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[14]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	  }
	  if((ctbtrig[3] && ctbtrig[16]) || (crttrig[13] && crttrig[14] && crttrig[17] && crttrig[18])){
	    for(size_t e=0;e<fOpChan.size();e++) CTBPairPE[15]->Fill(fOpChan[e],fPE[e]);
	    for(size_t f=0;f<fPando_time.size();f++) {
	      PandaTrackProUp[15]->Fill(fTrkUSPos2x_Pando[f],fTrkUSPos2y_Pando[f]);
	      PandaTrackProDown[15]->Fill(fTrkDSPosx_Pando[f],fTrkDSPosy_Pando[f]);
	      TVector3 tempvec = CTBPixelCtr[16]-CTBPixelCtr[3];
	      CosTrigger[15]->Fill(tempvec.Unit()*PandoTracks[f].Unit());
	      //CosTrigger[15]->Fill(acos((tempvec*PandoTracks[f]) / (tempvec.Mag()*PandoTracks[f].Mag())));
	    }
	    }*/
	  
	  
      //for(size_t o=0;o<fCTBChan.size();o++) for(size_t p=0;p<fCRTChan.size();p++) CTBCRTMap[fCTBChan[o]]->Fill(fCRTChan[p]);
      if( pat||pmt) fTree->Fill();
	  //}
    }
  }
}
//}
  //ctb_time = status.timestamp;

  //Get CTB Stuff
  //const int64_t ctbtime = status.timestamp;
 


DEFINE_ART_MODULE(pdsp::PDSPmatch)
