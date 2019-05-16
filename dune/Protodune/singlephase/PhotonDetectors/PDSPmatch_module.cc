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

#include "dunetpc/dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"
//#include "dunetpc/dune/Geometry/ProtoDUNESPCRTSorter.h"

#include "dune-raw-data/Overlays/CRTFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "TTree.h"
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
  const int64_t fCRTCTBOffset;
 
  const uint64_t fCRTWindow;
  
  TTree *fTree;
  int run;
  int subrun;
  int event;
  int64_t timing_time;
  int64_t fCTB_time;
  
  std::vector<uint32_t> fCRTChan, fCTBChan;
  std::vector<int64_t> fCRT_time, fPDS_time, fPando_time, fPM_time;
  
  std::vector<double> fOpChan, fPE;//, ft0, fx0, fy0, fz0;
  
  std::vector<double> fTrkStartx_Pando, fTrkStarty_Pando, fTrkStartz_Pando, fTrkEndx_Pando, fTrkEndy_Pando, fTrkEndz_Pando;
  std::vector<double> fTrkStartx_PM, fTrkStarty_PM, fTrkStartz_PM, fTrkEndx_PM, fTrkEndy_PM, fTrkEndz_PM;
  
  // Declare member data here.

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
  fCRTWindow(p.get<uint64_t>("CRTWindow"))
  //: EDAnalyzer{p}  // ,
  // More initializers here.
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
  art::ServiceHandle < geo::Geometry > geom;
  for (unsigned int i=0; i<32; i++){
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
      //}
  }
  fTree = tFileService->make<TTree>("ProtoDUNE_Evt_Match","ProtoDUNE Matched Event Tree");
  
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("timing_time",&timing_time,"timing_time/L");
  
  fTree->Branch("fCRT_time",&fCRT_time);
  fTree->Branch("fCRTChan", &fCRTChan);
  
  fTree->Branch("fCTB_time",&fCTB_time,"fCTB_time/L");
  fTree->Branch("fCTBChan", &fCTBChan);
  
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
  
  fTree->Branch("fPM_time",&fPM_time);
  fTree->Branch("fTrkStartx_PM",&fTrkStartx_PM);
  fTree->Branch("fTrkStarty_PM",&fTrkStarty_PM);
  fTree->Branch("fTrkStartz_PM",&fTrkStartz_PM);
  fTree->Branch("fTrkEndx_PM",&fTrkEndx_PM);
  fTree->Branch("fTrkEndy_PM",&fTrkEndy_PM);
  fTree->Branch("fTrkEndz_PM",&fTrkEndz_PM);  
}

void pdsp::PDSPmatch::analyze(art::Event const& e){
  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  
  fCRT_time.clear();
  fCRTChan.clear();
  
  //fCTB_time.clear();
  fCTBChan.clear();
  
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
  
  fPM_time.clear();
  fTrkStartx_PM.clear();
  fTrkStarty_PM.clear();
  fTrkStartz_PM.clear();
  fTrkEndx_PM.clear();
  fTrkEndy_PM.clear();
  fTrkEndz_PM.clear();
 

  const auto timeStamps = e.getValidHandle<std::vector<raw::RDTimeStamp>>(fTimeLabel);
  int trigger = -1;
  if(timeStamps.isValid() && timeStamps->size() == 1){
    const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
    trigger = timeStamp.GetFlags();
    timing_time = timeStamp.GetTimeStamp();
    std::cout << "Timing: " << std::setw(20) << timing_time << std::endl;
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
  
  //bool us = false;
  //bool ds = false;
  
  for(const auto& status: statuses){
    fCRT_time.clear();
    fCRTChan.clear();
    
    //fCTB_time.clear();
    fCTBChan.clear();
    
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
   
    fPM_time.clear();
    fTrkStartx_PM.clear();
    fTrkStarty_PM.clear();
    fTrkStartz_PM.clear();
    fTrkEndx_PM.clear();
    fTrkEndy_PM.clear();
    fTrkEndz_PM.clear();
   
    if(status.crt != 0) {
      std::vector<CRT::Trigger> inWindow(crtHandle->size());
      const int64_t ctbetime = status.timestamp;
      const auto newEnd = std::copy_if(crtHandle->begin(), crtHandle->end(), inWindow.begin(),[this,ctbetime](const auto& trigger){
	return (uint64_t)(labs(ctbetime - ((int64_t)(trigger.Timestamp()) - fCRTCTBOffset))) < fCRTWindow; });
      inWindow.resize(std::distance(inWindow.begin(), newEnd));   
      const std::string binary = std::bitset<32>(status.crt).to_string();
      for(std::string::size_type l = binary.length();l>0;l--) {
	if(binary[l]=='1') {
	  std::cout << "CTB: " << status.timestamp << "fCTBChan: " << l << std::endl;
	  fCTB_time = status.timestamp;
	  fCTBChan.push_back(binary.length()-l);
	}
      }
      for(size_t k=0;k<inWindow.size();++k){
	std::cout << "CRT: " << inWindow[k].Timestamp() << "fCRTChan: " << inWindow[k].Channel() << std::endl;
	fCRT_time.push_back(inWindow[k].Timestamp());
	fCRTChan.push_back(inWindow[k].Channel());
      }
    }   
    for(const auto& OpHit: *OpHitHandle){
      // if((timing_time-OpHit.PeakTime == 112){
      std::cout << "SSP Peaktime: " << static_cast<int64_t>(OpHit.PeakTime())/3 << " SSP Peaktime Abs: " << static_cast<int64_t>(OpHit.PeakTimeAbs())/3 << " OpChannel: " << OpHit.OpChannel() << " PE: " << OpHit.PE() <<  std::endl;
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
      if((Trk->Vertex().Z() > 10) || (Trk->End().Z() < 650)) continue;
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
    }
    
    art::Handle<std::vector<recob::Track>> PMTrkHandle;
    std::vector<art::Ptr<recob::Track>> PMTrks;
    if(e.getByLabel(fPMLabel,PMTrkHandle)) art::fill_ptr_vector(PMTrks,PMTrkHandle);
    
    art::FindManyP<anab::T0> PPMT0(PMTrkHandle,e,fPMLabel);
    
    for(size_t r = 0;r<PMTrks.size();++r){
      auto & Trk = PMTrks[r];
      if((Trk->Vertex().Z() > 10) || (Trk->End().Z() < 650)) continue;
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
    }
  }
      
  
    //fCRTChan.push_back(inWindow[k].Channel());
    //fCTBChan.push_back(status.crt);
  fTree->Fill();
}
  

  //ctb_time = status.timestamp;

  //Get CTB Stuff
  //const int64_t ctbtime = status.timestamp;
 


DEFINE_ART_MODULE(pdsp::PDSPmatch)
