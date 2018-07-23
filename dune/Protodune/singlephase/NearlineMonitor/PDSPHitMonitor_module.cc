/////////////////////////////////////////////
// Hit Monitor Module                      //
// June 2018                               //
// georgios.christodoulou <at> cern.ch     //
/////////////////////////////////////////////

#ifndef PDSPHitMonitor_module
#define PDSPHitMonitor_module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/DataViewImpl.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

// Data type includes
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"

// ROOT includes
#include "TTree.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TProfile.h"

// C++ Includes
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace PDSPHitmonitor_module{
  
  class PDSPHitMonitorModule : public art::EDAnalyzer{
  public:
    
    explicit PDSPHitMonitorModule(fhicl::ParameterSet const& pset);
    virtual ~PDSPHitMonitorModule();
    
    void beginJob();
    void analyze(const art::Event& evt);
    void reconfigure(fhicl::ParameterSet const & p);
    
  private:
    art::InputTag fTPCHitTag;

    unsigned int fChansPerAPA;
    unsigned int fNofAPA;

    unsigned int fUChanMin;
    unsigned int fVChanMin;
    unsigned int fUChanMax;
    unsigned int fVChanMax;
    unsigned int fZ0ChanMax;
    unsigned int fZ0ChanMin;
    unsigned int fZ1ChanMin;
    unsigned int fZ1ChanMax;

    // Histograms
    TH1I *fTotalNHits;
    TH1F *fHitCharge;
    TH1F *fHitRMS;
    TH1F *fHitPeakTime;

    std::vector<TH1I*> fNHitsAPAViewU;
    std::vector<TH1I*> fNHitsAPAViewV;
    std::vector<TH1I*> fNHitsAPAViewZ;

    std::vector<TH1F*> fChargeAPAViewU;
    std::vector<TH1F*> fChargeAPAViewV;
    std::vector<TH1F*> fChargeAPAViewZ;

    std::vector<TH1F*> fRMSAPAViewU;
    std::vector<TH1F*> fRMSAPAViewV;
    std::vector<TH1F*> fRMSAPAViewZ;

    std::vector<TH1F*> fHitPeakTimeAPAViewU;
    std::vector<TH1F*> fHitPeakTimeAPAViewV;
    std::vector<TH1F*> fHitPeakTimeAPAViewZ;

    // Profiles
    std::vector<TProfile*> fNHitsAPAViewU_prof;
    std::vector<TProfile*> fNHitsAPAViewV_prof;
    std::vector<TProfile*> fNHitsAPAViewZ_prof;

    std::vector<TProfile*> fChargeAPAViewU_prof;
    std::vector<TProfile*> fChargeAPAViewV_prof;
    std::vector<TProfile*> fChargeAPAViewZ_prof;

    std::vector<TProfile*> fRMSAPAViewU_prof;
    std::vector<TProfile*> fRMSAPAViewV_prof;
    std::vector<TProfile*> fRMSAPAViewZ_prof;

    geo::GeometryCore const * fGeom = &*(art::ServiceHandle<geo::Geometry>());

    std::vector<unsigned int> fApaLabelNum;
    
  };
  
  //-----------------------------------------------------------------------
  void PDSPHitMonitorModule::reconfigure(fhicl::ParameterSet const & p){

    fTPCHitTag = p.get<art::InputTag>("TPCHitTag", "a:b:c");

  }

  //-----------------------------------------------------------------------
  PDSPHitMonitorModule::PDSPHitMonitorModule(fhicl::ParameterSet const& pset) : EDAnalyzer(pset), fApaLabelNum{3,5,2,6,1,4} {
    this->reconfigure(pset);
  }
  
  //-----------------------------------------------------------------------
  PDSPHitMonitorModule::~PDSPHitMonitorModule(){}
  
  //-----------------------------------------------------------------------
  void PDSPHitMonitorModule::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    // Accquiring geometry data
    fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2;
    fChansPerAPA = fGeom->Nchannels()/fNofAPA;
    mf::LogVerbatim("HitMonitor") << " NTPCs = " << fGeom->NTPC() << " , Ncryostats = " << fGeom->Ncryostats() << " , NofAPA = " << fNofAPA << " , ChansPerAPA = " << fChansPerAPA << std::endl;

    // loop through channels in the first APA to find the channel boundaries for each view
    // will adjust for desired APA after
    fUChanMin = 0;
    fVChanMin = 0;
    fUChanMax = 0;
    fVChanMax = 0;
    fZ0ChanMax = 0;
    fZ0ChanMin = 0;
    fZ1ChanMin = 0;
    fZ1ChanMax = fChansPerAPA - 1;
    for(unsigned int c = fUChanMin + 1; c < fZ1ChanMax; c++){
      if(fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU){
        fVChanMin = c;
        fUChanMax = c - 1;
      }
      if(fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV){
        fZ0ChanMin = c;
        fVChanMax = c-1;
      }
      if(fGeom->View(c) == geo::kZ && fGeom->ChannelToWire(c)[0].TPC == fGeom->ChannelToWire(c-1)[0].TPC + 1){
        fZ1ChanMin = c;
        fZ0ChanMax = c-1;
      }
    }

    //unsigned int fNUCh=fUChanMax-fUChanMin+1;
    //unsigned int fNVCh=fVChanMax-fVChanMin+1;
    //unsigned int fNZ0Ch=fZ0ChanMax-fZ0ChanMin+1;
    //unsigned int fNZ1Ch=fZ1ChanMax-fZ1ChanMin+1;

    //mf::LogVerbatim("HitMonitor") << " U: "<< fNUCh <<"  V:  "<< fNVCh << "  Z0:  "<< fNZ0Ch << "  Z1:  " << fNZ1Ch << std::endl;
    
    for(unsigned int i=0;i<fNofAPA;i++){

      unsigned int j=fApaLabelNum.at(i);

      unsigned int UChMin=fUChanMin + i*fChansPerAPA;
      unsigned int UChMax=fUChanMax + i*fChansPerAPA;
      unsigned int VChMin=fVChanMin + i*fChansPerAPA;
      unsigned int VChMax=fVChanMax + i*fChansPerAPA;
      unsigned int ZChMin=fZ0ChanMin + i*fChansPerAPA;
      //unsigned int ZChMax=fZ0ChanMax + i*fChansPerAPA;
      unsigned int ZChMax=fZ1ChanMax + i*fChansPerAPA; //including unused channels

      //std::cout<< "UCh: " << UChMin << " - " << UChMax << std::endl;
      //std::cout<< "VCh: " << VChMin << " - " << VChMax << std::endl;
      //std::cout<< "ZCh: " << ZChMin << " - " << ZChMax << std::endl;

      fNHitsAPAViewU.push_back(tfs->make<TH1I>(Form("NHitsAPA%d_U",j),Form("Number of hits APA%d-U",j),1000,0,20000));
      fNHitsAPAViewV.push_back(tfs->make<TH1I>(Form("NHitsAPA%d_V",j),Form("Number of hits APA%d-V",j),1000,0,20000));
      fNHitsAPAViewZ.push_back(tfs->make<TH1I>(Form("NHitsAPA%d_Z",j),Form("Number of hits APA%d-Z",j),1000,0,20000));

      fChargeAPAViewU.push_back(tfs->make<TH1F>(Form("HitChargeAPA%d_U",j),Form("Hit Charge APA%d-U",j),100,0,1500));
      fChargeAPAViewV.push_back(tfs->make<TH1F>(Form("HitChargeAPA%d_V",j),Form("Hit Charge APA%d-V",j),100,0,1500));
      fChargeAPAViewZ.push_back(tfs->make<TH1F>(Form("HitChargeAPA%d_Z",j),Form("Hit Charge APA%d-Z",j),100,0,1500));

      fRMSAPAViewU.push_back(tfs->make<TH1F>(Form("HitRMSAPA%d_U",j),Form("Hit RMS APA%d-U",j),100,0,10));
      fRMSAPAViewV.push_back(tfs->make<TH1F>(Form("HitRMSAPA%d_V",j),Form("Hit RMS APA%d-V",j),100,0,10));
      fRMSAPAViewZ.push_back(tfs->make<TH1F>(Form("HitRMSAPA%d_Z",j),Form("Hit RMS APA%d-Z",j),100,0,10));

      fHitPeakTimeAPAViewU.push_back(tfs->make<TH1F>(Form("HitPeakTimeAPA%d_U",j),Form("Hit Peak Time APA%d-U",j),100,0,10000));
      fHitPeakTimeAPAViewV.push_back(tfs->make<TH1F>(Form("HitPeakTimeAPA%d_V",j),Form("Hit Peak Time APA%d-V",j),100,0,10000));
      fHitPeakTimeAPAViewZ.push_back(tfs->make<TH1F>(Form("HitPeakTimeAPA%d_Z",j),Form("Hit Peak Time APA%d-Z",j),100,0,10000));

      // U view
      fNHitsAPAViewU_prof.push_back(tfs->make<TProfile>(Form("fNHitsViewU%d_prof", j),Form("Number of hits distribution vs Channel(Plane U, APA%d)", j),  UChMax - UChMin + 1, UChMin-0.5, UChMax+0.5, "s"));
      fChargeAPAViewU_prof.push_back(tfs->make<TProfile>(Form("fChargeViewU%d_prof", j),Form("Profiled Hit Charge distribution vs Channel(Plane U, APA%d)", j),  UChMax - UChMin + 1, UChMin-0.5, UChMax+0.5, "s"));
      fRMSAPAViewU_prof.push_back(tfs->make<TProfile>(Form("fRMSViewU%d_prof", j),Form("Profiled Hit RMS distribution vs Channel(Plane U, APA%d)", j),  UChMax - UChMin + 1, UChMin-0.5, UChMax+0.5, "s"));
      
      // V view
      fNHitsAPAViewV_prof.push_back(tfs->make<TProfile>(Form("fNHitsViewV%d_prof",j),Form("Number of hits distribution vs Channel(Plane V, APA%d)",j), VChMax - VChMin + 1, VChMin-0.5, VChMax+0.5, "s"));
      fChargeAPAViewV_prof.push_back(tfs->make<TProfile>(Form("fChargeViewV%d_prof",j),Form("Profiled Hit Charge distribution vs Channel(Plane V, APA%d)",j), VChMax - VChMin + 1, VChMin-0.5, VChMax+0.5, "s"));
      fRMSAPAViewV_prof.push_back(tfs->make<TProfile>(Form("fRMSViewV%d_prof",j),Form("Profiled Hit RMS distribution vs Channel(Plane V, APA%d)",j), VChMax - VChMin + 1, VChMin-0.5, VChMax+0.5, "s"));

      // Z view
      fNHitsAPAViewZ_prof.push_back(tfs->make<TProfile>(Form("fNHitsViewZ%d_prof",j),Form("Number of hits distribution vs Channel(Plane Z, APA%d)",j), ZChMax - ZChMin + 1, ZChMin-0.5, ZChMax+0.5, "s"));
      fChargeAPAViewZ_prof.push_back(tfs->make<TProfile>(Form("fChargeViewZ%d_prof",j),Form("Profiled Hit Charge distribution vs Channel(Plane Z, APA%d)",j),  ZChMax - ZChMin + 1, ZChMin-0.5, ZChMax+0.5, "s"));
      fRMSAPAViewZ_prof.push_back(tfs->make<TProfile>(Form("fRMSViewZ%d_prof",j),Form("Profiled Hit RMS distribution vs Channel(Plane Z, APA%d)",j),  ZChMax - ZChMin + 1, ZChMin-0.5, ZChMax+0.5, "s"));

      // Set titles
      fNHitsAPAViewU[i]->GetXaxis()->SetTitle("NHits");
      fNHitsAPAViewV[i]->GetXaxis()->SetTitle("NHits");
      fNHitsAPAViewZ[i]->GetXaxis()->SetTitle("NHits");

      fChargeAPAViewU[i]->GetXaxis()->SetTitle("Charge (ADC units)");
      fChargeAPAViewV[i]->GetXaxis()->SetTitle("Charge (ADC units)");
      fChargeAPAViewZ[i]->GetXaxis()->SetTitle("Charge (ADC units)");

      fRMSAPAViewU[i]->GetXaxis()->SetTitle("Hit RMS");
      fRMSAPAViewV[i]->GetXaxis()->SetTitle("Hit RMS");
      fRMSAPAViewZ[i]->GetXaxis()->SetTitle("Hit RMS");

      fHitPeakTimeAPAViewU[i]->GetXaxis()->SetTitle("Hit Peak Time");
      fHitPeakTimeAPAViewV[i]->GetXaxis()->SetTitle("Hit Peak Time");
      fHitPeakTimeAPAViewZ[i]->GetXaxis()->SetTitle("Hit Peak Time");

      fNHitsAPAViewU_prof[i]->GetXaxis()->SetTitle("Channel ID"); fNHitsAPAViewU_prof[i]->GetYaxis()->SetTitle("NHits");
      fNHitsAPAViewV_prof[i]->GetXaxis()->SetTitle("Channel ID"); fNHitsAPAViewV_prof[i]->GetYaxis()->SetTitle("NHits");
      fNHitsAPAViewZ_prof[i]->GetXaxis()->SetTitle("Channel ID"); fNHitsAPAViewZ_prof[i]->GetYaxis()->SetTitle("NHits");

      fChargeAPAViewU_prof[i]->GetXaxis()->SetTitle("Channel ID"); fChargeAPAViewU_prof[i]->GetYaxis()->SetTitle("Charge (ADC units)");
      fChargeAPAViewV_prof[i]->GetXaxis()->SetTitle("Channel ID"); fChargeAPAViewV_prof[i]->GetYaxis()->SetTitle("Charge (ADC units)");
      fChargeAPAViewZ_prof[i]->GetXaxis()->SetTitle("Channel ID"); fChargeAPAViewZ_prof[i]->GetYaxis()->SetTitle("Charge (ADC units)");

      fRMSAPAViewU_prof[i]->GetXaxis()->SetTitle("Channel ID"); fRMSAPAViewU_prof[i]->GetYaxis()->SetTitle("Hit RMS");
      fRMSAPAViewV_prof[i]->GetXaxis()->SetTitle("Channel ID"); fRMSAPAViewV_prof[i]->GetYaxis()->SetTitle("Hit RMS");
      fRMSAPAViewZ_prof[i]->GetXaxis()->SetTitle("Channel ID"); fRMSAPAViewZ_prof[i]->GetYaxis()->SetTitle("Hit RMS");
    }

    // Summary histograms from all APAs
    fTotalNHits   = tfs->make<TH1I>("fTotalNHits"   ,"Total number of hits" ,1000,0,1000000);
    fHitCharge    = tfs->make<TH1F>("fHitCharge"    ,"Hit Charge"           ,100 ,0,1500);
    fHitRMS       = tfs->make<TH1F>("fHitChargeRMS" ,"Hit RMS"              ,100 ,0,10);
    fHitPeakTime  = tfs->make<TH1F>("fHitPeakTime"  ,"Hit Peak Time"        ,100 ,0,10000);

    // Set titles
    fTotalNHits->GetXaxis()->SetTitle("NHits");
    fHitCharge->GetXaxis()->SetTitle("Charge (ADC units)");
    fHitRMS->GetXaxis()->SetTitle("Hit RMS");
    fHitPeakTime->GetXaxis()->SetTitle("Hit Peak Time");

  }
  
  //-----------------------------------------------------------------------
  void PDSPHitMonitorModule::analyze(const art::Event& evt){

    art::Handle<std::vector<recob::Hit> > hitHandle;
    // Many options here: gaushit, hitfd, linecluster, lineclusterdc, robusthit, fasthit, etc.
    bool retVal = evt.getByLabel(fTPCHitTag, hitHandle);
    if(retVal==true)
      ;
    else{
      mf::LogWarning("HitMonitor") << " Getting Hits FAIL: " << fTPCHitTag << std::endl;
      return;
    }

    std::vector<art::Ptr<recob::Hit> > hitlist;
    art::fill_ptr_vector(hitlist, hitHandle);

    int NHits = hitlist.size();
    
    // Arrays ot store number of hits per APA
    int NHitsPerApa[10][3];
    for(unsigned int j=0; j<10; j++){
      for(int k=0; k<3; k++){
	NHitsPerApa[j][k] = 0;
      }
    }

    int NHitsChannel[20000];
    for(unsigned int j=0; j<20000; j++){
      NHitsChannel[j] = 0;
    }
    
    mf::LogVerbatim("HitMonitor") << " Number of hits = " << NHits << std::endl;

    fTotalNHits->Fill(NHits);
    for(int i=0; i<NHits; i++){
      art::Ptr<recob::Hit> hit = hitlist[i];

      int hit_channel = hit->Channel();
      unsigned int apa = std::floor(hit->WireID().TPC/2);
      //unsigned int apa = std::floor(hit_channel/fChansPerAPA);

      // Protection
      if(apa > fNofAPA){
	mf::LogWarning("HitMonitor") << "APA number found (" << apa << ") larger than maximum (" << fNofAPA << "). Skipping hit!" << std::endl;
	continue;
      }

      float hit_charge = hit->Integral();
      //float hit_chargeSigma = hit->SigmaIntegral();
      float hit_rms = hit->RMS();
      float hit_peakT = hit->PeakTime();
      //float hit_wireID = hit->WireID().Wire;
      int hit_plane = (int)hit->WireID().Plane;

      fHitCharge->Fill(hit_charge);
      fHitRMS->Fill(hit_rms);
      fHitPeakTime->Fill(hit_peakT);

      NHitsPerApa[apa][hit_plane]++;
      NHitsChannel[hit_channel]++;

      if(hit_plane == 0){
	fChargeAPAViewU[apa]->Fill(hit_charge);
	fRMSAPAViewU[apa]->Fill(hit_rms);
	fHitPeakTimeAPAViewU[apa]->Fill(hit_peakT);

	fChargeAPAViewU_prof[apa]->Fill(hit_channel, hit_charge, 1);
	fRMSAPAViewU_prof[apa]->Fill(hit_channel, hit_rms, 1);
      }
      else if(hit_plane == 1){
        fChargeAPAViewV[apa]->Fill(hit_charge);
	fRMSAPAViewV[apa]->Fill(hit_rms);
	fHitPeakTimeAPAViewV[apa]->Fill(hit_peakT);

	fChargeAPAViewV_prof[apa]->Fill(hit_channel, hit_charge, 1);
	fRMSAPAViewV_prof[apa]->Fill(hit_channel, hit_rms, 1);
      }
      else if(hit_plane == 2){
        fChargeAPAViewZ[apa]->Fill(hit_charge);
	fRMSAPAViewZ[apa]->Fill(hit_rms);
	fHitPeakTimeAPAViewZ[apa]->Fill(hit_peakT);

	fChargeAPAViewZ_prof[apa]->Fill(hit_channel, hit_charge, 1);
	fRMSAPAViewZ_prof[apa]->Fill(hit_channel, hit_rms, 1);
      }
    }
    
    // Now fill th number of hits histograms
    for(unsigned int j=0; j<fNofAPA; j++){
      fNHitsAPAViewU[j]->Fill(NHitsPerApa[j][0]);
      fNHitsAPAViewV[j]->Fill(NHitsPerApa[j][1]);
      fNHitsAPAViewZ[j]->Fill(NHitsPerApa[j][2]);

      unsigned int UChMin=fUChanMin + j*fChansPerAPA;
      unsigned int UChMax=fUChanMax + j*fChansPerAPA;
      unsigned int VChMin=fVChanMin + j*fChansPerAPA;
      unsigned int VChMax=fVChanMax + j*fChansPerAPA;
      unsigned int ZChMin=fZ0ChanMin + j*fChansPerAPA;
      //unsigned int ZChMax=fZ0ChanMax + j*fChansPerAPA;
      unsigned int ZChMax=fZ1ChanMax + j*fChansPerAPA; //including unused channels

      for(unsigned int k=UChMin; k<UChMax; k++){
	fNHitsAPAViewU_prof[j]->Fill(k, NHitsChannel[k], 1);	
      }

      for(unsigned int k=VChMin; k<VChMax; k++){
        fNHitsAPAViewV_prof[j]->Fill(k, NHitsChannel[k], 1);
      }

      for(unsigned int k=ZChMin; k<ZChMax; k++){
        fNHitsAPAViewZ_prof[j]->Fill(k, NHitsChannel[k], 1);
      }
    }

  }
 
} // namespace

DEFINE_ART_MODULE(PDSPHitmonitor_module::PDSPHitMonitorModule)

#endif
