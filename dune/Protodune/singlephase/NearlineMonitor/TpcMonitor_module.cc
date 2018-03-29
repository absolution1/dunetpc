///////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Class:       TpcMonitor_module
// Module type: analyzer
// File:        TpcMonitor_module.cc
// Author:      Jingbo Wang (jiowang@ucdavis.edu), February 2018
//
///////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef TpcMonitor_module
#define TpcMonitor_module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Data type includes
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

// ROOT includes.
#include "TH1.h"
#include "TH2.h"
#include "TH2F.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TProfile.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<int> >+;
#endif

namespace tpc_monitor{

  class TpcMonitor : public art::EDAnalyzer{
  	
  public:
 
    explicit TpcMonitor(fhicl::ParameterSet const& pset);
    virtual ~TpcMonitor();

    void beginJob();
    void beginRun(const art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt); 

    void endJob();
    
  private:
    // Parameters in .fcl file
    std::string fRawDigitLabel;
    std::string fTPCInput;
    std::string fTPCInstance;

    // Branch variables for tree
    unsigned int fEvent;
    unsigned int fRun;
    unsigned int fSubRun;

    int NumberOfRCEs;

    // TPC
    unsigned int fNUCh;
    unsigned int fNVCh;
    unsigned int fNZ0Ch;
    unsigned int fNZ1Ch;
    // find channel boundaries for each view
    unsigned int fUChanMin;
    unsigned int fUChanMax;
    unsigned int fVChanMin;
    unsigned int fVChanMax;
    unsigned int fZ0ChanMin;
    unsigned int fZ0ChanMax;
    unsigned int fZ1ChanMin;
    unsigned int fZ1ChanMax;
    unsigned int fNticks;

    unsigned int fNofAPA;
    unsigned int fChansPerAPA;
    
    // sampling rate
    float fSampleRate;
    
    // bin width in kHz
    float fBinWidth;

    // Mean/RMS by offline channel
    std::vector<TH2F*> fChanRMSU;
    std::vector<TH2F*> fChanMeanU;
    std::vector<TH2F*> fChanFFTU;

    std::vector<TH2F*> fChanRMSV;
    std::vector<TH2F*> fChanMeanV;
    std::vector<TH2F*> fChanFFTV;

    std::vector<TH2F*> fChanRMSZ;
    std::vector<TH2F*> fChanMeanZ;
    std::vector<TH2F*> fChanFFTZ;
    	
    // profiled Mean/RMS by offline channel
    std::vector<TProfile*> fChanRMSU_pfx;
    std::vector<TProfile*> fChanMeanU_pfx;
    std::vector<TProfile*> fChanRMSV_pfx;
    std::vector<TProfile*> fChanMeanV_pfx;
    std::vector<TProfile*> fChanRMSZ_pfx;
    std::vector<TProfile*> fChanMeanZ_pfx;
    	
    // Mean/RMS by slot number
    std::vector<TH2F*> fSlotChanMean;
    std::vector<TH2F*> fSlotChanRMS;	
    std::vector<TProfile*> fSlotChanMean_pfx;
    std::vector<TProfile*> fSlotChanRMS_pfx;	
    
    // FFT by slot number
    std::vector<TH2F*> fSlotChanFFT;
    
    // Persistent and overlay wavefroms by fiber
    //std::vector<TH2I*> fPersistentWav_by_Fiber;
    std::vector<TH2F*> fPersistentFFT_by_Fiber;
    
    // Profiled fft by fiber
    std::vector<TProfile*> fFFT_by_Fiber_pfx;
    	
    // Rebin factor
    int fRebinX;
    int fRebinY; 

    // define nADC counts for uncompressed vs compressed
    unsigned int nADC_comp;
    unsigned int nADC_uncomp;
    unsigned int nADC_uncompPed;

    // define my RMS adcs in 3 views
    float fRMSadc_U[5];
    float fRMSadc_V[5];
    float fRMSadc_Z[5];
    float fMEANadc_U[5];
    float fMEANadc_V[5];
    float fMEANadc_Z[5];

    float maxMean_U[5];
    float maxRMS_U[5];

    float maxMean_V[5];
    float maxRMS_V[5];

    float maxMean_Z[5];
    float maxRMS_Z[5];

    // define functions
    float rmsADC(std::vector< short > & uncompressed);
    float meanADC(std::vector< short > & uncompressed);
    void calculateFFT(TH1D* hist_waveform, TH1D* graph_frequency);
    geo::GeometryCore const * fGeom = &*(art::ServiceHandle<geo::Geometry>());


 }; // TpcMonitor

  //-----------------------------------------------------------------------

  TpcMonitor::TpcMonitor(fhicl::ParameterSet const& parameterSet)
  : EDAnalyzer(parameterSet), fRebinX(1), fRebinY(1) {
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------


  TpcMonitor::~TpcMonitor() {}
   

  //-----------------------------------------------------------------------
  
  void TpcMonitor::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    unsigned int UChMin;
    unsigned int UChMax;
    unsigned int VChMin;
    unsigned int VChMax;
    unsigned int ZChMin;
    unsigned int ZChMax;

    // Accquiring geometry data
    fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2;
    fChansPerAPA = fGeom->Nchannels()/fNofAPA;

    // loop through channels in the first APA to find the channel boundaries for each view
    // will adjust for desired APA after
    fUChanMin = 0;
    fZ1ChanMax = fChansPerAPA - 1;
    for ( unsigned int c = fUChanMin + 1; c < fZ1ChanMax; c++ ){
      if ( fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU ){
        fVChanMin = c;
        fUChanMax = c - 1;
      }
      if ( fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV ){
        fZ0ChanMin = c;
        fVChanMax = c-1;
      }
      if ( fGeom->View(c) == geo::kZ && fGeom->ChannelToWire(c)[0].TPC == fGeom->ChannelToWire(c-1)[0].TPC + 1 ){
        fZ1ChanMin = c;
        fZ0ChanMax = c-1;
      }
    }
    

    fNUCh=fUChanMax-fUChanMin+1;
    fNVCh=fVChanMax-fVChanMin+1;
    fNZ0Ch=fZ0ChanMax-fZ0ChanMin+1;
    fNZ1Ch=fZ1ChanMax-fZ1ChanMin+1;

    mf::LogVerbatim("TpcMonitor")
      <<"U: "<< fNUCh<<"  V:  "<<fNVCh<<"  Z0:  "<<fNZ0Ch << "  Z1:  " <<fNZ1Ch << std::endl;
    
    //Mean/RMS by offline channel for each view in each APA
    for(unsigned int i=0;i<fNofAPA;i++){
      UChMin=fUChanMin + i*fChansPerAPA;
      UChMax=fUChanMax + i*fChansPerAPA;      
      VChMin=fVChanMin + i*fChansPerAPA;
      VChMax=fVChanMax + i*fChansPerAPA;      
      ZChMin=fZ0ChanMin + i*fChansPerAPA;
      //ZChMax=fZ0ChanMax + i*fChansPerAPA;
      ZChMax=fZ1ChanMax + i*fChansPerAPA; //including unused channels
      
      std::cout<<"UCh:"<<UChMin<<" - "<<UChMax<<std::endl;
      std::cout<<"VCh:"<<VChMin<<" - "<<VChMax<<std::endl;
      std::cout<<"ZCh:"<<ZChMin<<" - "<<ZChMax<<std::endl;
      
      // U view
      fChanRMSU.push_back(tfs->make<TH2F>(Form("fChanRMSU%d", i),Form("raw-ped RMS vs Channel(Plane U, APA%d)", i),  UChMax - UChMin + 1, UChMin, UChMax, 625,0,5000));      
      fChanRMSU_pfx.push_back(tfs->make<TProfile>(Form("fChanRMSU%d_pfx", i),Form("Profiled raw-ped RMS vs Channel(Plane U, APA%d)", i),  UChMax - UChMin + 1, UChMin, UChMax)); 
      fChanMeanU.push_back(tfs->make<TH2F>(Form("fChanMeanU%d",i),Form("raw-ped MEAN vs Channel(Plane U, APA%d)",i),  UChMax - UChMin + 1, UChMin, UChMax, 625,0,5000));  
      fChanMeanU_pfx.push_back(tfs->make<TProfile>(Form("fChanMeanU%d_pfx",i),Form("Profiled raw-ped MEAN vs Channel(Plane U, APA%d)",i),  UChMax - UChMin + 1, UChMin, UChMax)); 
      fChanFFTU.push_back(tfs->make<TH2F>(Form("fChanFFTU%d", i),Form("fChanFFT (Plane U, APA%d)", i),  UChMax - UChMin + 1, UChMin, UChMax, fNticks/2,0,fNticks/2*fBinWidth));
      
      // V view
      fChanRMSV.push_back(tfs->make<TH2F>(Form("fChanRMSV%d",i),Form("raw-ped RMS vs Channel(Plane V, APA%d)",i), VChMax - VChMin + 1, VChMin, VChMax, 625,0,5000));
      fChanRMSV_pfx.push_back(tfs->make<TProfile>(Form("fChanRMSV%d_pfx",i),Form("Profiled raw-ped RMS vs Channel(Plane V, APA%d)",i),  VChMax - VChMin + 1, VChMin, VChMax)); 
      fChanMeanV.push_back(tfs->make<TH2F>(Form("fChanMeanV%d",i),Form("raw-ped Mean vs Channel(Plane V, APA%d)",i),  VChMax - VChMin + 1, VChMin, VChMax, 625,0,5000));
      fChanMeanV_pfx.push_back(tfs->make<TProfile>(Form("fChanMeanV%d_pfx",i),Form("Profiled raw-ped Mean vs Channel(Plane V, APA%d)",i),  VChMax - VChMin + 1, VChMin, VChMax));   
      fChanFFTV.push_back(tfs->make<TH2F>(Form("fChanFFTV%d", i),Form("fChanFFT (Plane V, APA%d)", i),  VChMax - VChMin + 1, VChMin, VChMax, fNticks/2,0,fNticks/2*fBinWidth));
      
      // Z view                                                                                                                                                           
      fChanRMSZ.push_back(tfs->make<TH2F>(Form("fChanRMSZ%d",i),Form("raw-ped RMS vs Channel(Plane Z, APA%d)",i),  ZChMax - ZChMin + 1, ZChMin, ZChMax, 625,0,5000));
      fChanRMSZ_pfx.push_back(tfs->make<TProfile>(Form("fChanRMSZ%d_pfx",i),Form("Profiled raw-ped RMS vs Channel(Plane Z, APA%d)",i),  ZChMax - ZChMin + 1, ZChMin, ZChMax)); 
      fChanMeanZ.push_back(tfs->make<TH2F>(Form("fChanMeanZ%d",i),Form("raw-ped Mean vs Channel(Plane Z, APA%d)",i),  ZChMax - ZChMin + 1, ZChMin, ZChMax, 625,0,5000));
      fChanMeanZ_pfx.push_back(tfs->make<TProfile>(Form("fChanMeanZ%d_pfx",i),Form("Profiled raw-ped Mean vs Channel(Plane Z, APA%d)",i),  ZChMax - ZChMin + 1, ZChMin, ZChMax)); 
      fChanFFTZ.push_back(tfs->make<TH2F>(Form("fChanFFTZ%d", i),Form("fChanFFT (Plane Z, APA%d)", i),  ZChMax - ZChMin + 1, ZChMin, ZChMax, fNticks/2,0,fNticks/2*fBinWidth));
      
      // Set titles
      fChanRMSU[i]->GetXaxis()->SetTitle("Chan"); fChanRMSU[i]->GetYaxis()->SetTitle("raw RMS"); 
      fChanRMSU_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanRMSU_pfx[i]->GetYaxis()->SetTitle("raw RMS"); 
      fChanMeanU[i]->GetXaxis()->SetTitle("Chan"); fChanMeanU[i]->GetYaxis()->SetTitle("raw Mean");
      fChanMeanU_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanMeanU_pfx[i]->GetYaxis()->SetTitle("raw Mean"); 
      fChanRMSV[i]->GetXaxis()->SetTitle("Chan"); fChanRMSV[i]->GetYaxis()->SetTitle("raw RMS"); 
      fChanRMSV_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanRMSV_pfx[i]->GetYaxis()->SetTitle("raw RMS"); 
      fChanMeanV[i]->GetXaxis()->SetTitle("Chan"); fChanMeanV[i]->GetYaxis()->SetTitle("raw Mean");
      fChanMeanV_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanMeanV_pfx[i]->GetYaxis()->SetTitle("raw Mean"); 
      fChanRMSZ[i]->GetXaxis()->SetTitle("Chan"); fChanRMSZ[i]->GetYaxis()->SetTitle("raw RMS"); 
      fChanRMSZ_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanRMSZ_pfx[i]->GetYaxis()->SetTitle("raw RMS"); 
      fChanMeanZ[i]->GetXaxis()->SetTitle("Chan"); fChanMeanZ[i]->GetYaxis()->SetTitle("raw Mean"); 
      fChanMeanZ_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanMeanZ_pfx[i]->GetYaxis()->SetTitle("raw Mean"); 
      fChanFFTU[i]->GetXaxis()->SetTitle("Chan"); fChanFFTU[i]->GetYaxis()->SetTitle("kHz"); 
      fChanFFTV[i]->GetXaxis()->SetTitle("Chan"); fChanFFTV[i]->GetYaxis()->SetTitle("kHz"); 
      fChanFFTZ[i]->GetXaxis()->SetTitle("Chan"); fChanFFTZ[i]->GetYaxis()->SetTitle("kHz"); 
      
      //  Rebin histograms
      //std::cout<<"RebinX = "<<fRebinX<<"  RebinY = "<<fRebinY<<std::endl;
      fChanRMSU[i]->Rebin2D(fRebinX, fRebinY);
      fChanMeanU[i]->Rebin2D(fRebinX, fRebinY);
      fChanRMSV[i]->Rebin2D(fRebinX, fRebinY);
      fChanRMSZ[i]->Rebin2D(fRebinX, fRebinY);
      fChanMeanZ[i]->Rebin2D(fRebinX, fRebinY);
      fChanFFTU[i]->Rebin2D(fRebinX, fRebinY);
      fChanFFTV[i]->Rebin2D(fRebinX, fRebinY);
      fChanFFTZ[i]->Rebin2D(fRebinX, fRebinY);
    }
    
    // Mean/RMS by slot channel number for each slot
    for(int i=0;i<30;i++) {
  	  fSlotChanMean.push_back(tfs->make<TH2F>(Form("Slot%d_Mean", i), Form("Slot%d:Mean_vs_SlotChannel", i), 512, 0, 512, 625, .0, 5000)); 
  	  fSlotChanRMS.push_back(tfs->make<TH2F>(Form("Slot%d_RMS", i), Form("Slot%d:RMS_vs_SlotChannel", i), 512, 0, 512, 625, .0, 5000)); 
  	  fSlotChanMean_pfx.push_back(tfs->make<TProfile>(Form("Slot%d_Mean_pfx", i), Form("Slot%d:Mean_vs_SlotChannel_pfx", i), 512, 0, 512)); 
  	  fSlotChanRMS_pfx.push_back(tfs->make<TProfile>(Form("Slot%d_RMS_pfx", i), Form("Slot%d:RMS_vs_SlotChannel_pfx", i), 512, 0, 512)); 
  	  fSlotChanFFT.push_back(tfs->make<TH2F>(Form("Slot%d_FFT", i), Form("Slot%d:FFT_vs_SlotChannel", i), 512, 0, 512, fNticks/2, 0, fNticks/2*fBinWidth));
  	  
  	  fSlotChanMean[i]->GetXaxis()->SetTitle("Slot Channel"); fSlotChanMean[i]->GetYaxis()->SetTitle("Raw Mean"); 
  	  fSlotChanRMS[i]->GetXaxis()->SetTitle("Slot Channel"); fSlotChanRMS[i]->GetYaxis()->SetTitle("Raw RMS"); 
  	  fSlotChanMean_pfx[i]->GetXaxis()->SetTitle("Slot Channel"); fSlotChanMean_pfx[i]->GetYaxis()->SetTitle("Profiled Mean"); 
  	  fSlotChanRMS_pfx[i]->GetXaxis()->SetTitle("Slot Channel"); fSlotChanRMS_pfx[i]->GetYaxis()->SetTitle("Profiled RMS"); 
  	  fSlotChanFFT[i]->GetXaxis()->SetTitle("Slot Channel"); fSlotChanFFT[i]->GetYaxis()->SetTitle("kHz");
    }
  
    // Persistent waveform for each fiber
    for(int i=0;i<120;i++) {
      // Takes too much RAM: fPersistentWav_by_Fiber.push_back(tfs->make<TH2I>(Form("Persistent_Waveform_Fiber#%d", i), Form("Persistent_Waveform_Fiber#%d", i), fNticks, 0, fNticks, 5000, 0, 5000));
      fPersistentFFT_by_Fiber.push_back(tfs->make<TH2F>(Form("Persistent_FFT_Fiber#%d", i), Form("Persistent_FFT_Fiber#%d", i), fNticks/2, 0, fNticks/2*fBinWidth, 150, -100, 50));
      fFFT_by_Fiber_pfx.push_back(tfs->make<TProfile>(Form("Profiled_FFT_Fiber#%d", i), Form("Profiled_FFT_Fiber#%d", i), fNticks/2, 0, fNticks/2*fBinWidth, -100, 50));
      //fPersistentWav_by_Fiber[i]->GetXaxis()->SetTitle("Tick"); fPersistentWav_by_Fiber[i]->GetYaxis()->SetTitle("ADC"); 
      fPersistentFFT_by_Fiber[i]->GetXaxis()->SetTitle("Frequency [kHz]"); fPersistentFFT_by_Fiber[i]->GetYaxis()->SetTitle("Amplitude [dB]"); 
      fFFT_by_Fiber_pfx[i]->GetXaxis()->SetTitle("Frequency [kHz]"); fFFT_by_Fiber_pfx[i]->GetYaxis()->SetTitle("Amplitude [dB]"); 
    }
  }

  //-----------------------------------------------------------------------
  void TpcMonitor::beginRun(const art::Run& run) {
    // place to read databases or run independent info
  }
  //-----------------------------------------------------------------------

  void TpcMonitor::reconfigure(fhicl::ParameterSet const& p){

    // reconfigure without recompiling
    // read in the parameters from the .fcl file
    // allows for interactive changes of the parameter values

    fTPCInput       = p.get< std::string >("TPCInputModule");
    fTPCInstance    = p.get< std::string >("TPCInstanceName");
    fRebinX         = p.get<int>("RebinFactorX");
    fRebinY         = p.get<int>("RebinFactorY");
    auto const *fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fNticks         = fDetProp->NumberTimeSamples();
    
    std::cout<< "Number of Ticks = "<< fNticks <<std::endl; 
    //get sampling rate
    fSampleRate = fDetProp->SamplingRate();
    std::cout<<"Sampling rate  = "<<fSampleRate<<std::endl;
    // width of frequencyBin in kHz
    fBinWidth = 1.0/(fNticks*fSampleRate*1.0e-6);
    std::cout<<"Bin Width (kHz)  = "<<fBinWidth<<std::endl;
    return;
  }

  //-----------------------------------------------------------------------

  void TpcMonitor::analyze(const art::Event& event) {
  	// Get channel map
  	art::ServiceHandle<dune::PdspChannelMapService> channelMap;
  	// TODO Use LOG_DEBUG
    LOG_INFO("TpcMonitor")
      << "-------------------- TPC TpcMonitor -------------------";

    // called once per event

    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    std::cout << "EventNumber = " << fEvent << std::endl;

    // Get the objects holding raw information: RawDigit for TPC data
    art::Handle< std::vector<raw::RawDigit> > RawTPC;
    event.getByLabel(fTPCInput, fTPCInstance, RawTPC);

    // Fill pointer vectors - more useful form for the raw data
    // a more usable form
    std::vector< art::Ptr<raw::RawDigit> > RawDigits;
    art::fill_ptr_vector(RawDigits, RawTPC);

    // Loop over all RawRCEDigits (entire channels)                                                                                                        
    for(auto const & dptr : RawDigits) {
      const raw::RawDigit & digit = *dptr;
      
      // Get the channel number for this digit
      uint32_t chan = digit.Channel();
      // number of samples in uncompressed ADC
      int nSamples = digit.Samples();
      unsigned int apa = std::floor( chan/fChansPerAPA );	  
      int pedestal = (int)digit.GetPedestal();
      
      std::vector<short> uncompressed(nSamples);
      // with pedestal	  
      raw::Uncompress(digit.ADCs(), uncompressed, pedestal, digit.Compression());

      // number of compressed ADCs
      nADC_comp=digit.NADC();
      // number of ADC uncompressed
      nADC_uncomp=uncompressed.size();	  
      // subtract pedestals
      std::vector<short> uncompPed(nSamples);
      
      for (int i=0; i<nSamples; i++) uncompPed.at(i)=uncompressed.at(i)-pedestal;
            
      // number of ADC uncompressed without pedestal
      nADC_uncompPed=uncompPed.size();	 
      
      // wavefrom histogram   
      int FiberID = channelMap->FiberIdFromOfflineChannel(chan);
      TH1D* histwav=new TH1D(Form("wav%d",(int)chan),Form("wav%d",(int)chan),nSamples,0,nSamples); 
      
      
      
      
      
      for(int k=0;k<(int)nADC_uncompPed;k++) {
	      histwav->SetBinContent(k+1, uncompPed.at(k));
	      // Fill Persistent waveform by fiber -- skip as it takes too much RAM
	      //fPersistentWav_by_Fiber.at(FiberID)->Fill(k+1, uncompPed.at(k));
	    }
	    
	    
	    
	    
	    // Do FFT for single waveforms
	    TH1D* histfft=new TH1D(Form("fft%d",(int)chan),Form("fft%d",(int)chan),nSamples,0,nSamples*fBinWidth); 
	    calculateFFT(histwav, histfft);
	    // Fill persistent/overlay FFT for each fiber/FEMB
	    for(int k=0;k<(int)nADC_uncompPed/2;k++) {
	      fPersistentFFT_by_Fiber.at(FiberID % 5)->Fill((k+0.5)*fBinWidth, histfft->GetBinContent(k+1));
	      fFFT_by_Fiber_pfx.at(FiberID % 5)->Fill((k+0.5)*fBinWidth, histfft->GetBinContent(k+1));
	    }
 
	    // Mean and RMS
	    float mean = meanADC(uncompPed);
	    float rms = rmsADC(uncompPed);
	     
      // U View, induction Plane	  
      if( fGeom->View(chan) == geo::kU){	
	      fRMSadc_U[apa]=0;
	      fMEANadc_U[apa]=0;
	      maxMean_U[apa] =0;
	      maxRMS_U[apa] =0;

	      fMEANadc_U[apa] = mean;
	      fRMSadc_U[apa] = rms;

	      if (maxMean_U[apa]<fMEANadc_U[apa]) {
	        maxMean_U[apa]=fMEANadc_U[apa];
	      }
        
	      if (maxRMS_U[apa]<fRMSadc_U[apa]) {
	        maxRMS_U[apa]=fRMSadc_U[apa];
	      }

	      fChanRMSU[apa]->Fill(chan,fRMSadc_U[apa]);
	      fChanMeanU[apa]->Fill(chan,fMEANadc_U[apa]);
	      fChanRMSU_pfx[apa]->Fill(chan, fRMSadc_U[apa], 1);
	      fChanMeanU_pfx[apa]->Fill(chan, fMEANadc_U[apa], 1);
	      
	      //fft
	      for(int l=0;l<nSamples/2;l++) { 
	        //for the 2D histos
	        fChanFFTU[apa]->Fill(chan, (l+0.5)*fBinWidth, histfft->GetBinContent(l+1));
	      }
      }// end of U View

      // V View, induction Plane
      if( fGeom->View(chan) == geo::kV){
        fRMSadc_V[apa]=0;
        fMEANadc_V[apa]=0;

	      maxMean_V[apa] =0;
        maxRMS_V[apa] =0;

        fMEANadc_V[apa] = mean;
        fRMSadc_V[apa] = rms;
        
        if (maxMean_V[apa]<fMEANadc_V[apa]) {
            maxMean_V[apa]=fMEANadc_V[apa];
        }

        if (maxRMS_V[apa]<fRMSadc_V[apa]) {
            maxRMS_V[apa]=fRMSadc_V[apa];
        }

        fChanRMSV[apa]->Fill(chan,fRMSadc_V[apa]);
        fChanMeanV[apa]->Fill(chan,fMEANadc_V[apa]);
        fChanRMSV_pfx[apa]->Fill(chan, fRMSadc_V[apa], 1);
	      fChanMeanV_pfx[apa]->Fill(chan, fMEANadc_V[apa], 1);
	      
	      //fft
	      for(int l=0;l<nSamples/2;l++) { 
	        //for the 2D histos
	        fChanFFTV[apa]->Fill(chan, (l+0.5)*fBinWidth, histfft->GetBinContent(l+1));
	      }
      }// end of V View               

      // Z View, collection Plane
      if( fGeom->View(chan) == geo::kZ){
        fRMSadc_Z[apa]=0;
	      fMEANadc_Z[apa]=0;
	      maxMean_Z[apa] =0;
        maxRMS_Z[apa] =0;
        fMEANadc_Z[apa] = mean;
        fRMSadc_Z[apa] = rms;
        if (maxMean_Z[apa]<fMEANadc_Z[apa]) {
            maxMean_Z[apa]=fMEANadc_Z[apa];
        }

        if (maxRMS_Z[apa]<fRMSadc_Z[apa]) {
            maxRMS_Z[apa]=fRMSadc_Z[apa];
        }

        fChanRMSZ[apa]->Fill(chan,fRMSadc_Z[apa]);
	      fChanMeanZ[apa]->Fill(chan,fMEANadc_Z[apa]);
	      fChanRMSZ_pfx[apa]->Fill(chan, fRMSadc_Z[apa], 1);
	      fChanMeanZ_pfx[apa]->Fill(chan, fMEANadc_Z[apa], 1);
	      
	      //fft
	      for(int l=0;l<nSamples/2;l++) {
	        //for the 2D histos
	        fChanFFTZ[apa]->Fill(chan, (l+0.5)*fBinWidth, histfft->GetBinContent(l+1));
	      }
      }// end of Z View
      
      // Mean/RMS by slot
      int SlotID = channelMap->SlotIdFromOfflineChannel(chan);
      int FiberNumber = channelMap->FEMBFromOfflineChannel(chan);
      int FiberChannelNumber = channelMap->FEMBChannelFromOfflineChannel(chan);
      uint32_t SlotChannelNumber = FiberNumber*128 + FiberChannelNumber; //128 channels per fiber
      fSlotChanMean.at(SlotID)->Fill(SlotChannelNumber, mean);
      fSlotChanRMS.at(SlotID)->Fill(SlotChannelNumber, rms);
      fSlotChanMean_pfx.at(SlotID)->Fill(SlotChannelNumber, mean, 1);
      fSlotChanRMS_pfx.at(SlotID)->Fill(SlotChannelNumber, rms, 1);
      
      // FFT by slot
      for(int l=0;l<nSamples/2;l++) {
	      //for the 2D histos
	      fSlotChanFFT.at(SlotID)->Fill(SlotChannelNumber, (l+0.5)*fBinWidth, histfft->GetBinContent(l+1));
	    }
      
      histwav->Delete(); 
      histfft->Delete();                                                                                                                    
      
    } // RawDigits
    
    return;
  }
  
  //-----------------------------------------------------------------------   
  // define RMS
  float TpcMonitor::rmsADC(std::vector< short > &uncomp)
  {
    int n = uncomp.size();
    float sum = 0.;
    for(int i = 0; i < n; i++){
      if(uncomp[i]!=0) sum += uncomp[i];
    }
    float mean = sum / n;
    sum = 0;
    for(int i = 0; i < n; i++)
      {
	   if (uncomp[i]!=0)     sum += (uncomp[i]-mean)*(uncomp[i]-mean);
      }
    return sqrt(sum / n);
  }

  //-----------------------------------------------------------------------  
  //define Mean
  float TpcMonitor::meanADC(std::vector< short > &uncomp)
  {
    int n = uncomp.size();
    float sum = 0.;
    for(int i = 0; i < n; i++)
      {
	if (uncomp[i]!=0) sum += abs(uncomp[i]);
      }
    return sum / n;
  }
  
  //-----------------------------------------------------------------------  
  //calculate FFT
  void TpcMonitor::calculateFFT(TH1D* hist_waveform, TH1D* hist_frequency) {
  
  int n_bins = hist_waveform->GetNbinsX();
  TH1* hist_transform = 0;

  // Create hist_transform from the input hist_waveform
  hist_transform = hist_waveform->FFT(hist_transform, "MAG");
  hist_transform -> Scale (1.0 / float(n_bins));
  int nFFT=hist_transform->GetNbinsX();
  
  double frequency;
  double amplitude;
  double amplitudeLog;
  
  // Loop on the hist_transform to fill the hist_transform_frequency                                                                                        
  for (int k = 0; k < nFFT/2; k++){

    frequency =  (k+0.5)*fBinWidth; // kHz
    amplitude = hist_transform->GetBinContent(k+1); 
    amplitudeLog = 20*log10(amplitude); // dB
    hist_frequency->Fill(frequency, amplitudeLog);
  }

  hist_transform->Delete();
  
}
  
 
  //-----------------------------------------------------------------------  
  void TpcMonitor::endJob() {

//    myfileU.close();
//    myfileV.close();
//    myfileZ.close();
    return;
  }
  
}

DEFINE_ART_MODULE(tpc_monitor::TpcMonitor)
  


#endif // TpcMonitore_module

