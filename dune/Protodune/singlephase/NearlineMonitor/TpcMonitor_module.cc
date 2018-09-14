///////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Class:       TpcMonitor_module
// Module type: analyzer
// File:        TpcMonitor_module.cc
// Author:      Jingbo Wang (jiowang@ucdavis.edu), February 2018.  Modifications by Tom Junk
//
// Modification: Maggie Greenwood July, 2018
//               Added large summary histograms.
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
#include "dune/Protodune/singlephase/RawDecoding/data/RDStatus.h"

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
#include "TProfile2D.h"

// C++ Includes
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
    int FEMBchanToHistogramMap(int, int);

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

    //int NumberOfRCEs; // unused

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

    // Mean and RMS distributions  by offline channel in each plane -- the vector indexes over APA number
    std::vector<TH1F*> fChanMeanDistU;
    std::vector<TH1F*> fChanRMSDistU;
    std::vector<TH1F*> fChanMeanDistV;
    std::vector<TH1F*> fChanRMSDistV;
    std::vector<TH1F*> fChanMeanDistZ;
    std::vector<TH1F*> fChanRMSDistZ;

    // stuck code fraction histograms by APA
    std::vector<TH1F*> fStuckCodeOffFrac;
    std::vector<TH1F*> fStuckCodeOnFrac;

    // stuck code fraction histograms by channel and plane
    std::vector<TProfile*> fChanStuckCodeOffFracU;
    std::vector<TProfile*> fChanStuckCodeOnFracU;
    std::vector<TProfile*> fChanStuckCodeOffFracV;
    std::vector<TProfile*> fChanStuckCodeOnFracV;
    std::vector<TProfile*> fChanStuckCodeOffFracZ;
    std::vector<TProfile*> fChanStuckCodeOnFracZ;

    // FFT by channel in each plane -- the vector indexes over APA number
    std::vector<TH2F*> fChanFFTU;
    std::vector<TH2F*> fChanFFTV;
    std::vector<TH2F*> fChanFFTZ;
    	
    // profiled Mean/RMS by offline channel
    std::vector<TProfile*> fChanRMSU_pfx;
    std::vector<TProfile*> fChanMeanU_pfx;
    std::vector<TProfile*> fChanRMSV_pfx;
    std::vector<TProfile*> fChanMeanV_pfx;
    std::vector<TProfile*> fChanRMSZ_pfx;
    std::vector<TProfile*> fChanMeanZ_pfx;

    // 2D histograms of all Mean/RMS by offline channel number
    // Intended as a color map with each bin to represent a single channel
    TProfile2D* fAllChanMean;
    TProfile2D* fAllChanRMS;

    //2Dhistograms of bits, using same mapping as the fAllChan histos
    //vector indexes over 0-11 bit numbers
    std::vector<TProfile2D*> fBitValue;

    
    // profiled over events Mean/RMS by slot number
    std::vector<TProfile*> fSlotChanMean_pfx;
    std::vector<TProfile*> fSlotChanRMS_pfx;	
    
    // FFT by slot number
    //std::vector<TH2F*> fSlotChanFFT;
    
    // Persistent and overlay wavefroms by fiber
    //std::vector<TH2F*> fPersistentFFT_by_Fiber;
    // change to only saving these for each APA
    std::vector<TH2F*> fPersistentFFT_by_APA;
    
    // Profiled fft by fiber
    std::vector<TProfile*> fFFT_by_Fiber_pfx;
    	
    // Rebin factor
    int fRebinX;
    int fRebinY; 

    // define nADC counts for uncompressed vs compressed
    unsigned int nADC_comp;
    unsigned int nADC_uncomp;
    unsigned int nADC_uncompPed;

    TH1F *fNTicksTPC;

    // Noise level cut parameters
    int fNoiseLevelMinNCountsU;
    int fNoiseLevelMinNCountsV;
    int fNoiseLevelMinNCountsZ;
    double fNoiseLevelNSigma;

    // Histograms to save dead/noisy channels
    TH1F* fNDeadChannelsHisto;
    TH1F* fNNoisyChannelsHistoFromNSigma;
    TH1F* fNNoisyChannelsHistoFromNCounts;
    TH1F* fNDeadChannelsList;
    TH1F* fNNoisyChannelsListFromNSigma;
    TH1F* fNNoisyChannelsListFromNCounts;

    TH1F* fNDeadChannelsHistoU;
    TH1F* fNNoisyChannelsHistoFromNSigmaU;
    TH1F* fNNoisyChannelsHistoFromNCountsU;

    TH1F* fNDeadChannelsHistoV;
    TH1F* fNNoisyChannelsHistoFromNSigmaV;
    TH1F* fNNoisyChannelsHistoFromNCountsV;

    TH1F* fNDeadChannelsHistoZ;
    TH1F* fNNoisyChannelsHistoFromNSigmaZ;
    TH1F* fNNoisyChannelsHistoFromNCountsZ;

    // define functions
    float rmsADC(std::vector< short > & uncompressed);
    float meanADC(std::vector< short > & uncompressed);
    void calculateFFT(TH1D* hist_waveform, TH1D* graph_frequency);
    void FillChannelHistos(TProfile* h1, double mean, double sigma, int& ndeadchannels, int& nnoisychannels_sigma, int& nnoisychannels_counts);
    geo::GeometryCore const * fGeom = &*(art::ServiceHandle<geo::Geometry>());

    std::vector<unsigned int> fApaLabelNum;


  }; // TpcMonitor

  //-----------------------------------------------------------------------

  TpcMonitor::TpcMonitor(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet), fRebinX(1), fRebinY(1), fApaLabelNum{3,5,2,6,1,4} {
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
    for(unsigned int i=0;i<fNofAPA;i++)
      {
	UChMin=fUChanMin + i*fChansPerAPA;
	UChMax=fUChanMax + i*fChansPerAPA;      
	VChMin=fVChanMin + i*fChansPerAPA;
	VChMax=fVChanMax + i*fChansPerAPA;      
	ZChMin=fZ0ChanMin + i*fChansPerAPA;
	//ZChMax=fZ0ChanMax + i*fChansPerAPA;
	ZChMax=fZ1ChanMax + i*fChansPerAPA; //including unused channels
      
	//std::cout<<"UCh:"<<UChMin<<" - "<<UChMax<<std::endl;
	//std::cout<<"VCh:"<<VChMin<<" - "<<VChMax<<std::endl;
	//std::cout<<"ZCh:"<<ZChMin<<" - "<<ZChMax<<std::endl;

	// summaries for all views

	unsigned int j=fApaLabelNum.at(i);

	fStuckCodeOffFrac.push_back(tfs->make<TH1F>(Form("fStuckCodeOffFrac%d",j),Form("Stuck-Off Code Fraction APA%d",j),100,0,1));
	fStuckCodeOnFrac.push_back(tfs->make<TH1F>(Form("fStuckCodeOnFrac%d",j),Form("Stuck-On Code Fraction APA%d",j),100,0,1));

	// U view
	fChanRMSU_pfx.push_back(tfs->make<TProfile>(Form("fChanRMSU%d_pfx", j),Form("Profiled raw-ped RMS vs Channel(Plane U, APA%d)", j),  UChMax - UChMin + 1, UChMin-0.5, UChMax+0.5, "s")); 
	fChanMeanU_pfx.push_back(tfs->make<TProfile>(Form("fChanMeanU%d_pfx",j),Form("Profiled raw-ped MEAN vs Channel(Plane U, APA%d)",j),  UChMax - UChMin + 1, UChMin-0.5, UChMax+0.5, "s")); 
	fChanFFTU.push_back(tfs->make<TH2F>(Form("fChanFFTU%d", j),Form("fChanFFT (Plane U, APA%d)", j),  UChMax - UChMin + 1, UChMin-0.5, UChMax+0.5, fNticks/2,0,fNticks/2*fBinWidth));
	fChanMeanDistU.push_back(tfs->make<TH1F>(Form("fChanMeanDistU%d",j),Form("Means of Channels in (Plane U, APA%d)",j), 4096, -0.5, 4095.5));
	fChanRMSDistU.push_back(tfs->make<TH1F>(Form("fChanRMSDistU%d",j),Form("RMSs of Channels in (Plane U, APA%d)",j), 100, 0, 50));
	fChanStuckCodeOffFracU.push_back(tfs->make<TProfile>(Form("fChanStuckCodeOffFracU%d",j),Form("Stuck-Off Code Fraction (Plane U, APA%d)",j), UChMax - UChMin + 1, UChMin-0.5, UChMax+0.5, "s"));
	fChanStuckCodeOnFracU.push_back(tfs->make<TProfile>(Form("fChanStuckCodeOnFracU%d",j),Form("Stuck-On Code Fraction (Plane U, APA%d)",j), UChMax - UChMin + 1, UChMin-0.5, UChMax+0.5, "s"));
      
	// V view
	fChanRMSV_pfx.push_back(tfs->make<TProfile>(Form("fChanRMSV%d_pfx",j),Form("Profiled raw-ped RMS vs Channel(Plane V, APA%d)",j),  VChMax - VChMin + 1, VChMin-0.5, VChMax+0.5, "s")); 
	fChanMeanV_pfx.push_back(tfs->make<TProfile>(Form("fChanMeanV%d_pfx",j),Form("Profiled raw-ped Mean vs Channel(Plane V, APA%d)",j),  VChMax - VChMin + 1, VChMin-0.5, VChMax+0.5, "s"));   
	fChanFFTV.push_back(tfs->make<TH2F>(Form("fChanFFTV%d", j),Form("fChanFFT (Plane V, APA%d)", j),  VChMax - VChMin + 1, VChMin-0.5, VChMax+0.5, fNticks/2,0,fNticks/2*fBinWidth));
	fChanMeanDistV.push_back(tfs->make<TH1F>(Form("fChanMeanDistV%d",j),Form("Means of Channels in (Plane V, APA%d)",j), 4096, -0.5, 4095.5));
	fChanRMSDistV.push_back(tfs->make<TH1F>(Form("fChanRMSDistV%d",j),Form("RMSs of Channels in (Plane V, APA%d)",j), 100, 0, 50));
	fChanStuckCodeOffFracV.push_back(tfs->make<TProfile>(Form("fChanStuckCodeOffFracV%d",j),Form("Stuck-Off Code Fraction (Plane V, APA%d)",j), VChMax - VChMin + 1, VChMin-0.5, VChMax+0.5, "s"));
	fChanStuckCodeOnFracV.push_back(tfs->make<TProfile>(Form("fChanStuckCodeOnFracV%d",j),Form("Stuck-On Code Fraction (Plane V, APA%d)",j), VChMax - VChMin + 1, VChMin-0.5, VChMax+0.5, "s"));
      
	// Z view                                                                                                                                                           
	fChanRMSZ_pfx.push_back(tfs->make<TProfile>(Form("fChanRMSZ%d_pfx",j),Form("Profiled raw-ped RMS vs Channel(Plane Z, APA%d)",j),  ZChMax - ZChMin + 1, ZChMin-0.5, ZChMax+0.5, "s")); 
	fChanMeanZ_pfx.push_back(tfs->make<TProfile>(Form("fChanMeanZ%d_pfx",j),Form("Profiled raw-ped Mean vs Channel(Plane Z, APA%d)",j),  ZChMax - ZChMin + 1, ZChMin-0.5, ZChMax+0.5, "s")); 
	fChanFFTZ.push_back(tfs->make<TH2F>(Form("fChanFFTZ%d", j),Form("fChanFFT (Plane Z, APA%d)", j),  ZChMax - ZChMin + 1, ZChMin-0.5, ZChMax+0.5, fNticks/2,0,fNticks/2*fBinWidth));
	fChanMeanDistZ.push_back(tfs->make<TH1F>(Form("fChanMeanDistZ%d",j),Form("Means of Channels in (Plane Z, APA%d)",j), 4096, -0.5, 4095.5));
	fChanRMSDistZ.push_back(tfs->make<TH1F>(Form("fChanRMSDistZ%d",j),Form("RMSs of Channels in (Plane Z, APA%d)",j), 100, 0, 50));
	fChanStuckCodeOffFracZ.push_back(tfs->make<TProfile>(Form("fChanStuckCodeOffFracZ%d",j),Form("Stuck-Off Code Fraction (Plane Z, APA%d)",j), ZChMax - ZChMin + 1, ZChMin-0.5, ZChMax+0.5, "s"));
	fChanStuckCodeOnFracZ.push_back(tfs->make<TProfile>(Form("fChanStuckCodeOnFracZ%d",j),Form("Stuck-On Code Fraction (Plane Z, APA%d)",j), ZChMax - ZChMin + 1, ZChMin-0.5, ZChMax+0.5, "s"));
      
	// Set titles
	fChanRMSU_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanRMSU_pfx[i]->GetYaxis()->SetTitle("raw RMS"); 
	fChanMeanU_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanMeanU_pfx[i]->GetYaxis()->SetTitle("raw Mean"); 
	fChanRMSV_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanRMSV_pfx[i]->GetYaxis()->SetTitle("raw RMS"); 
	fChanMeanV_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanMeanV_pfx[i]->GetYaxis()->SetTitle("raw Mean"); 
	fChanRMSZ_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanRMSZ_pfx[i]->GetYaxis()->SetTitle("raw RMS"); 
	fChanMeanZ_pfx[i]->GetXaxis()->SetTitle("Chan"); fChanMeanZ_pfx[i]->GetYaxis()->SetTitle("raw Mean"); 
	fChanFFTU[i]->GetXaxis()->SetTitle("Chan"); fChanFFTU[i]->GetYaxis()->SetTitle("kHz"); 
	fChanFFTV[i]->GetXaxis()->SetTitle("Chan"); fChanFFTV[i]->GetYaxis()->SetTitle("kHz"); 
	fChanFFTZ[i]->GetXaxis()->SetTitle("Chan"); fChanFFTZ[i]->GetYaxis()->SetTitle("kHz"); 
	fChanStuckCodeOffFracU[i]->GetXaxis()->SetTitle("Chan"); fChanStuckCodeOffFracZ[i]->GetYaxis()->SetTitle("Fraction"); 
	fChanStuckCodeOnFracU[i]->GetXaxis()->SetTitle("Chan"); fChanStuckCodeOnFracZ[i]->GetYaxis()->SetTitle("Fraction"); 
	fChanStuckCodeOffFracV[i]->GetXaxis()->SetTitle("Chan"); fChanStuckCodeOffFracZ[i]->GetYaxis()->SetTitle("Fraction"); 
	fChanStuckCodeOnFracV[i]->GetXaxis()->SetTitle("Chan"); fChanStuckCodeOnFracZ[i]->GetYaxis()->SetTitle("Fraction"); 
	fChanStuckCodeOffFracZ[i]->GetXaxis()->SetTitle("Chan"); fChanStuckCodeOffFracZ[i]->GetYaxis()->SetTitle("Fraction"); 
	fChanStuckCodeOnFracZ[i]->GetXaxis()->SetTitle("Chan"); fChanStuckCodeOnFracZ[i]->GetYaxis()->SetTitle("Fraction"); 

	fChanMeanDistU[i]->GetXaxis()->SetTitle("Mean (ADC counts)");
	fChanRMSDistU[i]->GetXaxis()->SetTitle("RMS (ADC counts)");
	fChanMeanDistV[i]->GetXaxis()->SetTitle("Mean (ADC counts)");
	fChanRMSDistV[i]->GetXaxis()->SetTitle("RMS (ADC counts)");
	fChanMeanDistZ[i]->GetXaxis()->SetTitle("Mean (ADC counts)");
	fChanRMSDistZ[i]->GetXaxis()->SetTitle("RMS (ADC counts)");

	//  Rebin histograms
	//std::cout<<"RebinX = "<<fRebinX<<"  RebinY = "<<fRebinY<<std::endl;
	fChanFFTU[i]->Rebin2D(fRebinX, fRebinY);
	fChanFFTV[i]->Rebin2D(fRebinX, fRebinY);
	fChanFFTZ[i]->Rebin2D(fRebinX, fRebinY);
      }
    
    //All in one view
    //make the histograms
    fAllChanMean = tfs->make<TProfile2D>("fAllChanMean", "Means for all channels", 240, -0.5, 239.5, 64, -0.5, 63.5);
    fAllChanRMS = tfs->make<TProfile2D>("fAllChanRMS", "RMS for all channels", 240, -0.5, 239.5, 64, -0.5, 63.5);
    //set titles and bin labels
    fAllChanMean->GetXaxis()->SetTitle("APA Number (online)"); fAllChanMean->GetYaxis()->SetTitle("Plane"); fAllChanMean->GetZaxis()->SetTitle("Raw Mean");
    fAllChanRMS->GetXaxis()->SetTitle("APA Number (online)"); fAllChanRMS->GetYaxis()->SetTitle("Plane"); fAllChanRMS->GetZaxis()->SetTitle("Raw RMS");
    fAllChanMean->GetXaxis()->SetLabelSize(.075); fAllChanMean->GetYaxis()->SetLabelSize(.05);
    fAllChanRMS->GetXaxis()->SetLabelSize(.075); fAllChanRMS->GetYaxis()->SetLabelSize(.05);
    fAllChanMean->GetXaxis()->SetBinLabel(40, "3"); fAllChanMean->GetXaxis()->SetBinLabel(120, "2"); fAllChanMean->GetXaxis()->SetBinLabel(200, "1");
    fAllChanRMS->GetXaxis()->SetBinLabel(40, "3"); fAllChanRMS->GetXaxis()->SetBinLabel(120, "2"); fAllChanRMS->GetXaxis()->SetBinLabel(200, "1");
    fAllChanMean->GetYaxis()->SetBinLabel(5, "U"); fAllChanMean->GetYaxis()->SetBinLabel(15, "V"); fAllChanMean->GetYaxis()->SetBinLabel(26, "Z");
    fAllChanMean->GetYaxis()->SetBinLabel(37, "U"); fAllChanMean->GetYaxis()->SetBinLabel(47, "V"); fAllChanMean->GetYaxis()->SetBinLabel(58, "Z");
    fAllChanRMS->GetYaxis()->SetBinLabel(5, "U"); fAllChanRMS->GetYaxis()->SetBinLabel(15, "V"); fAllChanRMS->GetYaxis()->SetBinLabel(26, "Z");
    fAllChanRMS->GetYaxis()->SetBinLabel(37, "U"); fAllChanRMS->GetYaxis()->SetBinLabel(47, "V"); fAllChanRMS->GetYaxis()->SetBinLabel(58, "Z");

    for(int i=0;i<12;i++)
      {
	fBitValue.push_back(tfs->make<TProfile2D>(Form("fBitValue%d",i),Form("Values for bit %d",i),240,-0.5,239.5,64,-0.5,63.5,0,1));
	fBitValue[i]->SetStats(false);
	fBitValue[i]->GetXaxis()->SetTitle("APA Number (online)"); fBitValue[i]->GetYaxis()->SetTitle("Plane"); fBitValue[i]->GetZaxis()->SetTitle("Bit Fraction On");
	fBitValue[i]->GetXaxis()->SetLabelSize(.075); fBitValue[i]->GetYaxis()->SetLabelSize(.05);
	fBitValue[i]->GetXaxis()->SetBinLabel(40, "3"); fBitValue[i]->GetXaxis()->SetBinLabel(120, "2"); fBitValue[i]->GetXaxis()->SetBinLabel(200, "1");
	fBitValue[i]->GetYaxis()->SetBinLabel(5, "U"); fBitValue[i]->GetYaxis()->SetBinLabel(15, "V"); fBitValue[i]->GetYaxis()->SetBinLabel(26, "Z");
	fBitValue[i]->GetYaxis()->SetBinLabel(37, "U"); fBitValue[i]->GetYaxis()->SetBinLabel(47, "V"); fBitValue[i]->GetYaxis()->SetBinLabel(58, "Z");
      }

    // Mean/RMS by slot channel number for each slot
    for(int i=0;i<30;i++) {
      int apaloc = fApaLabelNum[i/5];
      int slotloc = i % 5;
      
      fSlotChanMean_pfx.push_back(tfs->make<TProfile>(Form("APA%d_Slot%d_Mean", apaloc, slotloc), Form("APA %d Slot%d Mean_vs_SlotChannel", apaloc, slotloc), 512, 0, 512, "s")); 
      fSlotChanRMS_pfx.push_back(tfs->make<TProfile>(Form("APA%d_Slot%d_RMS", apaloc, slotloc), Form("APA %d Slot %d  RMS_vs_SlotChannel", apaloc, slotloc), 512, 0, 512, "s")); 
      //fSlotChanFFT.push_back(tfs->make<TH2F>(Form("APA%d_Slot%d_FFT", apaloc, slotloc), Form("APA %d Slot %d FFT_vs_SlotChannel", apaloc, slotloc), 512, 0, 512, fNticks/2, 0, fNticks/2*fBinWidth));
  	  
      fSlotChanMean_pfx[i]->GetXaxis()->SetTitle("Slot Channel"); fSlotChanMean_pfx[i]->GetYaxis()->SetTitle("Profiled Mean"); 
      fSlotChanRMS_pfx[i]->GetXaxis()->SetTitle("Slot Channel"); fSlotChanRMS_pfx[i]->GetYaxis()->SetTitle("Profiled RMS"); 
      //fSlotChanFFT[i]->GetXaxis()->SetTitle("Slot Channel"); fSlotChanFFT[i]->GetYaxis()->SetTitle("kHz");
    }

    unsigned int fembmap_by_fiberID[120] =
      {
	320,315,310,305,319,314,309,304,318,313,308,303,317,312,307,302,316,311,306,301,505,510,515,520,504,509,514,519,503,508,513,518,502,507,512,517,501,506,511,516,220,215,210,205,219,
	214,209,204,218,213,208,203,217,212,207,202,216,211,206,201,605,610,615,620,604,609,614,619,603,608,613,618,602,607,612,617,601,606,611,616,120,115,110,105,119,114,109,104,118,113,
        108,103,117,112,107,102,116,111,106,101,405,410,415,420,404,409,414,419,403,408,413,418,402,407,412,417,401,406,411,416
      };

  
    // FFT's by fiber
    for(int i=0;i<120;i++) {
      unsigned int imb = fembmap_by_fiberID[i];
      //fPersistentFFT_by_Fiber.push_back(tfs->make<TH2F>(Form("Persistent_FFT_FEMB_%d", imb), Form("FFT FEMB%d WIB%d", imb, ( (i/4) % 5)+1), fNticks/2, 0, fNticks/2*fBinWidth, 150, -100, 50));
      //fPersistentFFT_by_Fiber[i]->GetXaxis()->SetTitle("Frequency [kHz]"); fPersistentFFT_by_Fiber[i]->GetYaxis()->SetTitle("Amplitude [dB]"); 
      // still keep the profiled FFT's by FEMB
      fFFT_by_Fiber_pfx.push_back(tfs->make<TProfile>(Form("Profiled_FFT_FEMB_%d", imb), Form("Profiled FFT FEMB_%d WIB%d", imb, ( (i/4) %5)+1), fNticks/2, 0, fNticks/2*fBinWidth, -100, 50));
      fFFT_by_Fiber_pfx[i]->GetXaxis()->SetTitle("Frequency [kHz]"); fFFT_by_Fiber_pfx[i]->GetYaxis()->SetTitle("Amplitude [dB]"); 
    }
    // persistent FFT now by APA
    for (int i=0;i<6;++i)
      {
	fPersistentFFT_by_APA.push_back(tfs->make<TH2F>(Form("Persistent_FFT_APA_%d", fApaLabelNum[i]), Form("FFT APA%d ", fApaLabelNum[i]), fNticks/2, 0, fNticks/2*fBinWidth, 150, -100, 50));
        fPersistentFFT_by_APA[i]->GetXaxis()->SetTitle("Frequency [kHz]"); 
	fPersistentFFT_by_APA[i]->GetYaxis()->SetTitle("Amplitude [dB]"); 
      }

    fNTicksTPC = tfs->make<TH1F>("NTicksTPC","NTicks in TPC Channels",100,0,20000);

    // Dead/noisy channels
    fNDeadChannelsHisto = tfs->make<TH1F>("fNDeadChannelsHisto","Number of dead channels",fNofAPA+1,0,fNofAPA+1);
    fNDeadChannelsHisto->GetYaxis()->SetTitle("Number of dead channels");
    fNNoisyChannelsHistoFromNSigma = tfs->make<TH1F>("fNNoisyChannelsHistoFromNSigma","Number of noisy channels",fNofAPA+1,0,fNofAPA+1);
    fNNoisyChannelsHistoFromNSigma->GetYaxis()->SetTitle("Number of noisy channels");
    fNNoisyChannelsHistoFromNCounts = tfs->make<TH1F>("fNNoisyChannelsHistoFromNCounts",Form("Number of noisy channels above counts %i-%i-%i (U-V-Z)", fNoiseLevelMinNCountsU, fNoiseLevelMinNCountsV, fNoiseLevelMinNCountsZ), fNofAPA+1,0,fNofAPA+1);
    fNNoisyChannelsHistoFromNCounts->GetYaxis()->SetTitle("Number of noisy channels");

    fNDeadChannelsHistoU = tfs->make<TH1F>("fNDeadChannelsHistoU","Number of dead channels (Plane U)",fNofAPA+1,0,fNofAPA+1);
    fNDeadChannelsHistoU->GetYaxis()->SetTitle("Number of dead channels (Plane U)");
    fNNoisyChannelsHistoFromNSigmaU = tfs->make<TH1F>("fNNoisyChannelsHistoFromNSigmaU","Number of noisy channels (Plane U)",fNofAPA+1,0,fNofAPA+1);
    fNNoisyChannelsHistoFromNSigmaU->GetYaxis()->SetTitle("Number of noisy channels (Plane U)");
    fNNoisyChannelsHistoFromNCountsU = tfs->make<TH1F>("fNNoisyChannelsHistoFromNCountsU",Form("Number of noisy channels above %i counts  (Plane U)", fNoiseLevelMinNCountsU), fNofAPA+1,0,fNofAPA+1);
    fNNoisyChannelsHistoFromNCountsU->GetYaxis()->SetTitle("Number of noisy channels (Plane U)");

    fNDeadChannelsHistoV = tfs->make<TH1F>("fNDeadChannelsHistoV","Number of dead channels (Plane V)",fNofAPA+1,0,fNofAPA+1);
    fNDeadChannelsHistoV->GetYaxis()->SetTitle("Number of dead channels (Plane V)");
    fNNoisyChannelsHistoFromNSigmaV = tfs->make<TH1F>("fNNoisyChannelsHistoFromNSigmaV","Number of noisy channels (Plane V)",fNofAPA+1,0,fNofAPA+1);
    fNNoisyChannelsHistoFromNSigmaV->GetYaxis()->SetTitle("Number of noisy channels (Plane V)");
    fNNoisyChannelsHistoFromNCountsV = tfs->make<TH1F>("fNNoisyChannelsHistoFromNCountsV",Form("Number of noisy channels above %i counts  (Plane V)", fNoiseLevelMinNCountsV), fNofAPA+1,0,fNofAPA+1);
    fNNoisyChannelsHistoFromNCountsV->GetYaxis()->SetTitle("Number of noisy channels (Plane V)");

    fNDeadChannelsHistoZ = tfs->make<TH1F>("fNDeadChannelsHistoZ","Number of dead channels (Plane Z)",fNofAPA+1,0,fNofAPA+1);
    fNDeadChannelsHistoZ->GetYaxis()->SetTitle("Number of dead channels (Plane Z)");
    fNNoisyChannelsHistoFromNSigmaZ = tfs->make<TH1F>("fNNoisyChannelsHistoFromNSigmaZ","Number of noisy channels (Plane Z)",fNofAPA+1,0,fNofAPA+1);
    fNNoisyChannelsHistoFromNSigmaZ->GetYaxis()->SetTitle("Number of noisy channels (Plane Z)");
    fNNoisyChannelsHistoFromNCountsZ = tfs->make<TH1F>("fNNoisyChannelsHistoFromNCountsZ",Form("Number of noisy channels above %i counts  (Plane Z)", fNoiseLevelMinNCountsZ), fNofAPA+1,0,fNofAPA+1);
    fNNoisyChannelsHistoFromNCountsZ->GetYaxis()->SetTitle("Number of noisy channels (Plane Z)");

    fNDeadChannelsList = tfs->make<TH1F>("fNDeadChannelsList","List of dead channels",fGeom->Nchannels()+1,fUChanMin,fGeom->Nchannels()+1);
    fNDeadChannelsList->GetXaxis()->SetTitle("Channel ID");
    fNNoisyChannelsListFromNSigma = tfs->make<TH1F>("fNNoisyChannelsListFromNSigma","List of noisy channels",fGeom->Nchannels()+1,fUChanMin,fGeom->Nchannels()+1);
    fNNoisyChannelsListFromNSigma->GetXaxis()->SetTitle("Channel ID");
    fNNoisyChannelsListFromNCounts = tfs->make<TH1F>("fNNoisyChannelsListFromNCounts",Form("Number of noisy channels above counts %i-%i-%i (U-V-Z)", fNoiseLevelMinNCountsU, fNoiseLevelMinNCountsV, fNoiseLevelMinNCountsZ),fGeom->Nchannels()+1,fUChanMin,fGeom->Nchannels()+1);
    fNNoisyChannelsListFromNCounts->GetXaxis()->SetTitle("Channel ID");

    for(unsigned int i=0;i<fNofAPA;i++){
      unsigned int j=fApaLabelNum.at(i);
      TString apastring = Form("APA %i", j);
      fNDeadChannelsHisto->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNNoisyChannelsHistoFromNSigma->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNNoisyChannelsHistoFromNCounts->GetXaxis()->SetBinLabel(j+1, apastring.Data());

      fNDeadChannelsHistoU->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNNoisyChannelsHistoFromNSigmaU->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNNoisyChannelsHistoFromNCountsU->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNDeadChannelsHistoV->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNNoisyChannelsHistoFromNSigmaV->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNNoisyChannelsHistoFromNCountsV->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNDeadChannelsHistoZ->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNNoisyChannelsHistoFromNSigmaZ->GetXaxis()->SetBinLabel(j+1, apastring.Data());
      fNNoisyChannelsHistoFromNCountsZ->GetXaxis()->SetBinLabel(j+1, apastring.Data());
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
    fNoiseLevelMinNCountsU = p.get<int>("NoiseLevelMinNCountsU");
    fNoiseLevelMinNCountsV = p.get<int>("NoiseLevelMinNCountsV");
    fNoiseLevelMinNCountsZ = p.get<int>("NoiseLevelMinNCountsZ");
    fNoiseLevelNSigma     = p.get<double>("NoiseLevelNSigma");
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

    // Get RDStatus handle

    art::Handle< std::vector<raw::RDStatus> > RDStatusHandle;
    event.getByLabel(fTPCInput, fTPCInstance, RDStatusHandle);

    // Fill pointer vectors - more useful form for the raw data
    // a more usable form
    std::vector< art::Ptr<raw::RawDigit> > RawDigits;
    art::fill_ptr_vector(RawDigits, RawTPC);

    //for the large all channel summary histograms these are key points for bin mapping
    //for each offline numbered apa, the left most bin should be at the x value:
    int xEdgeAPA[6] = {0,0,80,80,160,160}; //these numbers may be adjusted to horizontally space out the histogram
    //for each of the apas, the bottom most bin should be at the y value:
    int yEdgeAPA[2] = {0,32}; //these numbers may be adjusted to vertically space out the histograms

    // example of retrieving RDStatus word and flags

    for ( auto const& rdstatus : (*RDStatusHandle) )
      {
	if (rdstatus.GetCorruptDataDroppedFlag())
	  {
	    LOG_INFO("TpcMonitor_module: ") << "Corrupt Data Dropped Flag set in RDStatus";
	  }
	//std::cout << "RDStatus:  Corrupt Data dropped " << rdstatus.GetCorruptDataDroppedFlag() << std::endl; 
	//std::cout << "RDStatus:  Corrupt Data kept " << rdstatus.GetCorruptDataKeptFlag() << std::endl; 
	//std::cout << "RDStatus:  Status Word " << rdstatus.GetStatWord() << std::endl; 
      }

    // Loop over all RawRCEDigits (entire channels)                                                                                                        
    for(auto const & dptr : RawDigits) {
      const raw::RawDigit & digit = *dptr;
      
      // Get the channel number for this digit
      uint32_t chan = digit.Channel();
      // number of samples in uncompressed ADC
      int nSamples = digit.Samples();
      fNTicksTPC->Fill(nSamples);
      unsigned int apa = std::floor( chan/fChansPerAPA );	  
      //int pedestal = (int)digit.GetPedestal();
      int pedestal = 0;  

      std::vector<short> uncompressed(nSamples);
      // with pedestal	  
      raw::Uncompress(digit.ADCs(), uncompressed, pedestal, digit.Compression());

      // number of compressed ADCs
      nADC_comp=digit.NADC();
      // number of ADC uncompressed
      nADC_uncomp=uncompressed.size();	  
      // subtract pedestals
      std::vector<short> uncompPed(nSamples);

      int nstuckoff=0;
      int nstuckon=0;
      for (int i=0; i<nSamples; i++) 
	{ 
	  auto adc=uncompressed.at(i);
	  auto adcl6b = adc & 0x3F;
	  if (adcl6b == 0) ++nstuckoff;
	  if (adcl6b == 0x3F) ++nstuckon;
	  uncompPed.at(i) = adc - pedestal;
        }
      float fracstuckoff = ((float) nstuckoff)/((float) nSamples);
      float fracstuckon = ((float) nstuckon)/((float) nSamples);

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
	fPersistentFFT_by_APA.at(apa)->Fill((k+0.5)*fBinWidth, histfft->GetBinContent(k+1));    // offline apa number.  Plot labels are online
	fFFT_by_Fiber_pfx.at(FiberID % 120)->Fill((k+0.5)*fBinWidth, histfft->GetBinContent(k+1));
      }

      // summary stuck code fraction distributions by APA -- here the APA is the offline APA number.  The plot labels contain the mapping

      fStuckCodeOffFrac[apa]->Fill(fracstuckoff);
      fStuckCodeOnFrac[apa]->Fill(fracstuckon);
 
      // Mean and RMS
      float mean = meanADC(uncompPed);
      float rms = rmsADC(uncompPed);

      //get ready to fill the summary plots
      //get the channel's FEMB and WIB
      int WIB = channelMap->WIBFromOfflineChannel(chan); //0-4
      int FEMB = channelMap->FEMBFromOfflineChannel(chan); //1-4
      int FEMBchan = channelMap->FEMBChannelFromOfflineChannel(chan);
      int iFEMB = ((WIB*4)+(FEMB-1)); //index of the FEMB 0-19
      //Get the location of any FEMBchan in the hitogram
      //put as a function for clenliness.
      int xBin = ((FEMBchanToHistogramMap(FEMBchan,0))+(iFEMB*4)+xEdgeAPA[apa]); // (fembchan location on histogram) + shift from mobo + shift from apa
      int yBin = ((FEMBchanToHistogramMap(FEMBchan,1))+yEdgeAPA[(apa%2)]); //(fembchan location on histogram) + shift from apa 

      fAllChanMean->Fill(xBin,yBin,mean); //histogram the mean
      fAllChanRMS->Fill(xBin,yBin,rms); //histogram the rms

      for (int i=0; i<nSamples; i++) //histogram the 12 bits
	{ 
	  auto adc=uncompressed.at(i);
	  int bitstring = adc;
	  for(int mm=0;mm<12;mm++)
	    {
	      // get the bit value from the adc
	      int bit = (bitstring%2);
	      fBitValue[mm]->Fill(xBin,yBin,bit);
	      bitstring = (bitstring/2);
	    }
	}

	     
      // U View, induction Plane	  
      if( fGeom->View(chan) == geo::kU){	
	fChanMeanU_pfx[apa]->Fill(chan, mean, 1);
	fChanRMSU_pfx[apa]->Fill(chan, rms, 1);
	fChanMeanDistU[apa]->Fill(mean);
	fChanRMSDistU[apa]->Fill(rms);
	fChanStuckCodeOffFracU[apa]->Fill(chan,fracstuckoff,1);
	fChanStuckCodeOnFracU[apa]->Fill(chan,fracstuckon,1);
	      
	//fft
	for(int l=0;l<nSamples/2;l++) { 
	  //for the 2D histos
	  fChanFFTU[apa]->Fill(chan, (l+0.5)*fBinWidth, histfft->GetBinContent(l+1));
	}

      }// end of U View

      // V View, induction Plane
      if( fGeom->View(chan) == geo::kV){
        fChanRMSV_pfx[apa]->Fill(chan, rms, 1);
	fChanMeanV_pfx[apa]->Fill(chan, mean, 1);
	fChanMeanDistV[apa]->Fill(mean);
	fChanRMSDistV[apa]->Fill(rms);
	fChanStuckCodeOffFracV[apa]->Fill(chan,fracstuckoff,1);
	fChanStuckCodeOnFracV[apa]->Fill(chan,fracstuckon,1);
	      
	//fft
	for(int l=0;l<nSamples/2;l++) { 
	  //for the 2D histos
	  fChanFFTV[apa]->Fill(chan, (l+0.5)*fBinWidth, histfft->GetBinContent(l+1));
	}

      }// end of V View               

      // Z View, collection Plane
      if( fGeom->View(chan) == geo::kZ){
	fChanMeanZ_pfx[apa]->Fill(chan, mean, 1);
	fChanRMSZ_pfx[apa]->Fill(chan, rms, 1);
	fChanMeanDistZ[apa]->Fill(mean);
	fChanRMSDistZ[apa]->Fill(rms);
	fChanStuckCodeOffFracZ[apa]->Fill(chan,fracstuckoff,1);
	fChanStuckCodeOnFracZ[apa]->Fill(chan,fracstuckon,1);
	      
	//fft
	for(int l=0;l<nSamples/2;l++) {
	  //for the 2D histos
	  fChanFFTZ[apa]->Fill(chan, (l+0.5)*fBinWidth, histfft->GetBinContent(l+1));
	}

      }// end of Z View
      
      // Mean/RMS by slot
      int SlotID = channelMap->SlotIdFromOfflineChannel(chan);
      int FiberNumber = channelMap->FEMBFromOfflineChannel(chan) - 1;
      int FiberChannelNumber = channelMap->FEMBChannelFromOfflineChannel(chan);
      uint32_t SlotChannelNumber = FiberNumber*128 + FiberChannelNumber; //128 channels per fiber
      fSlotChanMean_pfx.at(SlotID)->Fill(SlotChannelNumber, mean, 1);
      fSlotChanRMS_pfx.at(SlotID)->Fill(SlotChannelNumber, rms, 1);
      
      // FFT by slot
      for(int l=0;l<nSamples/2;l++) {
	//for the 2D histos
	// fSlotChanFFT.at(SlotID)->Fill(SlotChannelNumber, (l+0.5)*fBinWidth, histfft->GetBinContent(l+1));
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
  // Fill dead/noisy channels tree
  void TpcMonitor::FillChannelHistos(TProfile* h1, double mean, double sigma, int& ndeadchannels, int& nnoisychannels_sigma, int& nnoisychannels_counts){

    double rms_threshold = mean + fNoiseLevelNSigma*sigma;
    TString htitle = h1->GetTitle();

    for(Int_t j=1; j <= h1->GetNbinsX(); j++){

      int fChannelID = h1->GetBinCenter(j);
      double fChannelValue = h1->GetBinContent(j);

      if(fChannelValue == 0){ // dead channel
        ndeadchannels++;
        fNDeadChannelsList->SetBinContent(fChannelID, 1.0);
      }
      else{
        if(fChannelValue > rms_threshold){ // noisy channel far away from mean
          nnoisychannels_sigma++;
          fNNoisyChannelsListFromNSigma->SetBinContent(fChannelID, 1.0);
        }
        if(htitle.Contains("Plane U"))
	  {
	    if (fChannelValue > fNoiseLevelMinNCountsU)
	      { // noisy U channel above count threshold
		nnoisychannels_counts++;
		fNNoisyChannelsListFromNCounts->SetBinContent(fChannelID, 1.0);
	      }
	  }
	else if(htitle.Contains("Plane V"))
	  {
	    if (fChannelValue > fNoiseLevelMinNCountsV)
	      { // noisy V channel above count threshold
		nnoisychannels_counts++;
		fNNoisyChannelsListFromNCounts->SetBinContent(fChannelID, 1.0);
	      }
	  }
	else if(htitle.Contains("Plane Z"))
	  {
	    if (fChannelValue > fNoiseLevelMinNCountsZ)
	      { // noisy Z channel above count threshold
		nnoisychannels_counts++;
		fNNoisyChannelsListFromNCounts->SetBinContent(fChannelID, 1.0);
	      }
	  }
	else{
	  mf::LogVerbatim("TpcMonitor::FillChannelHistos")
	    << " Unknown histogram title: " << htitle.Data() << std::endl;
	}
      }
    }

    return;
  }

  //----------------------------------------------------------------------
  //define the mapping of FEMBchans to the histogram.
  int TpcMonitor::FEMBchanToHistogramMap(int FEMBchan, int coord){
    //to see the reason for this channel mapping, check DocDB 4064 Table 5
    //for one FEMB, this dictates the coordinates on the histogram as a 4X32 block.
    int FEMBchanToHistogram[128][2] = { {0,0},{0,1},{0,2},{0,3},{0,4},//for U
                                        {0,10},{0,11},{0,12},{0,13},{0,14},//for V
                                        {0,20},{0,21},{0,22},{0,23},{0,24},{0,25},//for Z
                                        {0,5},{0,6},{0,7},{0,8},{0,9},//for U
                                        {0,15},{0,16},{0,17},{0,18},{0,19},//for V
                                        {0,26},{0,27},{0,28},{0,29},{0,30},{0,31},//for Z
                                        {1,20},{1,21},{1,22},{1,23},{1,24},{1,25},//for Z
                                        {1,10},{1,11},{1,12},{1,13},{1,14},//for V
                                        {1,0},{1,1},{1,2},{1,3},{1,4},//for U
                                        {1,26},{1,27},{1,28},{1,29},{1,30},{1,31},//for Z
                                        {1,15},{1,16},{1,17},{1,18},{1,19},//for V
                                        {1,5},{1,6},{1,7},{1,8},{1,9},//for U
                                        {2,0},{2,1},{2,2},{2,3},{2,4},//for U
                                        {2,10},{2,11},{2,12},{2,13},{2,14},//for V
                                        {2,20},{2,21},{2,22},{2,23},{2,24},{2,25},//for Z
                                        {2,5},{2,6},{2,7},{2,8},{2,9},//for U
                                        {2,15},{2,16},{2,17},{2,18},{2,19},//for V
                                        {2,26},{2,27},{2,28},{2,29},{2,30},{2,31},//for Z
                                        {3,20},{3,21},{3,22},{3,23},{3,24},{3,25},//for Z
                                        {3,10},{3,11},{3,12},{3,13},{3,14},//for V
                                        {3,0},{3,1},{3,2},{3,3},{3,4},//for U
                                        {3,26},{3,27},{3,28},{3,29},{3,30},{3,31},//for Z
                                        {3,15},{3,16},{3,17},{3,18},{3,19},//for V
                                        {3,5},{3,6},{3,7},{3,8},{3,9} };//for U
    return FEMBchanToHistogram[FEMBchan][coord];
  }

  //-----------------------------------------------------------------------  
  void TpcMonitor::endJob() {

    // Find dead/noisy channels. Do this separately for each APA and for each view.
    std::vector<double> fURMS_mean; std::vector<double> fURMS_sigma;
    std::vector<double> fVRMS_mean; std::vector<double> fVRMS_sigma;
    std::vector<double> fZRMS_mean; std::vector<double> fZRMS_sigma;
    for(unsigned int i = 0; i < fNofAPA; i++){
      // U plane
      TH1F* h1 = (TH1F*)fChanRMSDistU.at(i);
      fURMS_mean.push_back(h1->GetMean());
      fURMS_sigma.push_back(h1->GetRMS());
      // V plane
      TH1F* h2 = (TH1F*)fChanRMSDistV.at(i);
      fVRMS_mean.push_back(h2->GetMean());
      fVRMS_sigma.push_back(h2->GetRMS());
      // Z plane
      TH1F* h3 = (TH1F*)fChanRMSDistZ.at(i);
      fZRMS_mean.push_back(h3->GetMean());
      fZRMS_sigma.push_back(h3->GetRMS());
    }

    std::vector<int> fUdch_vec; std::vector<int> fUnch_vec; std::vector<int> fUcch_vec;
    std::vector<int> fVdch_vec; std::vector<int> fVnch_vec; std::vector<int> fVcch_vec;
    std::vector<int> fZdch_vec; std::vector<int> fZnch_vec; std::vector<int> fZcch_vec;

    for(unsigned int i = 0; i < fNofAPA; i++){
      int ndeadchannels = 0; int nnoisychannels = 0; int nnoisychannels_counts = 0;

      // U plane
      TProfile* h1 = (TProfile*)fChanRMSU_pfx.at(i);
      FillChannelHistos(h1, fURMS_mean.at(i), fURMS_sigma.at(i), ndeadchannels, nnoisychannels, nnoisychannels_counts);
      fUdch_vec.push_back(ndeadchannels);
      fUnch_vec.push_back(nnoisychannels);
      fUcch_vec.push_back(nnoisychannels_counts);

      // V plane
      ndeadchannels = 0; nnoisychannels = 0; nnoisychannels_counts = 0;
      TProfile* h2 = (TProfile*)fChanRMSV_pfx.at(i);
      FillChannelHistos(h2, fVRMS_mean.at(i), fVRMS_sigma.at(i), ndeadchannels, nnoisychannels, nnoisychannels_counts);
      fVdch_vec.push_back(ndeadchannels);
      fVnch_vec.push_back(nnoisychannels);
      fVcch_vec.push_back(nnoisychannels_counts);

      // Z plane
      ndeadchannels = 0; nnoisychannels = 0; nnoisychannels_counts = 0;
      TProfile* h3 = (TProfile*)fChanRMSZ_pfx.at(i);
      FillChannelHistos(h3, fZRMS_mean.at(i), fZRMS_sigma.at(i), ndeadchannels, nnoisychannels, nnoisychannels_counts);
      fZdch_vec.push_back(ndeadchannels);
      fZnch_vec.push_back(nnoisychannels);
      fZcch_vec.push_back(nnoisychannels_counts);
    }

    // Fill summary histograms
    for(unsigned int i = 0; i < fNofAPA; i++){
      unsigned int j=fApaLabelNum.at(i);
      int nch = fUdch_vec.at(i) + fVdch_vec.at(i) + fZdch_vec.at(i);
      fNDeadChannelsHisto->SetBinContent(j+1, nch);
      nch = fUnch_vec.at(i) + fVnch_vec.at(i) + fZnch_vec.at(i);
      fNNoisyChannelsHistoFromNSigma->SetBinContent(j+1, nch);
      nch = fUcch_vec.at(i) + fVcch_vec.at(i) + fZcch_vec.at(i);
      fNNoisyChannelsHistoFromNCounts->SetBinContent(j+1, nch);

      fNDeadChannelsHistoU->SetBinContent(j+1, fUdch_vec.at(i));
      fNDeadChannelsHistoV->SetBinContent(j+1, fVdch_vec.at(i));
      fNDeadChannelsHistoZ->SetBinContent(j+1, fZdch_vec.at(i));

      fNNoisyChannelsHistoFromNSigmaU->SetBinContent(j+1, fUnch_vec.at(i));
      fNNoisyChannelsHistoFromNSigmaV->SetBinContent(j+1, fVnch_vec.at(i));
      fNNoisyChannelsHistoFromNSigmaZ->SetBinContent(j+1, fZnch_vec.at(i));

      fNNoisyChannelsHistoFromNCountsU->SetBinContent(j+1, fUcch_vec.at(i));
      fNNoisyChannelsHistoFromNCountsV->SetBinContent(j+1, fVcch_vec.at(i));
      fNNoisyChannelsHistoFromNCountsZ->SetBinContent(j+1, fZcch_vec.at(i));
    }

    //    myfileU.close();
    //    myfileV.close();
    //    myfileZ.close();
    return;
  }
  
}

DEFINE_ART_MODULE(tpc_monitor::TpcMonitor)
  


#endif // TpcMonitore_module

