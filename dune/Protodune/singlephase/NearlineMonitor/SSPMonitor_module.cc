////////////////////////////////////////////////////////////////////////
//
//  PlotOpticalDetails_module.cc
//
//  Author: Tom Junk, based on code that was in the SSP raw decoder by Antonino Sergi
//
//  Module to provide basic photon detector information for protoDUNE
//  Nearline Monitoring
//
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TGraph.h"

// LArSoft includes
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

// C++ Includes
#include <memory>
#include <map>
#include <vector>
#include <set>

namespace nlana {
 
  class SSPMonitor : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    SSPMonitor(const fhicl::ParameterSet&);
    virtual ~SSPMonitor();
    void beginJob();
    void analyze (const art::Event&);     

  private:

    void calculateFFT(TH1D* hist_waveform, TH1D* graph_frequency);

    // The parameters we'll read from the .fcl file.
    std::string fOpDetWaveformModuleLabel;       // Input tag for OpDetWaveform
    std::string fOpHitModuleLabel;               // Input tag for OpHit
    double ADC_max;                            // axis boundaries for persistent_waveforms
    double ADC_min;
    int min_time;
    int max_time;

    int timesamples;

    double startTime;
    bool haveStartTime;

    // summary histograms

    TH1D *adc_values_;
    TH1D *peaks_;
    TH1D *areas_;
    TH1D *hit_times_;
    //TH1D *channels_;

    // histograms, one per channel.  Create them as we discover channels as we read in the data.

    std::map<size_t,TH1D*> chan_peaks_;  // peaks distribuion (in all hits)
    std::map<size_t,TH1D*> chan_areas_;  // hit areas
    std::map<size_t,TH2D*> persistent_waveform_;
    std::map<size_t,TH1D*> fft_;
    std::set<size_t> has_waveform;  // indexed by channel number

    TH1I *fHEventNumber;

  };

} 

//-----------------------------------------------------------------------
// Constructor
nlana::SSPMonitor::SSPMonitor(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
{

  // Get channel map
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;

  fOpDetWaveformModuleLabel = pset.get<std::string>("OpDetWaveformLabel");
  fOpHitModuleLabel = pset.get<std::string>("OpHitLabel");
  ADC_min=pset.get<int>("SSP_ADC_min");
  ADC_max=pset.get<int>("SSP_ADC_max");
  timesamples=pset.get<int>("SSP_TIMESAMPLES");
  min_time=pset.get<int>("SSP_min_time");
  max_time=pset.get<int>("SSP_max_time");

  haveStartTime = false;

  // summary histogram creation

  art::ServiceHandle<art::TFileService> tFileService;
  fHEventNumber = tFileService->make<TH1I>("EventNumber","SSP: EventNumber;Event Number",  100, 0, 10000);
  adc_values_ = tFileService->make<TH1D>("ssp_adc_values","SSP: ADC_Values;ADC Value",4096,-0.5,4095.5);  
  peaks_ = tFileService->make<TH1D>("peaks","Peak Amplitudes;Peak Amplitude",100,-30,50);
  areas_ = tFileService->make<TH1D>("areas","Hit Areas;Hit Area",100,-10000,15000);

  hit_times_ = tFileService->make<TH1D>("ssp_hit_times","Hit Times",1000,min_time,max_time);  
  hit_times_->GetYaxis()->SetTitle("Number of hits");
  hit_times_->GetXaxis()->SetTitle("Time [s]");
}

//-----------------------------------------------------------------------
// Destructor
nlana::SSPMonitor::~SSPMonitor() 
{}
 
//-----------------------------------------------------------------------
void nlana::SSPMonitor::beginJob()
{}

//-----------------------------------------------------------------------
void nlana::SSPMonitor::analyze(const art::Event& evt) 
{

  art::ServiceHandle<art::TFileService> tFileService;

  art::Handle< std::vector< recob::OpHit > > OpHitHandle;
  evt.getByLabel(fOpHitModuleLabel, OpHitHandle);

  art::Handle< std::vector< raw::OpDetWaveform > > OpDetWaveformHandle;
  evt.getByLabel(fOpDetWaveformModuleLabel, OpDetWaveformHandle);

  fHEventNumber->Fill(evt.event());

  for(unsigned int i = 0; i < OpHitHandle->size(); ++i)
    {
      art::Ptr< recob::OpHit > ohp(OpHitHandle, i);
      //recob::OpHit TheOpHit = *TheOpHitPtr;  -- see if we can use the ptr straight up.
      auto channel = ohp->OpChannel();
      if (chan_peaks_.find(channel) == chan_peaks_.end())
	{
	  chan_peaks_[channel] = tFileService->make<TH1D>(Form("peaks_%d",channel),Form("Peak Amplitudes_%d;Peak Amplitude",channel),100,-30,50);
	}
      if (chan_areas_.find(channel) == chan_areas_.end())
	{
	  chan_areas_[channel] = tFileService->make<TH1D>(Form("areas_%d",channel),Form("Peak Areas_%d;Peak Area",channel),100,-10000,15000);
	}
      peaks_->Fill(ohp->Amplitude());
      chan_peaks_[channel]->Fill(ohp->Amplitude());
      areas_->Fill(ohp->Area());
      chan_areas_[channel]->Fill(ohp->Area());

      double time = ohp->PeakTimeAbs()*1E-6;
      if (!haveStartTime) 
	{
	  haveStartTime = true;
	  startTime = time;
	}
      hit_times_->Fill(time - startTime);

    } // End loop over OpHits

  for(size_t iwaveform = 0; iwaveform < OpDetWaveformHandle->size(); ++iwaveform)
    {
      art::Ptr< raw::OpDetWaveform > odp(OpDetWaveformHandle, iwaveform);

      size_t nADC = odp->size();

      TH1D *hist = new TH1D("hist","hist",nADC,0,nADC);
      auto channel = odp->ChannelNumber();
      if (persistent_waveform_.find(channel) == persistent_waveform_.end())
	{
	  TH2D* pwave = tFileService->make<TH2D>(Form("persistent_waveform_%d",channel),Form("persistent_waveform_%d",channel), 500,0,timesamples, (int)(ADC_max-ADC_min),ADC_min,ADC_max);
          pwave->SetTitle(Form("Persistent waveform - Channel %d",channel));
          pwave->GetYaxis()->SetTitle("ADC value");
          pwave->GetXaxis()->SetTitle("Time sample");
	  persistent_waveform_[channel] = pwave;
	}
      if (fft_.find(channel) == fft_.end())
	{
	  TH1D *fftp = tFileService->make<TH1D>(Form("fft_channel_%d",channel),Form("fft_channel_%d",channel), 100,0,4);
          fftp->SetTitle(Form("FFT - Channel %d",channel));
          fftp->GetXaxis()->SetTitle("Frequency [MHz]");	  
	  fft_[channel] = fftp;
	}

      TH2D *pwavep = persistent_waveform_[channel];
      for (size_t iadc=0; iadc < nADC; ++iadc)
	{
	  auto adcval = odp->at(iadc);
	  adc_values_->Fill(adcval);

	  ///> Save the waveforms for all traces like an oscilloscope
	  pwavep->Fill(iadc,adcval,1);  // (x,y,weight=1)
	  hist->SetBinContent(iadc+1,adcval);
	}

      // save one waveform per channel

      if ( has_waveform.find(channel) == has_waveform.end() )
	{
	  has_waveform.insert(channel);
	  char histname[100];
	  sprintf(histname,"evt%i_channel%d",evt.event(), channel);
	  TH1D *htf = tFileService->make<TH1D>(histname,histname,nADC,0,nADC);
	  htf->Reset();
	  htf->Add(hist);
	}
	  
      // FFT on the single waveform, output divided by channel
      calculateFFT(hist, fft_[channel]);
      
      hist->Delete();

    } // End loop over OpDetWaveforms

}

void nlana::SSPMonitor::calculateFFT(TH1D* hist_waveform, TH1D* hist_frequency) {
  
  int n_bins = hist_waveform->GetNbinsX();
  TH1* hist_transform = 0;

  // Create hist_transform from the input hist_waveform
  hist_transform = hist_waveform->FFT(hist_transform, "MAG");
  hist_transform -> Scale (1.0 / float(n_bins));
  int nFFT=hist_transform->GetNbinsX();
  
  Double_t frequency;
  Double_t amplitude;
  
  // Loop on the hist_transform to fill the hist_transform_frequency                                                                                        
  for (int k = 2; k <= nFFT/40; ++k){

    frequency =  (k-1)/(n_bins/150.); // MHz
    amplitude = hist_transform->GetBinContent(k);

    hist_frequency->Fill(frequency, amplitude);
  }

  hist_transform->Delete();
  
}

namespace nlana {
  DEFINE_ART_MODULE(SSPMonitor)
}

