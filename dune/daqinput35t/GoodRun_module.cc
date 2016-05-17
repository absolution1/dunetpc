////////////////////////////////////////////////////////////////////////
// Class:       GoodRun
// Module Type: Analyzer
// File:        GoodRun_module.cc
//
// Module written by Karl Warburton to compare the timings for RCE, SSP
//  and PTB data streams. 
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/RawData/ExternalTrigger.h"
#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/raw.h"
#include "lardata/RawData/OpDetWaveform.h"
#include "lardata/RecoBase/OpHit.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"
#include "TH2D.h"

const int MaxSamples  = 15000;
const int MaxChannels = 2048;

namespace DAQToOffline {
  class GoodRun;
}

class DAQToOffline::GoodRun : public art::EDAnalyzer {
public:
  explicit GoodRun(fhicl::ParameterSet const & p);
  virtual ~GoodRun();

  void analyze(art::Event const & evt) override;
  void beginJob();
  void beginRun(const art::Run& r);
  void endRun(const art::Run& r);
  void endJob();
  void printParameterSet();

private:

  // Information for good run list histogram
  TTree* fTree;
  int RunNumber;
  int TotEvents;
  int NumberOfRCEs;
  double FracRCEs;
  double AvADCPedDiff;
  double FracWaveforms;
  double FracOpHits;
  int nPTBTrigsOn110;
  int nPTBTrigsOn111;
  int nPTBTrigsOn112;
  int nPTBTrigsOn113;
  int nPTBTrigsOn114;
  int nPTBTrigsOn115;
  int EvJustSSPs;

  //TH2D* ADCDiff;
  //TH2D* ADCPed;

  /*
  TTree* FlatTree;
  int Event;
  int DigSize;
  int NSamples;
  float Channel[MaxChannels];
  float ADCs[MaxChannels][MaxSamples];
  float Pedestal[MaxChannels];
  */
  std::string fCounterModuleLabel, fWaveformModuleLabel, fRawDigitModuleLabel, fOpHitModuleLabel;
  
};

DAQToOffline::GoodRun::GoodRun(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
    //---------------------------RCE--------------------------------------
  , fCounterModuleLabel (pset.get<std::string>("CounterModuleLabel"))
  , fWaveformModuleLabel(pset.get<std::string>("WaveformModuleLabel"))
  , fRawDigitModuleLabel(pset.get<std::string>("RawDigitModuleLabel"))
  , fOpHitModuleLabel   (pset.get<std::string>("OpHitModuleLabel"))
    //---------------------------SSP--------------------------------------
//  ,fSSPFragType        ( pset.get<std::string>("SSPFragType"))
//  ,fSSPRawDataLabel    ( pset.get<std::string>("SSPRawDataLabel"))
    //---------------------------PTB--------------------------------------
//  ,fPTBFragType        ( pset.get<std::string>("PTBFragType"))
//  ,fPTBRawDataLabel    ( pset.get<std::string>("PTBRawDataLabel"))
    //--------------------------------------------------------------------
{
}

DAQToOffline::GoodRun::~GoodRun() {
}

void DAQToOffline::GoodRun::beginJob() {
  
  art::ServiceHandle<art::TFileService> tfs;
  ///*  
  fTree = tfs->make<TTree>("RunList","RunList Information");
  fTree->Branch("RunNumber"     ,&RunNumber     );
  fTree->Branch("TotEvents"     ,&TotEvents     );
  fTree->Branch("NumberOfRCEs"  ,&NumberOfRCEs  );
  fTree->Branch("FracRCEs"      ,&FracRCEs      );
  fTree->Branch("AvADCPedDiff"  ,&AvADCPedDiff  );
  fTree->Branch("FracWaveforms" ,&FracWaveforms );
  fTree->Branch("FracOpHits"    ,&FracOpHits    );
  fTree->Branch("nPTBTrigsOn110",&nPTBTrigsOn110);
  fTree->Branch("nPTBTrigsOn111",&nPTBTrigsOn111);
  fTree->Branch("nPTBTrigsOn112",&nPTBTrigsOn112);
  fTree->Branch("nPTBTrigsOn113",&nPTBTrigsOn113);
  fTree->Branch("nPTBTrigsOn114",&nPTBTrigsOn114);
  fTree->Branch("nPTBTrigsOn115",&nPTBTrigsOn115);
  fTree->Branch("EvJustSSPs"    ,&EvJustSSPs    );
  //*/

  //ADCDiff = tfs->make<TH2D>("ADCDiff","Difference in average ADC value and pedestal; Channel number; Difference (ADCs)", 2049,0,2048, 400,-50,50);
  //ADCPed  = tfs->make<TH2D>("ADCPed" ,"Comparison of Average ADC values, and pedestal; Average ADC; Pedestal", 1200, 400, 1000, 1200, 400, 1000);

  /*
  FlatTree = tfs->make<TTree>("FlatDigitTree","FlatDigitTree");
  FlatTree->Branch("Event"   ,&Event   ,"Event/I"   );
  FlatTree->Branch("NSamples",&NSamples,"NSamples/I");
  FlatTree->Branch("DigSize" ,&DigSize, "DigSize/I" );
  FlatTree->Branch("Channel" ,&Channel ,"Channel[DigSize]/F" );
  FlatTree->Branch("Pedestal",&Pedestal,"Pedestal[DigSize]/F");
  FlatTree->Branch("ADCs"    ,&ADCs    ,"ADCs[DigSize][15000]/F");
  */
}

void DAQToOffline::GoodRun::beginRun(const art::Run & r) {

  RunNumber = TotEvents = NumberOfRCEs = FracRCEs = AvADCPedDiff = FracWaveforms = FracOpHits = nPTBTrigsOn110 = nPTBTrigsOn111 = nPTBTrigsOn112 = nPTBTrigsOn113 = nPTBTrigsOn114 = nPTBTrigsOn115 = EvJustSSPs = 0;
  std::cout << "At the start of the run....going to reset all the variables. " << std::endl;

}

void DAQToOffline::GoodRun::analyze(art::Event const & evt)
{
  if (TotEvents==0) {
    RunNumber = evt.run();
    std::cout << "Got runNumber it is " << RunNumber << std::endl;
  }
  ++TotEvents;

  // get raw::RawDigits
  art::Handle< std::vector< raw::RawDigit> > externalRawDigitListHandle;
  std::vector< art::Ptr< raw::RawDigit> > digits;
  if (evt.getByLabel(fRawDigitModuleLabel, externalRawDigitListHandle) )
    art::fill_ptr_vector(digits,externalRawDigitListHandle);

  /*
  Event = evt.event();
  DigSize = digits.size();
  NSamples = digits[0]->Samples();
  for (int dig=0; dig<DigSize; ++dig ) {
    int Chan = digits[dig]->Channel();
    Channel[dig]  = Chan;
    Pedestal[dig] = digits[dig]->GetPedestal();
    for (int tick=0; tick<NSamples; ++tick) {
      ADCs[dig][tick] = digits[dig]->ADC(tick);
      //if (tick<5)std::cout << "Setting ADCs["<<dig<<"]["<<tick<<"] to " << ADCs[dig][tick] << " " << digits[dig]->ADC(tick) << std::endl;
    }
  }
  if (digits.size()) FlatTree->Fill();
  */

  // get raw::ExternalTriggers
  art::Handle< std::vector< raw::ExternalTrigger> > externalTriggerListHandle;
  std::vector< art::Ptr< raw::ExternalTrigger> > trigs;
  if (evt.getByLabel(fCounterModuleLabel, externalTriggerListHandle) )
    art::fill_ptr_vector(trigs,externalTriggerListHandle);

  // get raw::OpDetWaveforms
  art::Handle< std::vector< raw::OpDetWaveform> > externalWaveformListHandle;
  std::vector< art::Ptr< raw::OpDetWaveform> > waveforms;
  if (evt.getByLabel(fWaveformModuleLabel, externalWaveformListHandle) )
    art::fill_ptr_vector(waveforms,externalWaveformListHandle);

  // get recob::OpHits
  art::Handle< std::vector< recob::OpHit> > externalOpHitListHandle;
  std::vector< art::Ptr< recob::OpHit> > ophits;
  if (evt.getByLabel(fOpHitModuleLabel, externalOpHitListHandle) )
    art::fill_ptr_vector(ophits,externalOpHitListHandle);

  bool RCEsPresent = false;
  //bool PTBPresent  = false;
  bool WavePresent = false;
  //bool OHitPresent = false;

  // RCEs 
  if (digits.size()) {
    RCEsPresent = true;
    ++FracRCEs;
    if (NumberOfRCEs == 0) {
      NumberOfRCEs = digits.size() / 128;
      
      //GET THE LIST OF BAD CHANNELS.
      lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
      lariov::ChannelStatusProvider::ChannelSet_t const BadChannels = channelStatus.BadChannels();
      
      for (size_t dig=0; dig<digits.size(); ++dig) {
	double AvADC = 0;
	unsigned int Channel = digits[dig]->Channel();
	double Pedestal = digits[dig]->GetPedestal();
	for (size_t samp=0; samp<digits[dig]->Samples(); ++samp)
	  AvADC += digits[dig]->ADC(samp);
	AvADC = AvADC / digits[dig]->Samples();
	// Do I want to use this channel? Both AvADC and Pedestal, plus not a 'bad' channel...
	bool UseChan = true;
	if ( AvADC == 0 && Pedestal == 0)
	  UseChan = false;
	for (auto it = BadChannels.begin(); it != BadChannels.end(); it++) {
	  if(Channel==*it) {
	    UseChan = false;
	    break;
	  }
	} // Loop through bad chans.
	if ( UseChan ) {
	  double ADCPedDiff = AvADC - Pedestal;
	  //ADCDiff->Fill( Channel, ADCPedDiff );
	  //ADCPed ->Fill( AvADC  , Pedestal   );
	  AvADCPedDiff += fabs(ADCPedDiff);
	} // If have data for this tick, ie not turned off.
      } // Loop over digits
      AvADCPedDiff = AvADCPedDiff / digits.size();
    } // NumberOfRCEs == 0
  }
  
  // PTB
  if (trigs.size()) {
    //PTBPresent = true;
    for (size_t tr=0; tr<trigs.size(); ++tr) {
      if (trigs[tr]->GetTrigID() < 100) continue;
      if (trigs[tr]->GetTrigID() == 110) ++nPTBTrigsOn110;
      if (trigs[tr]->GetTrigID() == 111) ++nPTBTrigsOn111;
      if (trigs[tr]->GetTrigID() == 112) ++nPTBTrigsOn112;
      if (trigs[tr]->GetTrigID() == 113) ++nPTBTrigsOn113;
      if (trigs[tr]->GetTrigID() == 114) ++nPTBTrigsOn114;
      if (trigs[tr]->GetTrigID() == 115) ++nPTBTrigsOn115;
    }
  }

  // Waveforms
  if (waveforms.size()) {
    WavePresent = true;
    ++FracWaveforms;
    if (!RCEsPresent) ++EvJustSSPs;
  }

  // OpHits
  if (ophits.size()) {
    //OHitPresent = true;
    ++FracOpHits;
    if (!RCEsPresent && !WavePresent) ++EvJustSSPs;
  }
  
}

void DAQToOffline::GoodRun::endRun(const art::Run & r)
{
  FracRCEs      = FracRCEs / TotEvents;
  FracWaveforms = FracWaveforms / TotEvents;
  FracOpHits    = FracOpHits / TotEvents;
  std::cout << "At the end of the run, I am putting the following into the TTree:\n"
	    << RunNumber << " " << TotEvents << " " << NumberOfRCEs << " " << FracRCEs << " " << AvADCPedDiff << " " << FracWaveforms << " " << FracOpHits << " "
	    << nPTBTrigsOn110 << " " << nPTBTrigsOn111 << " " << nPTBTrigsOn112 << " " << nPTBTrigsOn113 << " " << nPTBTrigsOn114 << " " << nPTBTrigsOn115 << " " << EvJustSSPs 
	    << std::endl;
  fTree->Fill();
}

void DAQToOffline::GoodRun::endJob() {
}

DEFINE_ART_MODULE(DAQToOffline::GoodRun)
