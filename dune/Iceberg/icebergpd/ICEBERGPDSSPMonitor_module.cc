////////////////////////////////////////////////////////////////////////
// Class:       ICEBERGPDSSPMonitor
// Plugin Type: analyzer (art v3_02_06)
// File:        ICEBERGPDSSPMonitor_module.cc
//
// Generated at Wed Oct  9 23:18:50 2019 by Biswaranjan Behera using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

// C++ Includes
#include <memory>
#include <map>
#include <vector>
#include <set>


namespace icebergpd {
  class ICEBERGPDSSPMonitor;
}


class icebergpd::ICEBERGPDSSPMonitor : public art::EDAnalyzer {
public:
  explicit ICEBERGPDSSPMonitor(fhicl::ParameterSet const& pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICEBERGPDSSPMonitor(ICEBERGPDSSPMonitor const&) = delete;
  ICEBERGPDSSPMonitor(ICEBERGPDSSPMonitor&&) = delete;
  ICEBERGPDSSPMonitor& operator=(ICEBERGPDSSPMonitor const&) = delete;
  ICEBERGPDSSPMonitor& operator=(ICEBERGPDSSPMonitor&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

private:

  // Declare member data here.
  // The parameters we'll read from the .fcl file.

  std::string fOpDetWaveformModuleLabel;       // Input tag for OpDetWaveform
  std::string fOpHitModuleLabel;               // Input tag for OpHit

  double fSampleFreq;                          // Sampling frequency in MHz 

  // Map to store how many waveforms are on one optical channel
  std::map< int, TH1D* > avgWaveforms;
  std::map< int, int   > countWaveform;
  std::map<size_t,TH2D*> persistent_waveform_;

  std::map< size_t, TH2D* > Waveforms;

  std::map< int, TH1F* > maxadchist;
  std::map< int, int > ADC;
  std::map< int, int > maxadc;
  std::map< int, int > threshold;
  //  std::map< int, int > PEdist;
  std::map<int, int>::iterator adcit;
  std::map<int, int>::iterator itr; 
  std::map<int, int>::iterator i; 
  //std::map<int, int>::iterator it;

  // Implementation of required member function here.
  art::ServiceHandle< art::TFileService > tfs;
  
  
  //  TH2F* smear = tfs->make<TH2F>("smear", ";S-arapuca; X-arapuca;", 100, 1800, 18000, 100, 1800,18000);
  //  TH2F* smear = new TH2F("smear", ";S-arapuca; X-arapuca;", 100, 1800, 18000, 100, 1800,18000);

  TTree * fADCTree;
  Float_t fmaxadc3;
  Float_t fmaxadc7;
  Float_t fdistpe;
  Float_t fpedestal;
  Float_t fthres3;
  Float_t fthres7;
  Float_t fpe3;
  Float_t fpe7;

  Float_t fadc3;
  Float_t fadc7;

  float fbaseline3;
  float fbaseline7;

  std::vector< Int_t >   fadcval7;
  std::vector< Int_t >   fadcval3;
};


icebergpd::ICEBERGPDSSPMonitor::ICEBERGPDSSPMonitor(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fOpDetWaveformModuleLabel = pset.get<std::string>("OpDetWaveformLabel");
  fOpHitModuleLabel = pset.get<std::string>("OpHitLabel");

  fADCTree = tfs->make<TTree>("ADCTree","ADCTree");
  fADCTree->Branch("Channel7",                     &fmaxadc7,   "Channel7/F");
  fADCTree->Branch("Channel3",                     &fmaxadc3,   "Channel3/F");
  fADCTree->Branch("distpe",                       &fdistpe,   "distpe/F");
  fADCTree->Branch("pedestal",                     &fpedestal,   "pedestal/F");
  fADCTree->Branch("thres3",                       &fthres3,   "thres3/F");
  fADCTree->Branch("thres7",                     &fthres7,   "thres7/F");
  fADCTree->Branch("PE7",                     &fpe7,   "PE7/F");
  fADCTree->Branch("PE3",                     &fpe3,   "PE3/F");
  fADCTree->Branch("ADC7",                     &fadc7,   "ADC7/F");
  fADCTree->Branch("ADC3",                     &fadc3,   "ADC3/F");
  fADCTree->Branch("Baseline3",                     &fbaseline3,   "Baseline3/F");
  fADCTree->Branch("Baseline7",                     &fbaseline7,   "Baseline7/F");
  fADCTree->Branch("adcval7",               &fadcval7);
  fADCTree->Branch("adcval3",               &fadcval3);
}

void icebergpd::ICEBERGPDSSPMonitor::analyze(art::Event const& evt)
{

  // Get Ophit from the event
  auto OpHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fOpHitModuleLabel);

  // Get OpDetWaveforms from the event
  auto OpDetWaveformHandle = evt.getHandle< std::vector< raw::OpDetWaveform > >(fOpDetWaveformModuleLabel);

  //  std::cout<< "Event #" << evt.id().event() <<"\t" << OpDetWaveformHandle->size() << std::endl;

  // Obtain parameters from DetectorClocksService
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  fSampleFreq = clockData.OpticalClock().Frequency();

  int n = 0; int ntick = 0; int three = 0; int seven = 0;// float adc = 0;
  long int adcmax = 0;  long int adcval = 0; long int adcval3 = 0, adcval7 = 0, sig_adcval3 = 0, sig_adcval7 = 0; 
  
  double sumthres7 = 0, pedist3 = 0, pedist7 = 0, thres = 0, thres3 = 0, thres7 = 0; 
  maxadc[3] = 0.0;
  maxadc[7] = 0.0;
  for (size_t i = 0; i < OpDetWaveformHandle->size(); i++)
    {
      // This was probably required to overcome the "const" problem 
      // with OpDetPulse::Waveform()
      art::Ptr< raw::OpDetWaveform > waveformPtr(OpDetWaveformHandle, i);
      //	raw::OpDetWaveform pulse = *waveformPtr;
      
      int channel = waveformPtr->ChannelNumber();
      //std::cout<< "channel # -------" << channel  << std::endl;    
      if (channel ==3) three++;
      if (channel ==7) seven++;
      // channel 3 is S-ARAPUCA and channel 7 is X-ARAPUCA
      if(channel != 3 && channel != 7) continue;
	   //if(channel != 3) continue;
      //      std::cout<< "Event #" << evt.id().event() << "\t ............"<<std::endl;
  
      
      /*
      // Create the TH1 if it doesn't exist
	if (avgWaveforms.find(channel) == avgWaveforms.end() ) {
	  TString histName = TString::Format("avgwaveform_channel_%03i", channel);
	  avgWaveforms[channel] =  tfs->make< TH1D >(histName, ";t (us);", waveformPtr->size(), 0, double(waveformPtr->size()) / fSampleFreq);
	  //  avgWaveforms[channel] =  tfs->make< TH1D >(histName, ";time sample;", waveformPtr->size(), 0, double(waveformPtr->size()));
	}
	// Add this waveform to this histogram
	for (size_t tick = 0; tick < waveformPtr->size(); tick++) {
	  avgWaveforms[channel]->Fill(double(tick)/fSampleFreq, waveformPtr->at(tick));
	  //avgWaveforms[channel]->Fill(waveformPtr->at(tick));
	}
	*/

	// Count number of waveforms on each channel
	countWaveform[channel]++;
	

	if (Waveforms.find(channel) == Waveforms.end())
	  {
	    //TString histName = TString::Format("waveform_channel_%03i", channel);
	    Waveforms[channel] = tfs->make<TH2D>(Form("waveform_%d",channel),Form("waveform_%d",channel), 2000,0,2000, 20000, 0, 20000);
	    //pwave->SetTitle(Form("Persistent waveform - Channel %d",channel));
	    // pwave->GetYaxis()->SetTitle("ADC value");
	    // pwave->GetXaxis()->SetTitle("Time sample");
	    // persistent_waveform_[channel] = pwave;
	    //      std::cout << persistent_waveform_.find(channel)->second << std::endl;
	  }

	

	for (size_t tick = 0; tick < waveformPtr->size(); tick++) {
          //avgWaveforms[channel]->Fill(double(tick)/fSampleFreq, waveformPtr->at(tick));
          Waveforms[channel]->Fill(tick+1, waveformPtr->at(tick));
	  adcval =  waveformPtr->at(tick);
	  ADC[channel] = adcval;
	  adcmax = std::max(adcmax, adcval);

	  if (channel == 3) {
	    fadcval3    .emplace_back(waveformPtr->at(tick));	    
	    if (tick <= 700){thres += adcval; n++;  adcval3 = waveformPtr->at(tick) - 1587; 
	      thres3 += adcval3;
	    }
	    if (tick < 1300 && tick > 790){
	      sig_adcval3 = waveformPtr->at(tick) - 1587;
	      pedist3 += sig_adcval3;}}
	  //  if (channel == 3) {
	  //if (adcval < 1586) continue; 
	  //}

	  if (channel == 7) {
	    fadcval7    .emplace_back(waveformPtr->at(tick));	    
	    if (tick <= 700){sumthres7 += adcval; ntick++;  adcval7 = waveformPtr->at(tick) - 1524; thres7 += adcval7;}//if (adcval < 1586) continue;
	    if (tick < 1300 && tick > 790){
	      sig_adcval7 = waveformPtr->at(tick) -1524;
	      pedist7 += sig_adcval7;}}
	  //}
          //	    std::cout << "adcvalue \t"<< adcval <<"\t adcmax \t"<< adcmax <<std::endl; 
	 
	}
	
	//if (channel == 3) 	
	fthres3 = thres3;
	//if (channel == 7)
	fthres7 = thres7;
//	  std::cout <<"\t thres3 \t"<<  thres3 <<std::endl;
	fpe3 = pedist3;
	fpe7 = pedist7;
	fbaseline3 = thres/n;
	fbaseline7 = sumthres7/ntick;

	

	//	std::cout <<  "Event #" << evt.id().event() <<"\t" << n << "\t"<< thres <<"\t"<< thres/n <<std::endl;
	if (maxadchist.find(channel) == maxadchist.end()){
	  TString histname = TString::Format("Maxadc_channel_%03i", channel);
	  maxadchist[channel] = tfs->make<TH1F>(histname,";Maximum ADC; Events", 50, 0,20000);
	}

	maxadchist[channel]->Fill(adcmax);	

	for (adcit = ADC.begin(); adcit != ADC.end(); ++adcit) {
	  //std::cout << '\t' << adcit->first
	  //	    << '\t' << adcit->second << std::endl;
	  if (adcit->first == 3 )     fadc3 = adcit->second;
	  if (adcit->first == 7 )     fadc7 = adcit->second;
	}	 

	maxadc[channel] = adcmax;
	//threshold[channel] = thres;
	//PEdist[channel] = pedist; 

	//std::cout<<channel<<
	//	smear->Fill( maxadc(3).second, maxadc(7).second); 
	//


    
	//	std::cout<< "Event #........." << evt.id().event() <<"\t" << OpDetWaveformHandle->size() << "\t"<<channel<< std::endl;
    }
  //  fdistpe = pedist;
  //fpedestal = sumthres;
  //  std::cout << sumthres << "\t" << pedist <<std::endl;  
  for (itr = maxadc.begin(); itr != maxadc.end(); ++itr) { 
    //    std::cout << '\t' << itr->first 
    //	      << '\t' << itr->second << std::endl; 
    if (itr->first == 3 )     fmaxadc3 = itr->second;
    if (itr->first == 7 )     fmaxadc7 = itr->second;
  }   



  //for (i = threshold.begin(); i != threshold.end(); ++i) {
    //    std::cout << '\t' << i->first
    //		  << '\t' << i->second << std::endl;
    //if (i->first == 3 )     fthres3 = i->second;
    //if (i->first == 7 )     fthres7 = i->second;
  //}

  //  for (it = PEdist.begin(); it != PEdist.end(); ++it) {
    //     std::cout << '\t' << it->first
    //		  << '\t' << it->second << std::endl;
    //if (it->first == 3 )     fpe3 = it->second;
    //if (it->first == 7 )     fpe7 = it->second;
  //}



  fADCTree->Fill();
  fadcval3.clear();
  fadcval7.clear();
  //  std::cout << n << "\t" << three <<  "\t" << seven << "\t" <<"\t adcmax \t"<< adcmax << std::endl;		

	/* 
  if (persistent_waveform_.find(channel) == persistent_waveform_.end())
    {
      TH2D* pwave = tfs->make<TH2D>(Form("persistent_waveform_%d",channel),Form("persistent_waveform_%d",channel), 500,0, 2000, 20000, 0, 20000); 
      pwave->SetTitle(Form("Persistent waveform - Channel %d",channel));
      pwave->GetYaxis()->SetTitle("ADC value");
      pwave->GetXaxis()->SetTitle("Time sample");
      persistent_waveform_[channel] = pwave;
      //      std::cout << persistent_waveform_.find(channel)->second << std::endl;
    }

   TH2D *pwavep = persistent_waveform_[channel];

   for (size_t iadc=0; iadc < waveformPtr->size(); ++iadc)
    {

      auto adcval = waveformPtr->at(iadc);
      //adc_values_->Fill(adcval);
      //      std::cout << iadc <<"\t" << adcval << "\t" << nADC  << std::endl;
      ///> Save the waveforms for all traces like an oscilloscope
      pwavep->Fill(iadc,adcval,1);  // (x,y,weight=1)
      //      hist->SetBinContent(iadc+1,adcval);
    }

    }
*/

  
}


DEFINE_ART_MODULE(icebergpd::ICEBERGPDSSPMonitor)
