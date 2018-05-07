#ifndef SIGMOIDFILTER_H
#define SIGMOIDFILTER_H

////////////////////////////////////////////////////////////////////////
//
// Sigmoidfilter_module
//
// k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "dune/RunHistory/DetPedestalDUNE.h"

// lbne-raw-data
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"

// C++
#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

// ROOT
#include <TTree.h>
#include <TMath.h>
#include "TF1.h"
#include "TH2.h"
#include "TH1.h"
#include "TVirtualFFT.h"
#include "TStyle.h"

namespace lbne {
  class Sigmoidfilter;
}

class lbne::Sigmoidfilter : public art::EDProducer {

public:

  explicit Sigmoidfilter(fhicl::ParameterSet const& pset);
  virtual ~Sigmoidfilter();

  virtual void produce(art::Event& evt) override;
  void reconfigure(fhicl::ParameterSet const& pset);
  void beginRun(art::Run & run) override;
  void beginJob() override;
  void endRun(art::Run & run) override;

private:

  std::string fDigitModuleLabel;
  std::string fDigitModuleInstance;
  bool fMakeTree;

  TF1* fColFilterFunc;  ///< Parameterized collection filter function.
  TF1* fIndUFilterFunc; ///< Parameterized induction filter function.
  TF1* fIndVFilterFunc; ///< Parameterized induction filter function.
  
  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<lbne::ChannelMapService> fChannelMap;

  TH2F* RawFFT;
  TH2F* FixFFT;
};

//-------------------------------------------------------------------
lbne::Sigmoidfilter::Sigmoidfilter(fhicl::ParameterSet const& pset) {
  std::string colFilt               = pset.get<std::string>("ColFilter");
  fColFilterFunc = new TF1("colFilter", colFilt.c_str());

  std::string indUFilt               = pset.get<std::string>("IndUFilter");
  fIndUFilterFunc = new TF1("indUFilter", indUFilt.c_str());

  std::string indVFilt               = pset.get<std::string>("IndVFilter");
  fIndVFilterFunc = new TF1("indVFilter", indVFilt.c_str());


  this->reconfigure(pset);
  produces<std::vector<raw::RawDigit> >();

}

//-------------------------------------------------------------------
lbne::Sigmoidfilter::~Sigmoidfilter(){
}

//-------------------------------------------------------------------
void lbne::Sigmoidfilter::reconfigure(fhicl::ParameterSet const& pset) {
  fDigitModuleLabel    = pset.get<std::string>("DigitModuleLabel");
  fDigitModuleInstance = pset.get<std::string>("DigitModuleInstance");
  fMakeTree            = pset.get<bool>("MakeTree");

  // Make the filter functions.
  //std::string colFilt               = pset.get<std::string>("ColFilter");
  std::vector<double> colFiltParams = pset.get<std::vector<double> >("ColFilterParams");
  //fColFilterFunc = new TF1("colFilter", colFilt.c_str());
  for(unsigned int i=0; i<colFiltParams.size(); ++i)
    fColFilterFunc->SetParameter(i, colFiltParams[i]);
    
  //std::string indUFilt               = pset.get<std::string>("IndUFilter");
  std::vector<double> indUFiltParams = pset.get<std::vector<double> >("IndUFilterParams");
  //fIndUFilterFunc = new TF1("indUFilter", indUFilt.c_str());
  for(unsigned int i=0; i<indUFiltParams.size(); ++i)
    fIndUFilterFunc->SetParameter(i, indUFiltParams[i]);

  //std::string indVFilt               = pset.get<std::string>("IndVFilter");
  std::vector<double> indVFiltParams = pset.get<std::vector<double> >("IndVFilterParams");
  //fIndVFilterFunc = new TF1("indVFilter", indVFilt.c_str());
  for(unsigned int i=0; i<indVFiltParams.size(); ++i)
    fIndVFilterFunc->SetParameter(i, indVFiltParams[i]);
}

//-------------------------------------------------------------------
void lbne::Sigmoidfilter::beginRun(art::Run & run) {
  return;
}
//-------------------------------------------------------------------
void lbne::Sigmoidfilter::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  if (fMakeTree) {
    RawFFT = tfs->make<TH2F>("RawFFT", "Raw FFT for all channels less than 1000 KHz; Channel Number; Frequency (KHz)"            , 2048, 0, 2048, 1000, 0, 1000);
    FixFFT = tfs->make<TH2F>("FixFFT", "FFT for all channels less than 1000 KHz after filtering; Channel Number; Frequency (KHz)", 2048, 0, 2048, 1000, 0, 1000);
  }
}

//-------------------------------------------------------------------
void lbne::Sigmoidfilter::endRun(art::Run & run) {
  return;
}
  
//-------------------------------------------------------------------
void lbne::Sigmoidfilter::produce(art::Event& evt) {

  art::ServiceHandle<geo::Geometry> geo;

  art::Handle<std::vector<raw::RawDigit> > rawDigitHandle;
  evt.getByLabel(fDigitModuleLabel, fDigitModuleInstance, rawDigitHandle);
  std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);

  std::vector<raw::RawDigit> filterRawDigitVector; // The new Vector of RawDigits
  
  // If I want to ignore any frequency ranges....
  std::vector< std::pair<int,int> > ZeroFreq;
  //ZeroFreq.push_back( std::make_pair(276 , 285 ) );
  //ZeroFreq.push_back( std::make_pair(558 , 568 ) );
  //ZeroFreq.push_back( std::make_pair(837 , 849 ) );
  //ZeroFreq.push_back( std::make_pair(1116, 1127) );
  //ZeroFreq.push_back( std::make_pair(4340, 5205) );

  for (size_t DigLoop=0; DigLoop < rawDigitVector.size(); ++DigLoop) {
    int Channel     = rawDigitVector[DigLoop].Channel();
    size_t NADC     = rawDigitVector[DigLoop].NADC();
    double Pedestal = rawDigitVector[DigLoop].GetPedestal();
    const geo::View_t view = geo->View(Channel);
        
    //std::cout << "Looking at rawDigitVector["<<DigLoop<<"] it was on channel " << rawDigitVector[DigLoop].Channel() << "("<<Channel<<") it is in View " << view
    //	      << ", NADC is " << rawDigitVector[DigLoop].NADC() << " ("<<NADC<<")"
    //	      << ", pedestal is " << rawDigitVector[DigLoop].GetPedestal() << " ("<<Pedestal<<")"
    //	      << std::endl;
    
    // Fill the RawDigit histogram for this histogram.
    TH1F hRawDigit("hRawDigit","",NADC,0,NADC/2);
    TH1F hRawFFT("hRawFFT"  ,"",NADC,0,NADC);
    for (size_t ADCs=0; ADCs < NADC; ++ADCs) {
      hRawDigit.SetBinContent( ADCs+1, rawDigitVector[DigLoop].ADC(ADCs)-Pedestal );
    }
    for (size_t ww=NADC; ww<NADC; ++ww)
      hRawFFT.SetBinContent( ww, 0 );
    // Make the FFT for this channel.
    hRawDigit.FFT( (&hRawFFT) ,"MAG");
    for (size_t bin = 0; bin < NADC; ++bin) {
      double BinVal = hRawFFT.GetBinContent(bin+1);
      double freq = 2000. * bin / (double)NADC;
      if (freq < 1000 && BinVal < 1e5 && fMakeTree) {
	RawFFT->Fill( (Channel-0.5) , freq, BinVal );
      }
    }
      
    // I want to do an inverse FFT, so need to convert the tranformed FFT into an array....
    //double Re[NADC], Im[NADC];
    std::unique_ptr<double[]> Re( new double[NADC]);
    std::unique_ptr<double[]> Im( new double[NADC]);
    TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
    fft->GetPointsComplex(Re.get(),Im.get());
    delete fft;
    
    // Set the noisy frequency range bins to an average value.
    for (size_t aa=0; aa<ZeroFreq.size(); ++aa) {
      for (int bb=ZeroFreq[aa].first; bb<ZeroFreq[aa].second; ++bb) {
	double ReMeanVal=0;
	double ImMeanVal=0;
	int Range = 50;
	for (int cc=0; cc<Range; ++cc) {
	  ReMeanVal += Re[ZeroFreq[aa].first-cc] + Re[ZeroFreq[aa].second+cc];
	  ImMeanVal += Im[ZeroFreq[aa].first-cc] + Im[ZeroFreq[aa].second+cc];
	}
	ReMeanVal = ReMeanVal / Range;
	Re[bb]    = Re[1500-bb] = ReMeanVal;
	ImMeanVal = ImMeanVal / Range;
	Im[bb]    = Im[1500-bb] = ImMeanVal;
      }
    }
    
    // Apply the filter...    
    for (size_t bin = 0; bin < NADC; ++bin) {
      double freq = 2000. * bin / NADC;
      if (view == geo::kU) { // U plane 
	Re[bin] = Re[bin]*fIndUFilterFunc->Eval(freq);
	Im[bin] = Im[bin]*fIndUFilterFunc->Eval(freq);
      } else if ( view == geo::kV) { // V plane
	Re[bin] = Re[bin]*fIndVFilterFunc->Eval(freq);
	Im[bin] = Im[bin]*fIndVFilterFunc->Eval(freq);
      } else if ( view == geo::kZ) { // Collection plane
	Re[bin] = Re[bin]*fColFilterFunc->Eval(freq);
	Im[bin] = Im[bin]*fColFilterFunc->Eval(freq);
      }
      
      double MagVal =  pow ( Re[bin]*Re[bin] + Im[bin]*Im[bin], 0.5);
      if (TMath::IsNaN(MagVal)) MagVal = 0;
	            
      // Now do the big histograms...
      if (freq < 1000 && MagVal < 1e5 && fMakeTree) {
	FixFFT -> Fill( (Channel-0.5) , freq, MagVal );
      }
    }

    // I have applied the filter so now transform back....
    int NBins = NADC;
    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &NBins, "C2R");
    fft_back->SetPointsComplex(Re.get(),Im.get());
    fft_back->Transform();
    TH1 *hb=0;
    hb = TH1::TransformHisto(fft_back, hb, "Re");
    delete fft_back;
    std::vector<short> NewADC;
    for (int BinNum=0; BinNum<NBins; ++BinNum) {
      short Val = rawDigitVector[DigLoop].GetPedestal() + hb->GetBinContent(BinNum+1) / NBins;
      NewADC.push_back( Val);
    }
    delete hb;
    raw::RawDigit theRawDigit( rawDigitVector[DigLoop].Channel(), NewADC.size(), NewADC );
    theRawDigit.SetPedestal( rawDigitVector[DigLoop].GetPedestal(), rawDigitVector[DigLoop].GetSigma());
    filterRawDigitVector.push_back( theRawDigit );
  } // RawDigit Loop
  
  // save filtered waveforms to event
  evt.put(std::make_unique<decltype(filterRawDigitVector)>(std::move(filterRawDigitVector)));
  return;
}

DEFINE_ART_MODULE(lbne::Sigmoidfilter)

#endif //SIGMOIDFILTER_H
