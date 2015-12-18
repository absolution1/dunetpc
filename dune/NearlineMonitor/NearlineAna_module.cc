////////////////////////////////////////////////////////////////////////
// Class:       NearlineAna
// Module Type: analyzer
// File:        NearlineAna_module.cc
//
// Generated at Wed Dec 16 08:25:59 2015 by Jonathan Davies using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h" 

//ROOT
#include "TH1I.h"
#include "TTimeStamp.h"
#include "TTree.h"

//larsoft
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Hit.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"


namespace nearline {
  class NearlineAna;
}

class nearline::NearlineAna : public art::EDAnalyzer {
public:
  explicit NearlineAna(fhicl::ParameterSet const & p);
  NearlineAna(NearlineAna const &) = delete;
  NearlineAna(NearlineAna &&) = delete;
  NearlineAna & operator = (NearlineAna const &) = delete;
  NearlineAna & operator = (NearlineAna &&) = delete;
  void reconfigure(fhicl::ParameterSet const & p);
  void printConfig();
  void analyze(art::Event const & e) override;
  void beginJob();
  void endJob();

  //NoiseSpectrum plots
  void makeNoiseSpectrumPlots();
  void fillNoiseSpectrumPlots(art::Event const & e);


private:

  bool fMakeNoiseSpectrumPlots;
  std::vector<unsigned int> fNoiseSpectrumChannels;
  std::vector<TH1I*> fVecNoiseSpectrumPlots;
  art::InputTag fRawDigitsTag;

  // Variables needed for the header info tree:
  TTree*       fHeader;
  unsigned int fRun;
  unsigned int fSubrun;
  int          fFirstEvent;
  int          fLastEvent;
  int          fNevents;
  unsigned int fStartYear;
  unsigned int fEndYear;
  unsigned int fStartMonth;
  unsigned int fEndMonth;
  unsigned int fStartDay;
  unsigned int fEndDay;
  double       fStartHour;
  double       fEndHour;
  unsigned long long int fStartTime;
  unsigned long long int fEndTime;

};

////////////////////////////////////////////////////////////////////////////////

nearline::NearlineAna::NearlineAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fHeader(0),
  fRun(0),
  fSubrun(0),
  fFirstEvent(1e9),
  fLastEvent(-1),
  fNevents(0),
  fStartYear(0),
  fEndYear(0),
  fStartMonth(0),
  fEndMonth(0),
  fStartDay(0),
  fEndDay(0),
  fStartHour(0.0),
  fEndHour(0.0),
  fStartTime(-1), // this is an unsigned int so it will default to a huge number
  fEndTime(0)
{

  {
    mf::LogInfo logInfo("NearlineAna::NearlineAna");
    logInfo << "NearlineAna" << "\n";
  }
  reconfigure(p);
  printConfig();

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::reconfigure(fhicl::ParameterSet const & p){

  mf::LogInfo logInfo("NearlineAna::reconfigure");
  logInfo << "reconfigure" << "\n";

  fMakeNoiseSpectrumPlots = p.get<bool>("MakeNoiseSpectrumPlots", true);
  fNoiseSpectrumChannels = p.get<std::vector<unsigned int>>("NoiseSpectrumChannels", {1,2,3,4});
  fRawDigitsTag = p.get<art::InputTag>("RawDigitsTag", "a:b:c");
}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::printConfig(){

  mf::LogInfo logInfo("NearlineAna::printConfig");

  logInfo << "fRawDigitsTag: " << fRawDigitsTag << "\n";
  logInfo << "fMakeNoiseSpectrumPlots: " << fMakeNoiseSpectrumPlots << "\n";
  if(fNoiseSpectrumChannels.size()){
    logInfo << "fNoiseSpectrumChannels:";
    for(size_t i=0;i<fNoiseSpectrumChannels.size();i++){
      logInfo << " " << fNoiseSpectrumChannels.at(i);
    }//fNoiseSpectrumChannels
    logInfo << "\n";
  }
}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::beginJob(){

  mf::LogInfo logInfo("NearlineAna::beginJob");
  logInfo << "beginJob" << "\n";

  // book the header object
  art::ServiceHandle<art::TFileService> tfs;
  fHeader = tfs->make<TTree>("Header","Subrun Information");
  
  fHeader->Branch("Run",&fRun);
  fHeader->Branch("Subrun",&fSubrun);
  fHeader->Branch("FirstEvent",&fFirstEvent);
  fHeader->Branch("LastEvent",&fLastEvent);
  fHeader->Branch("Nevents",&fNevents);
  fHeader->Branch("StartYear",&fStartYear);
  fHeader->Branch("StartMonth",&fStartMonth);
  fHeader->Branch("StartDay",&fStartDay);
  fHeader->Branch("StartHour",&fStartHour);
  fHeader->Branch("EndYear",&fEndYear);
  fHeader->Branch("EndMonth",&fEndMonth);
  fHeader->Branch("EndDay",&fEndDay);
  fHeader->Branch("EndHour",&fEndHour);

  if(fMakeNoiseSpectrumPlots) makeNoiseSpectrumPlots();

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::endJob(){

  //
  // Compute header info.
  //

  //
  // DISECTING the time from evt.time().value() into "human readable" format to display the date & time
  //
  unsigned int hour, minute, second;
  int          nano;

  // Get the time stamp.  art::Timestamp::value() returns a TimeValue_t which is a typedef to unsigned long long.
  // The conventional use is for the upper 32 bits to have the seconds since 1970 epoch and the lower 32 bits to be
  // the number of nanoseconds with the current second.
  //
  // NOTE: It seems that the above is NOT the convention for the 35t events. They only seem to use the lower 32 bits
  //       for the year/month/day/hour/second. For now, I have reversed the values of lup and llo to get the time right.
  //
  //       THESE VARIABLES WILL NEED TO BE SWITCHED BACK IF USING THE LOWER 32 BITS FOR NANOSECONDS IS IMPLEMENTED!!!


  const unsigned long int mask32 = 0xFFFFFFFFUL;

  // taking start time apart

  // unsigned long int lup = ( fStartTime >> 32 ) & mask32;
  // unsigned long int llo = fStartTime & mask32;
  unsigned long int llo = ( fStartTime >> 32 ) & mask32; // reversed value (see above comment)
  unsigned long int lup = fStartTime & mask32;           // reversed value (see above comment)
  TTimeStamp ts1(lup, (int)llo);
  ts1.GetDate(kTRUE,0,&fStartYear,&fStartMonth,&fStartDay);
  ts1.GetTime(kTRUE,0,&hour,&minute,&second);
  nano = ts1.GetNanoSec();
  double sec = ((double)second + (double)nano/1.0e9);
  fStartHour = (double)hour + (double)minute/60.0 + sec/3600.0;

  // taking end time apart
  // lup = ( fEndTime >> 32 ) & mask32;
  // llo = fEndTime & mask32;
  llo = ( fEndTime >> 32 ) & mask32; // reversed value (see above comment)
  lup = fEndTime & mask32;           // reversed value (see above comment)
  TTimeStamp ts2(lup, (int)llo);
  ts2.GetDate(kTRUE,0,&fEndYear,&fEndMonth,&fEndDay);
  ts2.GetTime(kTRUE,0,&hour,&minute,&second);
  nano = ts2.GetNanoSec();
  sec = ((double)second + (double)nano/1.0e9);
  fEndHour = (double)hour + (double)minute/60.0 + sec/3600.0;

  fHeader->Fill();

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::analyze(art::Event const & e)
{

  //
  // Extract event info for the header...
  //
  unsigned int           run    = e.run();
  unsigned int           subrun = e.subRun();
  unsigned int           event  = e.id().event();
  unsigned long long int time   = e.time().value();

  fNevents++;
  fRun    = run;
  fSubrun = subrun;

  // Don't assume first/last events are coorelated with start/end times...
  if(time < fStartTime)        fStartTime = time;
  if((int)event < fFirstEvent) fFirstEvent = event;
  if(time > fEndTime)          fEndTime = time;
  if((int)event > fLastEvent)  fLastEvent = event;



  //
  // Fill the desired plots.
  //
  if(fMakeNoiseSpectrumPlots) fillNoiseSpectrumPlots(e);

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::makeNoiseSpectrumPlots(){

  mf::LogInfo logInfo("NearlineAna::makeNoiseSpectrumPlots");
  logInfo << "fNoiseSpectrumChannels:" << "\n";

  art::ServiceHandle<art::TFileService> tfs;

  for(auto channel: fNoiseSpectrumChannels){
    int numBins = 100;
    int xmin = 0;
    int xmax = 2048;
    std::string hist_name = "hnoise_spectrum_chan_" + std::to_string(channel);
    std::string hist_title = "Average ADC Spectrum (Channel " + std::to_string(channel) + ")";
    TH1I* histTemp = tfs->make<TH1I>(hist_name.c_str(), hist_title.c_str(), numBins, xmin, xmax);
    fVecNoiseSpectrumPlots.push_back(histTemp);
    logInfo << "channel: " << channel << " hist_name: " << hist_name << "\n";
  }
  logInfo << "\n";

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::fillNoiseSpectrumPlots(art::Event const & e){

  mf::LogInfo logInfo("NearlineAna::fillNoiseSpectrumPlots");

  for(size_t index=0;index<fVecNoiseSpectrumPlots.size();index++){
    auto channel = fNoiseSpectrumChannels.at(index);
    auto hist = fVecNoiseSpectrumPlots.at(index);
    logInfo << "channel: " << channel << " hist_title: " << hist->GetTitle() << "\n";
  }


  art::Handle<std::vector<raw::RawDigit> > digitHandle;

  bool retVal = e.getByLabel(fRawDigitsTag, digitHandle);
  if(retVal!=true){
    mf::LogWarning("NearlineAna::fillNoiseSpectrumPlots") << "Getting RawDigits FAIL: " << fRawDigitsTag << std::endl;
    return;
  }
  
  try { digitHandle->size(); }
  catch(std::exception e) {
    mf::LogError("NearlineAna::fillNoiseSpectrumPlots") << "WARNING: Issue with digitHandle for RawDigits" << std::endl;
    return;
  }

  if(!digitHandle.isValid()){
    mf::LogError("NearlineAna::fillNoiseSpectrumPlots") << "Run: " << e.run()
                            << ", SubRun: " << e.subRun()
                            << ", Event: " << e.event()
                            << " is NOT VALID";
    throw cet::exception("RawDigit NOT VALID");
    return;
  }

  size_t numDigitChans = digitHandle->size();
//  logInfo << "Got RawDigit Handle - size: " << numDigitChans << "\n";

  //Loop through the vector of rawDigits and pick out channels that match our list of channels
  for(size_t rdIter=0;rdIter<numDigitChans;rdIter++){
    art::Ptr<raw::RawDigit> digitVec(digitHandle, rdIter);
    auto channel =  digitVec->Channel();    
    for(size_t index=0;index<fNoiseSpectrumChannels.size();index++){
      auto this_channel = fNoiseSpectrumChannels.at(index);

      //Only proceed if this is a channel of interest
      if(this_channel!=channel) continue; 

      auto numSamples = digitVec->Samples();
      auto compression = digitVec->Compression();

      //Only proceed if there is a non-zero number of samples
      if(numSamples==0) continue;
//      logInfo << "rdIter: " << rdIter << " channel: " << channel << " this_channel: " << this_channel << "\n";
//      logInfo << "numSamples: " << numSamples << " compression: " << compression << "\n";

      //We should uncompress the samples in case there is some form of compression applied
      raw::RawDigit::ADCvector_t ADCsUncompressed(numSamples);
      raw::Uncompress(digitVec->ADCs(), ADCsUncompressed, compression);

      //Calculate the channel ADC average
      double averageADC=0;
      for(unsigned int sample=0;sample<numSamples;sample++){
        averageADC+=ADCsUncompressed[sample];
       // logInfo << "ADC: " << ADCsUncompressed[sample]<< " cumulative ADC: " << averageADC << "\n";
      }//sample

      averageADC/=numSamples; 

//      logInfo << "averageADC: " << averageADC << "\n";

      TH1I* histTemp = fVecNoiseSpectrumPlots.at(index);
      histTemp->Fill(averageADC);

    }//index
  }//rdIter



//unsigned int  Nchannels () const 
}

////////////////////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(nearline::NearlineAna)
