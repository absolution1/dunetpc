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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h" 

//cpp
#include <sstream>
#include <fstream>

//ROOT
#include "TH1I.h"
#include "TTimeStamp.h"
#include "TTree.h"

//larsoft
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

//dunetpc
#include "dune/NearlineMonitor/NearlineVersion.h"

namespace nearline {
  class NearlineAna;
  void BuildTPCChannelMap(std::string channelMapFile, std::map<int,int>& channelMap);
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
  void beginJob() override;
  void endJob() override;

  size_t getRawDigits(art::Event const & e, art::Handle<std::vector<raw::RawDigit>> & digitHandle);
  size_t getHits(art::Event const & e, art::Handle<std::vector<recob::Hit>> & hitsHandle);

  //HitsPerEvent plots
  void makeHitsPerEventPlots();
  void fillHitsPerEventPlots(art::Event const & e);

  //PedestalPerEvent plots
  void makePedestalPerEventPlots();
  void fillPedestalPerEventPlots(art::Event const & e);
  void writePedestalPerEventSummaryFile(std::string fileName);
 
  //PedestalPerTick plots
  void makePedestalPerTickPlots();
  void fillPedestalPerTickPlots(art::Event const & e);
  void writePedestalPerTickSummaryFile(std::string fileName);

private:

  bool fVerboseOutput;

  bool fUseOnlineChannels;
  std::string fChannelMapFile;
  std::map<int,int> fChannelMap;

  bool fMakeHitsPerEventPlots;
  std::vector<unsigned int> fHitsPerEventChannels;
  std::vector<TH1I*> fVecHitsPerEventPlots;
  art::InputTag fHitsTag;

  bool fMakePedestalPerEventPlots;
  bool fWritePedestalPerEventFile;
  std::string fPedestalPerEventFileName;
  std::vector<unsigned int> fPedestalPerEventChannels;
  std::vector<TH1I*> fVecPedestalPerEventPlots;
  art::InputTag fRawDigitsTag;

  bool fMakePedestalPerTickPlots;
  bool fWritePedestalPerTickFile;
  std::string fPedestalPerTickFileName;
  std::vector<unsigned int> fPedestalPerTickChannels;
  std::vector<TH1I*> fVecPedestalPerTickPlots;
  

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

  //Histogram to store the 
  TH1I* fHistNearlineVersion;

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

  fVerboseOutput = p.get<bool>("VerboseOutput", false);

  fUseOnlineChannels = p.get<bool>("UseOnlineChannels", true);
  fChannelMapFile = p.get<std::string>("TPCChannelMapFile");

  if(fUseOnlineChannels) BuildTPCChannelMap(fChannelMapFile, fChannelMap);

  fMakeHitsPerEventPlots = p.get<bool>("MakeHitsPerEventPlots", true);
  fHitsPerEventChannels = p.get<std::vector<unsigned int>>("HitsPerEventChannels", {1,2,3,4});

  fHitsTag = p.get<art::InputTag>("HitsTag", "a:b:c");

  fMakePedestalPerEventPlots = p.get<bool>("MakePedestalPerEventPlots", true);
  fPedestalPerEventChannels = p.get<std::vector<unsigned int>>("PedestalPerEventChannels", {1,2,3,4});
  fWritePedestalPerEventFile = p.get<bool>("WritePedestalPerEventFile", false);
  fPedestalPerEventFileName = p.get<std::string>("PedestalPerEventFileName", "PedestalPerEvent.txt");

  fMakePedestalPerTickPlots = p.get<bool>("MakePedestalPerTickPlots", true);
  fPedestalPerTickChannels = p.get<std::vector<unsigned int>>("PedestalPerTickChannels", {1,2,3,4});
  fWritePedestalPerTickFile = p.get<bool>("WritePedestalPerTickFile", false);
  fPedestalPerTickFileName = p.get<std::string>("PedestalPerTickFileName", "PedestalPerTick.txt");

  fRawDigitsTag = p.get<art::InputTag>("RawDigitsTag", "a:b:c");
}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::printConfig(){

  mf::LogInfo logInfo("NearlineAna::printConfig");

  logInfo << "fVerboseOutput: " << (fVerboseOutput ? "true" : "false") << "\n";
  logInfo << "fRawDigitsTag: " << fRawDigitsTag << "\n";
  logInfo << "fHitsTag: " << fHitsTag << "\n";
  logInfo << "fUseOnlineChannels: " << (fUseOnlineChannels ? "true" : "false") << "\n";

  logInfo << "fMakeHitsPerEventPlots: " << (fMakeHitsPerEventPlots ? "true" : "false") << "\n";
  if(fHitsPerEventChannels.size()){
    logInfo << "fHitsPerEventChannels";
    if(fUseOnlineChannels) logInfo << " (online/offline): ";
    else logInfo << "(offline): ";
    for(size_t i=0;i<fHitsPerEventChannels.size();i++){
      logInfo << " " << fHitsPerEventChannels.at(i);
      if(fUseOnlineChannels) logInfo << "/" << fChannelMap.at(fHitsPerEventChannels.at(i));
    }//fHitsPerEventChannels
    logInfo << "\n";
  }

  logInfo << "fMakePedestalPerEventPlots: " << (fMakePedestalPerEventPlots ? "true" : "false") << "\n";
  if(fPedestalPerEventChannels.size()){
    logInfo << "fPedestalPerEventChannels";
    if(fUseOnlineChannels) logInfo << " (online/offline): ";
    else logInfo << "(offline): ";
    for(size_t i=0;i<fPedestalPerEventChannels.size();i++){
      logInfo << " " << fPedestalPerEventChannels.at(i);
      if(fUseOnlineChannels) logInfo << "/" << fChannelMap.at(fPedestalPerEventChannels.at(i));
    }//fPedestalPerEventChannels
    logInfo << "\n";
  }
  logInfo << "fWritePedestalPerEventFile: " << (fWritePedestalPerEventFile ? "true" : "false") << "\n";
  if(fWritePedestalPerEventFile) logInfo << "fPedestalPerEventFileName: " << fPedestalPerEventFileName << "\n";


  logInfo << "fMakePedestalPerTickPlots: " << (fMakePedestalPerTickPlots ? "true" : "false") << "\n";
  if(fPedestalPerTickChannels.size()){
    logInfo << "fPedestalPerTickChannels";
    if(fUseOnlineChannels) logInfo << " (online/offline): ";
    else logInfo << "(offline): ";
    for(size_t i=0;i<fPedestalPerTickChannels.size();i++){
      logInfo << " " << fPedestalPerTickChannels.at(i);
      if(fUseOnlineChannels) logInfo << "/" << fChannelMap.at(fPedestalPerTickChannels.at(i));
    }//fPedestalPerTickChannels
    logInfo << "\n";
  }
  logInfo << "fWritePedestalPerTickFile: " << (fWritePedestalPerTickFile ? "true" : "false") << "\n";
  if(fWritePedestalPerTickFile) logInfo << "fPedestalPerTickFileName: " << fPedestalPerTickFileName << "\n";


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

  // Set Nearline Version Number
  fHistNearlineVersion = tfs->make<TH1I>("hist_nearline_version", "hist_nearline_version", 2, 0, 2);
  fHistNearlineVersion->GetXaxis()->SetBinLabel(1,"NearlineMinorVersion");
  fHistNearlineVersion->GetXaxis()->SetBinLabel(2,"NearlineMajorVersion");
  fHistNearlineVersion->SetBinContent(1, NearlineMinorVersion);
  fHistNearlineVersion->SetBinContent(2, NearlineMajorVersion);
  
  if(fMakeHitsPerEventPlots) makeHitsPerEventPlots();
  if(fMakePedestalPerEventPlots) makePedestalPerEventPlots();
  if(fMakePedestalPerTickPlots) makePedestalPerTickPlots();

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

  //Write out pedestal plots summary to files
  if(fWritePedestalPerEventFile) writePedestalPerEventSummaryFile(fPedestalPerEventFileName);
  if(fWritePedestalPerTickFile) writePedestalPerTickSummaryFile(fPedestalPerTickFileName);


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
  if(time < fStartTime && e.time() != art::Timestamp::invalidTimestamp())        fStartTime = time;
  if((int)event < fFirstEvent) fFirstEvent = event;
  if(time > fEndTime)          fEndTime = time;
  if((int)event > fLastEvent)  fLastEvent = event;

  //
  // Fill the desired plots.
  //
  if(fMakeHitsPerEventPlots) fillHitsPerEventPlots(e);
  if(fMakePedestalPerEventPlots) fillPedestalPerEventPlots(e);
  if(fMakePedestalPerTickPlots) fillPedestalPerTickPlots(e);

}

////////////////////////////////////////////////////////////////////////////////

size_t nearline::NearlineAna::getRawDigits(art::Event const & e, art::Handle<std::vector<raw::RawDigit>> & digitHandle){

  bool retVal = e.getByLabel(fRawDigitsTag, digitHandle);
  if(retVal!=true){
    mf::LogWarning("NearlineAna::getRawDigits") << "Getting RawDigits FAIL: " << fRawDigitsTag << std::endl;
    return 0;
  }
  
  try { digitHandle->size(); }
  catch(std::exception e) {
    mf::LogError("NearlineAna::getRawDigits") << "WARNING: Issue with digitHandle for RawDigits" << std::endl;
    return 0;
  }

  if(!digitHandle.isValid()){
    mf::LogError("NearlineAna::getRawDigits") << "Run: " << e.run()
                            << ", SubRun: " << e.subRun()
                            << ", Event: " << e.event()
                            << " is NOT VALID";
    throw cet::exception("RawDigit NOT VALID");
    return 0;
  }

  return digitHandle->size();
}

////////////////////////////////////////////////////////////////////////////////

size_t nearline::NearlineAna::getHits(art::Event const & e, art::Handle<std::vector<recob::Hit>> & hitsHandle){

  bool retVal = e.getByLabel(fHitsTag, hitsHandle);
  if(retVal!=true){
    mf::LogWarning("NearlineAna::getHits") << "Getting Hits FAIL: " << fHitsTag << std::endl;
      return 0;
  }
  try { hitsHandle->size(); }
  catch(std::exception e) {
    mf::LogError("NearlineAna::getHits") << "WARNING: Issue with hitsHandle for Hits" << std::endl;
    return 0;
  }

  if(!hitsHandle.isValid()){
    mf::LogError("NearlineAna::getHits") << "Run: " << e.run()
                            << ", SubRun: " << e.subRun()
                            << ", Event: " << e.event()
                            << " is NOT VALID";
    throw cet::exception("Hit NOT VALID");
    return 0;
  }

  return hitsHandle->size();

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::makeHitsPerEventPlots(){

  mf::LogInfo logInfo("NearlineAna::makeHitsPerEventPlots");
  if(fVerboseOutput) logInfo << "fHitsPerEventChannels:" << "\n";

  art::ServiceHandle<art::TFileService> tfs;

  for(auto channel: fHitsPerEventChannels){
    int numBins = 100;
    int xmin = 0;
    int xmax = 3200;
    std::string hist_name = "hhits_per_event_chan_" + std::to_string(channel);
    std::string hist_title = "Hits Per Event - Channel " 
    + (fUseOnlineChannels ? 
       "(online/offline) " + std::to_string(channel) + "/" + std::to_string(fChannelMap.at(channel)) : 
       "(offline) " + std::to_string(channel));

    TH1I* histTemp = tfs->make<TH1I>(hist_name.c_str(), hist_title.c_str(), numBins, xmin, xmax);
    histTemp->GetXaxis()->SetTitle("Hits per Event");
    histTemp->GetYaxis()->SetTitle("Events");
    fVecHitsPerEventPlots.push_back(histTemp);
    if(fVerboseOutput) logInfo << "channel: " << channel << " hist_name: " << hist_name << " hist_title: " << hist_title << "\n";
  }
  if(fVerboseOutput) logInfo << "\n";

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::makePedestalPerEventPlots(){


  mf::LogInfo logInfo("NearlineAna::makePedestalPerEventPlots");
  if(fVerboseOutput) logInfo << "fPedestalPerEventChannels:" << "\n";

  art::ServiceHandle<art::TFileService> tfs;

  for(auto channel: fPedestalPerEventChannels){
    int numBins = 128;
    int xmin = 0;
    int xmax = 4096;
    std::string hist_name = "hped_per_event_chan_" + std::to_string(channel);
    std::string hist_title = "Average ADC Per Event - Channel " 
    + (fUseOnlineChannels ? 
       "(online/offline) " + std::to_string(channel) + "/" + std::to_string(fChannelMap.at(channel)) : 
       "(offline) " + std::to_string(channel));

    TH1I* histTemp = tfs->make<TH1I>(hist_name.c_str(), hist_title.c_str(), numBins, xmin, xmax);
    histTemp->GetXaxis()->SetTitle("ADC");
    histTemp->GetYaxis()->SetTitle("Events");
    fVecPedestalPerEventPlots.push_back(histTemp);
    if(fVerboseOutput) logInfo << "channel: " << channel << " hist_name: " << hist_name << " hist_title: " << hist_title << "\n";
  }
  if(fVerboseOutput) logInfo << "\n";

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::makePedestalPerTickPlots(){


  mf::LogInfo logInfo("NearlineAna::makePedestalPerTickPlots");
  if(fVerboseOutput) logInfo << "fPedestalPerTickChannels:" << "\n";

  art::ServiceHandle<art::TFileService> tfs;

  for(auto channel: fPedestalPerTickChannels){
    int numBins = 128;
    int xmin = 0;
    int xmax = 4096;
    std::string hist_name = "hped_per_tick_chan_" + std::to_string(channel);

    std::string hist_title = "ADC Per Tick - Channel " 
    + (fUseOnlineChannels ? 
       "(online/offline) " + std::to_string(channel) + "/" + std::to_string(fChannelMap.at(channel)) : 
       "(offline) " + std::to_string(channel));

    TH1I* histTemp = tfs->make<TH1I>(hist_name.c_str(), hist_title.c_str(), numBins, xmin, xmax);
    histTemp->GetXaxis()->SetTitle("ADC");
    histTemp->GetYaxis()->SetTitle("Events");
    fVecPedestalPerTickPlots.push_back(histTemp);
    if(fVerboseOutput) logInfo << "channel: " << channel << " hist_name: " << hist_name << " hist_title: " << hist_title << "\n";
  }
  if(fVerboseOutput) logInfo << "\n";

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::fillHitsPerEventPlots(art::Event const & e){


  mf::LogInfo logInfo("NearlineAna::fillHistPerEventPlots");
  art::Handle<std::vector<recob::Hit>> hitHandle;
  size_t numHits = getHits(e, hitHandle);

  //Count the number of hits on channels we are watching
  std::vector<unsigned int> numHitsPerChannel(fHitsPerEventChannels.size(), 0);

  if(numHits==0) return;

  for(size_t hitIter=0;hitIter<numHits;hitIter++){
    art::Ptr<recob::Hit> hitVec(hitHandle,hitIter);
    raw::ChannelID_t channel = hitVec->Channel();
    // geo::View_t      view = hitVec->View();
    // geo::SigType_t   signalType = hitVec->SignalType();
    // geo::WireID      wireID = hitVec->WireID();
    // logInfo << "Hit Channel: " << channel 
    //         << " View: " << view 
    //         << " SignalType: " << signalType 
    //         << " WireID: " << wireID << "\n";

    for(size_t index=0;index<fHitsPerEventChannels.size();index++){
      auto this_channel = (fUseOnlineChannels ? fChannelMap.at(fHitsPerEventChannels.at(index)) : fHitsPerEventChannels.at(index));
      if(this_channel!=channel) continue;
      numHitsPerChannel.at(index) += 1;
    }//channel index
  }//hitIter

  //Now fill histogram
  for(size_t index=0;index<fHitsPerEventChannels.size();index++){
    TH1I* histTemp = fVecHitsPerEventPlots.at(index);
    histTemp->Fill(numHitsPerChannel.at(index));
  }//channel index

  if(fVerboseOutput){
    logInfo << "Channel (Number of Hits): \n";
    for(size_t index=0;index<fHitsPerEventChannels.size();index++){
      logInfo << fHitsPerEventChannels.at(index) << " (" << numHitsPerChannel.at(index) << ") ";
    }
    logInfo << "\n";
  }


}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::fillPedestalPerEventPlots(art::Event const & e){

  mf::LogInfo logInfo("NearlineAna::fillPedestalPerEventPlots");

  // for(size_t index=0;index<fVecPedestalPerEventPlots.size();index++){
  //   auto channel = fPedestalPerEventChannels.at(index);
  //   auto offline_channel = (fUseOnlineChannels ? fChannelMap.at(channel) : channel );
  //   auto hist = fVecPedestalPerEventPlots.at(index);
  //   logInfo << "channel " << (fUseOnlineChannels ? "(online/offline): " + std::to_string(channel) + "/" + std::to_string(offline_channel) : "(offline): " + std::to_string(offline_channel) );
  //   logInfo << " hist_title: " << hist->GetTitle() << "\n";
  // }


  art::Handle<std::vector<raw::RawDigit> > digitHandle;
  size_t numDigitChans = getRawDigits(e, digitHandle);

  //Loop through the vector of rawDigits and pick out channels that match our list of channels

  for(size_t rdIter=0;rdIter<numDigitChans;rdIter++){
    art::Ptr<raw::RawDigit> digitVec(digitHandle, rdIter);
    auto channel =  digitVec->Channel();    
    for(size_t index=0;index<fPedestalPerEventChannels.size();index++){
      auto this_channel = (fUseOnlineChannels ? fChannelMap.at(fPedestalPerEventChannels.at(index)) : fPedestalPerEventChannels.at(index));

      //Only proceed if this is a channel of interest
      if(this_channel!=channel) continue; 

      //DEBUG 
      if(fUseOnlineChannels && fVerboseOutput) logInfo << "this_channel (online/offline): " << this_channel << " (" << fPedestalPerEventChannels.at(index) << "/" << fChannelMap.at(fPedestalPerEventChannels.at(index)) << ")\n";

      auto numSamples = digitVec->Samples();
      auto compression = digitVec->Compression();

      //Only proceed if there is a non-zero number of samples
      if(numSamples==0) continue;

      //We should uncompress the samples in case there is some form of compression applied
      raw::RawDigit::ADCvector_t ADCsUncompressed(numSamples);
      raw::Uncompress(digitVec->ADCs(), ADCsUncompressed, compression);

      //Calculate the channel ADC average
      double averageADC=0;
      for(unsigned int sample=0;sample<numSamples;sample++){
        averageADC+=ADCsUncompressed[sample];
      }//sample

      averageADC/=numSamples; 

      TH1I* histTemp = fVecPedestalPerEventPlots.at(index);
      histTemp->Fill(averageADC);

    }//index
  }//rdIter

}


////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::fillPedestalPerTickPlots(art::Event const & e){

  mf::LogInfo logInfo("NearlineAna::fillPedestalPerTickPlots");

  // for(size_t index=0;index<fVecPedestalPerTickPlots.size();index++){
  //   auto channel = fPedestalPerTickChannels.at(index);
  //   auto offline_channel = (fUseOnlineChannels ? fChannelMap.at(channel) : channel );
  //   auto hist = fVecPedestalPerTickPlots.at(index);
  //   logInfo << "channel " << (fUseOnlineChannels ? "(online/offline): " + std::to_string(channel) + "/" + std::to_string(offline_channel) : "(offline): " + std::to_string(offline_channel) );
  //   logInfo << " hist_title: " << hist->GetTitle() << "\n";
  // }


  art::Handle<std::vector<raw::RawDigit> > digitHandle;
  size_t numDigitChans = getRawDigits(e, digitHandle);

  //Loop through the vector of rawDigits and pick out channels that match our list of channels

  for(size_t rdIter=0;rdIter<numDigitChans;rdIter++){
    art::Ptr<raw::RawDigit> digitVec(digitHandle, rdIter);
    auto channel =  digitVec->Channel();    
    for(size_t index=0;index<fPedestalPerTickChannels.size();index++){
      auto this_channel = (fUseOnlineChannels ? fChannelMap.at(fPedestalPerTickChannels.at(index)) : fPedestalPerTickChannels.at(index));

      //Only proceed if this is a channel of interest
      if(this_channel!=channel) continue; 

      //DEBUG 
      if(fUseOnlineChannels && fVerboseOutput) logInfo << "this_channel (online/offline): " << this_channel << " (" << fPedestalPerTickChannels.at(index) << "/" << fChannelMap.at(fPedestalPerTickChannels.at(index)) << ")\n";

      auto numSamples = digitVec->Samples();
      auto compression = digitVec->Compression();

      //Only proceed if there is a non-zero number of samples
      if(numSamples==0) continue;

      //We should uncompress the samples in case there is some form of compression applied
      raw::RawDigit::ADCvector_t ADCsUncompressed(numSamples);
      raw::Uncompress(digitVec->ADCs(), ADCsUncompressed, compression);

      TH1I* histTemp = fVecPedestalPerTickPlots.at(index);

      for(unsigned int sample=0;sample<numSamples;sample++){
        histTemp->Fill(ADCsUncompressed[sample]);
      }//sample
    }//index
  }//rdIter

}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::writePedestalPerEventSummaryFile(std::string fileName){
  std::ostringstream my_ostream;
  my_ostream << "online_channel " << "offline_channel " << "pedestal_mean " << "pedestal_rms " << "\n";

  for(size_t index=0;index<fPedestalPerEventChannels.size();index++){
    auto online_channel = -1;
    auto offline_channel = fPedestalPerEventChannels.at(index);
    if(fUseOnlineChannels){
      online_channel = fPedestalPerEventChannels.at(index);
      offline_channel = fChannelMap.at(fPedestalPerEventChannels.at(index));
    } 
    TH1I* histTemp = fVecPedestalPerEventPlots.at(index);
    double mean = histTemp->GetMean();
    double rms = histTemp->GetRMS();
    my_ostream << online_channel << " " << offline_channel << " " << mean << " " << rms << "\n";
  }//index

  std::ofstream outFile(fileName);
  if(outFile.is_open()) outFile << my_ostream.str();
  else mf::LogWarning("writePedestalPerEventSummaryFile") << "FAILED to open file: " << fileName;
  outFile.close();

  //DEBUG
  mf::LogInfo("writePedestalPerEventSummaryFile") << my_ostream.str();
}

////////////////////////////////////////////////////////////////////////////////

void nearline::NearlineAna::writePedestalPerTickSummaryFile(std::string fileName){
  std::ostringstream my_ostream;
  my_ostream << "online_channel " << "offline_channel " << "pedestal_mean " << "pedestal_rms " << "\n";

  for(size_t index=0;index<fPedestalPerTickChannels.size();index++){
    auto online_channel = -1;
    auto offline_channel = fPedestalPerTickChannels.at(index);
    if(fUseOnlineChannels){
      online_channel = fPedestalPerTickChannels.at(index);
      offline_channel = fChannelMap.at(fPedestalPerTickChannels.at(index));
    } 
    TH1I* histTemp = fVecPedestalPerTickPlots.at(index);
    double mean = histTemp->GetMean();
    double rms = histTemp->GetRMS();
    my_ostream << online_channel << " " << offline_channel << " " << mean << " " << rms << "\n";
  }//index

  std::ofstream outFile(fileName);
  if(outFile.is_open()) outFile << my_ostream.str();
  else mf::LogWarning("writePedestalPerEventSummaryFile") << "FAILED to open file: " << fileName;
  outFile.close();

  //DEBUG
  mf::LogInfo("writePedestalPerTickSummaryFile") << my_ostream.str();

}

////////////////////////////////////////////////////////////////////////////////


void nearline::BuildTPCChannelMap(std::string channelMapFile, std::map<int,int>& channelMap) {

  /// Builds TPC channel map from the map txt file

  channelMap.clear();

  int onlineChannel;
  int offlineChannel;

  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(channelMapFile, fullname);

  if (fullname.empty())
    mf::LogWarning("DAQToOffline") << "Input TPC channel map file " << channelMapFile << " not found in FW_SEARCH_PATH.  Using online channel numbers!" << std::endl;

  else {
    mf::LogVerbatim("DAQToOffline") << "Build TPC Online->Offline channel Map from " << fullname;
    std::ifstream infile(fullname);
    while (infile.good()) {
      infile >> onlineChannel >> offlineChannel;
      channelMap.insert(std::make_pair(onlineChannel,offlineChannel));
      mf::LogVerbatim("DAQToOffline") << "   " << onlineChannel << " -> " << offlineChannel;
    }
    std::cout << "channelMap has size " << channelMap.size() << ". If this is 2048, then it's fine even if the above lines skipped a 'few' channels..." << std::endl;
  }

}

DEFINE_ART_MODULE(nearline::NearlineAna)
