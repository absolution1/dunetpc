////////////////////////////////////////////////////////////////////////
// Class:       Muoncounter
// Module Type: analyzer
// File:        Muoncounter_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
////////////////////////////////////////////////////////////////////////
#ifndef Muoncounter_Module
#define Muoncounter_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "dune/daqinput35t/PennToOffline.h"
#include "dune/daqinput35t/tpcFragmentToRawDigits.h"
#include "dune/daqinput35t/CheckTime.h"
#include "dune/daqinput35t/utilities/UnpackFragment.h"
#include "dune/daqinput35t/SSPReformatterAlgs.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TString.h"
#include "TText.h"
#include "TStyle.h"

#include "dune/NearlineMonitor/NearlineVersion.h"

namespace MyMuoncounter {

class Muoncounter : public art::EDAnalyzer {
public:
  explicit Muoncounter(fhicl::ParameterSet const & p);
  virtual ~Muoncounter();

    void beginJob();
    void endJob();

    void reconfigure(fhicl::ParameterSet const& pset);
    void analyze (const art::Event& evt); 

private:

  int l_TSU   = 0;   
  int u_TSU   = 43;  
  int l_BSU   = 44;  
  int u_BSU   = 92; 
  int l_Extra = 93;
  int u_Extra = 96;
  int l_Trig  = 100;
  int u_Trig  = 120;
  int event_Count = 0;
  double event_Length = 10e-3;

  TH1D *fHist1;
  TH1D *fHist2;
  TH1D *fHist3;

  double      fCombinedTimeDelay;

  // RCE Fragments
  std::string fRCEFragType, fRCERawDataLabel;
  // SSP Fragments
  std::string fSSPFragType, fSSPRawDataLabel;
  // PTB Fragments 
  std::string fPTBFragType, fPTBRawDataLabel;

  // Information for good run list histogram
  long long RCETime = 0, SSPTime = 0, PTBTime = 0;
  long long RCE_PTB_diff = 0, RCE_SSP_diff = 0, SSP_PTB_diff = 0;
  int nSynchronousEvents = 0;
  int nSSPPayloads = 0, nRCEPayloads = 0, nPTBPayloads = 0;
  int nConsistRCEPayloads = 0;
  int nPTBTrigsOn110 = 0, nPTBTrigsOn111 = 0, nPTBTrigsOn112 = 0, nPTBTrigsOn113 = 0, nPTBTrigsOn114 = 0, nPTBTrigsOn115 = 0;
  int sumNADCs = 0;
  TH1D *fGoodRunHisto;

  // Variables needed for the header info tree for Nearline:
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
  DAQToOffline::SSPReformatterAlgs sspReform;

  //Histogram to store the 
  TH1I* fHistNearlineVersion;

  //PTB Map
  std::string fPTBMapFile;
  std::string fPTBMapDir;
  
  std::map<int,int> fPTBMap;
};

Muoncounter::Muoncounter(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset),
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
    fEndTime(0),
    sspReform(pset.get<fhicl::ParameterSet>("SSPReformatter"))
{
   // Read in the parameters from the .fcl file.
    this->reconfigure(pset);
}
  
void Muoncounter::reconfigure(fhicl::ParameterSet const& p)
{   // Read parameters from the .fcl file. The names in the arguments
  // to p.get<TYPE> must match names in the .fcl file.
  fCombinedTimeDelay = p.get< double >("CombinedTimeDelay");
  fPTBMapFile        = p.get<std::string>("PTBMapFile");
  fPTBMapDir         = p.get<std::string>("PTBMapDir");
  fRCEFragType       = p.get<std::string>("RCEFragType");
  fRCERawDataLabel   = p.get<std::string>("RCERawDataLabel");
  fSSPFragType       = p.get<std::string>("SSPFragType");
  fSSPRawDataLabel   = p.get<std::string>("SSPRawDataLabel");
  fPTBFragType       = p.get<std::string>("PTBFragType");
  fPTBRawDataLabel   = p.get<std::string>("PTBRawDataLabel");
  DAQToOffline::BuildPTBChannelMap(fPTBMapDir, fPTBMapFile, fPTBMap);
  return;
}

Muoncounter::~Muoncounter()
{
  // Clean up dynamic memory and other resources here.
}

void Muoncounter::analyze(const art::Event& evt)
{

  // Get event number / timing information required by Nearline
  //
  // Extract event info for the header...
  //
  unsigned int           run    = evt.run();
  unsigned int           subrun = evt.subRun();
  unsigned int           event  = evt.id().event();
  unsigned long long int time   = evt.time().value();

  fNevents++;
  fRun    = run;
  fSubrun = subrun;

  // Don't assume first/last events are coorelated with start/end times...
  if(time < fStartTime && evt.time() != art::Timestamp::invalidTimestamp())        fStartTime = time;
  if((int)event < fFirstEvent) fFirstEvent = event;
  if(time > fEndTime)          fEndTime = time;
  if((int)event > fLastEvent)  fLastEvent = event;

  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

  // Reset some Good Events List parameters
  RCETime = SSPTime = PTBTime = 0;
  RCE_SSP_diff = RCE_PTB_diff = SSP_PTB_diff = 0;
  int NumADCs = 0, ConsistRCE = -1;
  int ThisEv110 = 0, ThisEv111 = 0, ThisEv112 = 0, ThisEv113 = 0, ThisEv114 = 0, ThisEv115 = 0;

  bool PTBPresent = true;
  art::Handle<artdaq::Fragments> PTBrawFragments;
  evt.getByLabel(fPTBRawDataLabel, fPTBFragType, PTBrawFragments);

  // Check if there is PTB data in this event
  // Don't crash code if not present, just don't save anything
  try { PTBrawFragments->size(); }
  catch(std::exception e) {
    mf::LogWarning("MuonCounter") << "WARNING: Raw PTB data not found in event " << evt.event();
    PTBPresent = false;
  }
  unsigned int total_Hits= 0;
  if (PTBPresent) {
    // Check that the data is valid
    if(!PTBrawFragments.isValid()){
      mf::LogError("MuonCounter") << "Run: " << evt.run() << ", SubRun: " << evt.subRun() << ", Event: " << evt.event() << " is NOT VALID";
      throw cet::exception("raw NOT VALID");
      return;
    }
    
    lbne::PennMicroSlice::Payload_Timestamp *FirstPTBTimestamp = nullptr;
    auto trigs = DAQToOffline::PennFragmentToExternalTrigger(*PTBrawFragments, fPTBMap, FirstPTBTimestamp);
    PTBTime = FirstPTBTimestamp->nova_timestamp;
    if (PTBTime) ++nPTBPayloads;
    //std::cout << "Got PTB start time, it is " << PTBTime << std::endl;
    
    total_Hits = trigs.size();
    
    for(unsigned int i = 0; i < total_Hits; i++) {
      int auxdetid = trigs.at(i).GetTrigID();
      if(auxdetid<=u_TSU) {
	fHist1->Fill(auxdetid);
      }	else if(auxdetid>=l_BSU && auxdetid<=u_BSU) {
	fHist2->Fill(auxdetid);
      }	else if(auxdetid>=l_Extra && auxdetid<=u_Extra) {
	fHist1->Fill(auxdetid-49);
      }	else {
	fHist3->Fill(auxdetid);
      }
      
      // Want to count the number of PTB trigs for in this event.
      if ( trigs.at(i).GetTrigID() == 110 ) ++ThisEv110;
      else if ( trigs.at(i).GetTrigID() == 111 ) ++ThisEv111;
      else if ( trigs.at(i).GetTrigID() == 112 ) ++ThisEv112;
      else if ( trigs.at(i).GetTrigID() == 113 ) ++ThisEv113;
      else if ( trigs.at(i).GetTrigID() == 114 ) ++ThisEv114;
      else if ( trigs.at(i).GetTrigID() == 115 ) ++ThisEv115;
      
      if ( trigs.at(i).GetTrigID() > 109 ) {
	std::cout << "Identifier:I Had a PTB trigger in event " << evt.event() << " on channel " << trigs.at(i).GetTrigID() << " at time " << trigs.at(i).GetTrigTime() << std::endl;
      }
    }
  }
  event_Count+=1;
  
  //---------------- SSP Good Event Timing stuff ----------------
  bool SSPPresent = true;
  art::Handle<artdaq::Fragments> SSPrawFragments;
  evt.getByLabel(fSSPRawDataLabel, fSSPFragType, SSPrawFragments);
  
  try { SSPrawFragments->size(); }
  catch(std::exception e) {
    mf::LogWarning("SSPToOffline") << "WARNING: Raw SSP data not found in event " << evt.event();
    SSPPresent = false;
  }
  
  if (SSPPresent) {
    if(!SSPrawFragments.isValid()){
      mf::LogError("SSPToOffline") << "Run: " << evt.run() << ", SubRun: " << evt.subRun() << ", Event: " << evt.event() << " is NOT VALID";
      throw cet::exception("raw NOT VALID");
      return;
    }
    
    artdaq::Fragments const& rawFragmentsSSP = *SSPrawFragments;
    DAQToOffline::GetSSPFirstTimestamp ( rawFragmentsSSP, nSSPPayloads, SSPTime );
    //std::cout << "Got SSP start time, it is " << SSPTime << std::endl;

    // Checking whether I have any waveforms......
    std::vector<raw::OpDetWaveform> waveforms = sspReform.SSPFragmentToOpDetWaveform(rawFragmentsSSP);
    if ( waveforms.size() ) {
      std::cout << "Identifier:Looking at event " << evt.event() << ", I have a vector of waveforms which has size " << waveforms.size() << std::endl; 
      std::cout << "Identifier:This event had " << total_Hits << " counter + trigger words " << std::endl;
    }
  } // SSP Present
  
  //---------------- RCE Good Event Timing stuff ----------------
  bool RCEPresent = true;
  art::Handle<artdaq::Fragments> RCErawFragments;
  evt.getByLabel(fRCERawDataLabel, fRCEFragType, RCErawFragments);
  
  try { RCErawFragments->size(); }
  catch(std::exception e) {
    mf::LogWarning("MuonCounter") << "WARNING: Raw RCE data not found in event " << evt.event() << std::endl;
    RCEPresent = false;
  }
  
  if (RCEPresent) {
    if(!RCErawFragments.isValid()){
      mf::LogError("SSPToOffline") << "Run: " << evt.run() << ", SubRun: " << evt.subRun() << ", Event: " << evt.event() << " is NOT VALID" << std::endl;
      throw cet::exception("RCErawFragments NOT VALID");
      return;
    }
    artdaq::Fragments const& rawFragmentsRCE = *RCErawFragments;
    DAQToOffline::GetRCEFirstTimestamp ( rawFragmentsRCE, ConsistRCE, NumADCs, RCETime );
    sumNADCs += NumADCs;
    if (NumADCs)
      std::cout << "Identifier:I had " << NumADCs << " ADCs in event " << evt.event() << ", first timestamp in this event is " << RCETime << "\n" << std::endl;
  } //RCEPresent
  //std::cout << "Got RCE start time, it is " << RCETime << std::endl;
  if (ConsistRCE > -1 ) ++nRCEPayloads;
  if (ConsistRCE == 1 ) ++nConsistRCEPayloads;
  
  // ---------------- Get the Good Event stuff ready for this event ----------------
  if (RCETime && PTBTime) RCE_PTB_diff = RCETime - PTBTime;
  if (RCETime && SSPTime) RCE_SSP_diff = RCETime - SSPTime;
  if (SSPTime && PTBTime) SSP_PTB_diff = SSPTime - PTBTime;
  
  if (RCETime && SSPTime && PTBTime) { // Check that all components are synchronous.
    RCE_PTB_diff = RCETime - PTBTime;
    RCE_SSP_diff = RCETime - SSPTime;
    if ( RCE_PTB_diff == 0 && RCE_SSP_diff == 0 ) ++nSynchronousEvents;
  }
  
  if (ConsistRCE == 1) {
    //std::cout << "!!!Got consistent RCEs so adding the trigger numbers" << std::endl;
    nPTBTrigsOn110 += ThisEv110;
    nPTBTrigsOn111 += ThisEv111;
    nPTBTrigsOn112 += ThisEv112;
    nPTBTrigsOn113 += ThisEv113;
    nPTBTrigsOn114 += ThisEv114;
    nPTBTrigsOn115 += ThisEv115;
  }
  if (fNevents%1000 == 0) 
    std::cout << "Looking at event " << evt.event() << " it had " << RCETime << " " << SSPTime << " " << PTBTime << " " << NumADCs 
	      << " RCE_SSP " << RCE_SSP_diff << " RCE_PTB " << RCE_PTB_diff << " SSP_PTB " << SSP_PTB_diff << ". "
	      << "So far have had " << nPTBTrigsOn110 << " " << nPTBTrigsOn111 << " " << nPTBTrigsOn112 << " " << nPTBTrigsOn113 << " " << nPTBTrigsOn114 << " " << nPTBTrigsOn115 << " trigs on each coincidence"
	      << ", and " << sumNADCs << " ADCs." << std::endl;
}

void Muoncounter::beginJob() {
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fHist1 = tfs->make<TH1D>("h1", "h1", u_TSU-l_TSU+1+4, l_TSU, u_TSU+1+4); 
  fHist2 = tfs->make<TH1D>("h2", "h2", u_BSU-l_BSU+1, l_BSU, u_BSU+1); 
  fHist3 = tfs->make<TH1D>("h3", "h3", u_Trig-l_Trig, l_Trig, u_Trig);
  
  fGoodRunHisto = tfs->make<TH1D>("GoddRunHisto","GoodRunHisto", 12, 0, 12);
  fGoodRunHisto->GetXaxis()->SetBinLabel(1 ,"PTB payload ratio");
  fGoodRunHisto->GetXaxis()->SetBinLabel(2 ,"SSP payload ratio");
  fGoodRunHisto->GetXaxis()->SetBinLabel(3 ,"RCE payload ratio");
  fGoodRunHisto->GetXaxis()->SetBinLabel(4 ,"Consistent RCE ratio");
  fGoodRunHisto->GetXaxis()->SetBinLabel(5 ,"Synchronous event ratio");
  fGoodRunHisto->GetXaxis()->SetBinLabel(6 ,"Trigs on Chan 110");
  fGoodRunHisto->GetXaxis()->SetBinLabel(7 ,"Trigs on Chan 111");
  fGoodRunHisto->GetXaxis()->SetBinLabel(8 ,"Trigs on Chan 112");
  fGoodRunHisto->GetXaxis()->SetBinLabel(9 ,"Trigs on Chan 113");
  fGoodRunHisto->GetXaxis()->SetBinLabel(10,"Trigs on Chan 114");
  fGoodRunHisto->GetXaxis()->SetBinLabel(11,"Trigs on Chan 115");
  fGoodRunHisto->GetXaxis()->SetBinLabel(12,"Total ADCs");
  

  //Making the Nearline header information tree
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
}

void Muoncounter::endJob()
{
  double total_Time = (double)(event_Length*event_Count);
  gStyle->SetOptStat(0);

  fHist1->Sumw2();
  fHist2->Sumw2();
  fHist3->Sumw2();
  if(total_Time > 0.0) fHist1->Scale(1/total_Time);
  if(total_Time > 0.0) fHist2->Scale(1/total_Time);
  if(total_Time > 0.0) fHist3->Scale(1/total_Time);

  TString fHist1_Title = Form("TSU Frequency");
  TString fHist2_Title = Form("BSU Frequency");
  TString fHist3_Title = Form("Trigger Frequency");
  TString fHist1_Name  = Form("TSUs");
  TString fHist2_Name  = Form("BSUs");
  TString fHist3_Name  = Form("Triggers");

  fHist1->SetTitle(fHist1_Title);
  fHist1->SetName(fHist1_Name);
  fHist2->SetTitle(fHist2_Title);
  fHist2->SetName(fHist2_Name);
  fHist3->SetTitle(fHist3_Title);
  fHist3->SetName(fHist3_Name);

  fHist1->GetXaxis()->SetTitle("Counter Number");
  fHist1->GetYaxis()->SetTitle("Frequency, [Hz]. (No Hits/Run Time)");
  fHist1->GetYaxis()->SetTitleOffset(1.3);

  fHist2->GetXaxis()->SetTitle("Counter Number");
  fHist2->GetYaxis()->SetTitle("Frequency, [Hz]. (No Hits/Run Time)");
  fHist2->GetYaxis()->SetTitleOffset(1.3);

  fHist3->GetXaxis()->SetTitle("Counter Number");
  fHist3->GetYaxis()->SetTitle("Frequency, [Hz]. (No Hits/Run Time)");
  fHist3->GetYaxis()->SetTitleOffset(1.3);

  for(int i = 28; i < 38; i++)
  {
    TString label = Form("WU%i", i-27);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 22; i < 28; i++)
  {
    TString label = Form("NU%i", i-21);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 0; i < 6; i++)
  {
    TString label = Form("SL%i", i+1);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 16; i < 22; i++)
  {
    TString label = Form("NL%i", i-15);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 1; i < 5; i++)
  {
    TString label = Form("XX%i", i);
    fHist1->GetXaxis()->SetBinLabel(i+44, label);
  }

  fHist1->GetXaxis()->SetBinLabel(15+1, "EL1");
  fHist1->GetXaxis()->SetBinLabel(14+1, "EL2");
  fHist1->GetXaxis()->SetBinLabel(13+1, "EL3");
  fHist1->GetXaxis()->SetBinLabel(12+1, "EL4");
  fHist1->GetXaxis()->SetBinLabel(11+1, "EL5");
  fHist1->GetXaxis()->SetBinLabel(10+1, "EL6");
  fHist1->GetXaxis()->SetBinLabel(9+1,  "EL7");
  fHist1->GetXaxis()->SetBinLabel(8+1,  "EL8");
  fHist1->GetXaxis()->SetBinLabel(7+1,  "EL9");
  fHist1->GetXaxis()->SetBinLabel(6+1,  "EL10");

  fHist1->GetXaxis()->SetBinLabel(43+1, "SU1");
  fHist1->GetXaxis()->SetBinLabel(42+1, "SU2");
  fHist1->GetXaxis()->SetBinLabel(41+1, "SU3");
  fHist1->GetXaxis()->SetBinLabel(40+1, "SU4");
  fHist1->GetXaxis()->SetBinLabel(39+1, "SU5");
  fHist1->GetXaxis()->SetBinLabel(38+1, "SU6");
  fHist1->GetXaxis()->SetLabelSize(0.018);

  for(int i = 67; i < 83; i++)
  {
    TString label = Form("RM%i", i-66);
    fHist2->GetXaxis()->SetBinLabel(i-43, label);
  }
  for(int i = 83; i < 93; i++)
  {
    TString label = Form("RL%i", i-82);
    fHist2->GetXaxis()->SetBinLabel(i-43, label);
  }
  for(int i = 44; i < 57; i++)
  {
    TString label = Form("CL%i", i-43);
    fHist2->GetXaxis()->SetBinLabel(i-43, label);
  }

  fHist2->GetXaxis()->SetBinLabel(66+1-44, "CU1");
  fHist2->GetXaxis()->SetBinLabel(65+1-44, "CU2");
  fHist2->GetXaxis()->SetBinLabel(64+1-44, "CU3");
  fHist2->GetXaxis()->SetBinLabel(63+1-44, "CU4");
  fHist2->GetXaxis()->SetBinLabel(62+1-44, "CU5");
  fHist2->GetXaxis()->SetBinLabel(61+1-44, "CU6");
  fHist2->GetXaxis()->SetBinLabel(60+1-44, "CU7");
  fHist2->GetXaxis()->SetBinLabel(59+1-44, "CU8");
  fHist2->GetXaxis()->SetBinLabel(58+1-44, "CU9");
  fHist2->GetXaxis()->SetBinLabel(57+1-44, "CU10");
  fHist2->GetXaxis()->SetLabelSize(0.025);

  // -------------------- Stuff for Good Events List ---------------------
  double PTBPayloadRat = nPTBPayloads / (double)fNevents;
  double SSPPayloadRat = nSSPPayloads / (double)fNevents;
  double RCEPayloadRat = nRCEPayloads / (double)fNevents;
  double ConsistRCERat = nConsistRCEPayloads / (double)nRCEPayloads;
  double SynchronRat   = nSynchronousEvents / (double)fNevents;

  fGoodRunHisto->SetBinContent(1 , PTBPayloadRat);
  fGoodRunHisto->SetBinContent(2 , SSPPayloadRat);
  fGoodRunHisto->SetBinContent(3 , RCEPayloadRat);
  fGoodRunHisto->SetBinContent(4 , ConsistRCERat);
  fGoodRunHisto->SetBinContent(5 , SynchronRat);
  fGoodRunHisto->SetBinContent(6 , nPTBTrigsOn110);
  fGoodRunHisto->SetBinContent(7 , nPTBTrigsOn111);
  fGoodRunHisto->SetBinContent(8 , nPTBTrigsOn112);
  fGoodRunHisto->SetBinContent(9 , nPTBTrigsOn113);
  fGoodRunHisto->SetBinContent(10, nPTBTrigsOn114);
  fGoodRunHisto->SetBinContent(11, nPTBTrigsOn115);
  fGoodRunHisto->SetBinContent(12, sumNADCs);
  

  std::cout << "IDENTIFIER: Run: " << fRun << " SubRun " << fSubrun << " has " << fNevents << " events in total" 
	    << ", ratio that have PTB payloads " << PTBPayloadRat 
	    << ", ratio that have SSP payloads " << SSPPayloadRat 
	    << ", ratio that have RCE payloads " << RCEPayloadRat << ", ratio which were consistent " << ConsistRCERat << ", I had a total of " << sumNADCs <<  " ADCs"
	    << ", ratio of synchronous events " << SynchronRat
	    << ". I had " << nPTBTrigsOn110 << " " << nPTBTrigsOn111 << " " << nPTBTrigsOn112 << " " << nPTBTrigsOn113 << " " << nPTBTrigsOn114 << " " << nPTBTrigsOn115 << " trigs on each coincidence."
	    << std::endl;

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


DEFINE_ART_MODULE(Muoncounter)

}  // namespace Muoncounter

#endif // Muoncounter_Module
