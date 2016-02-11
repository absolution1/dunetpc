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
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/Utilities/LArProperties.h"
#include "lardata/Utilities/DetectorProperties.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/TimeService.h"
#include "larsim/Simulation/AuxDetSimChannel.h"
#include "lardata/RawData/ExternalTrigger.h"
#include "dune/daqinput35t/PennToOffline.h"
#include "dune/daqinput35t/utilities/UnpackFragment.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TString.h"
#include "TText.h"

#include "dune/NearlineMonitor/NearlineVersion.h"

const int kMaxHits       = 10000; //maximum number of hits
const int kMaxAuxDets = 100;

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

  
  int total_Hits;
  int l_TSU = 0;   
  int u_TSU = 47;  
  int l_BSU = 48;  
  int u_BSU = 100; 
  int l_Trig = 101;
  int u_Trig = 120;
  int event_Count = 0;
  double event_Length = 10e-3;

  TH1D *fHist1;
  TH1D *fHist2;
  TH1D *fHist3;

  double      fCombinedTimeDelay;
  std::string fCounterModuleLabel;
  std::string fFragType;
  std::string fRawDataLabel;

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
  fEndTime(0)
{
   // Read in the parameters from the .fcl file.
    this->reconfigure(pset);
}
  
void Muoncounter::reconfigure(fhicl::ParameterSet const& p)
{   // Read parameters from the .fcl file. The names in the arguments
  // to p.get<TYPE> must match names in the .fcl file.

  fCounterModuleLabel = p.get< std::string >("CounterModuleLabel");
  //    fLArG4ModuleLabel(p.get< std::string >("LArGeantModuleLabel", "largeant"));
  fCombinedTimeDelay = p.get< double >("CombinedTimeDelay");
  fFragType          = p.get<std::string>("FragType");
  fRawDataLabel      = p.get<std::string>("DataLabel");
  fPTBMapFile        = p.get<std::string>("PTBMapFile");
  fPTBMapDir         = p.get<std::string>("PTBMapDir");

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

  // New way to get counter hits
  /*
  art::Handle< std::vector< raw::ExternalTrigger> > externalTriggerListHandle;
  evt.getByLabel(fCounterModuleLabel, externalTriggerListHandle);
  std::vector< art::Ptr< raw::ExternalTrigger> > trigs;
  art::fill_ptr_vector(trigs,externalTriggerListHandle);
 */

  art::Handle< artdaq::Fragments > fragment;

  evt.getByLabel(fRawDataLabel, fFragType, fragment);

  // Check if there is PTB data in this event
  // Don't crash code if not present, just don't save anything
  try { fragment->size(); }
  catch(std::exception e) {
    mf::LogWarning("MuonCounter") << "WARNING: Raw PTB data not found in event " << evt.event();
    return;
  }

  // Check that the data is valid
  if(!fragment.isValid()){
    mf::LogError("MuonCounter") << "Run: " << evt.run()
				 << ", SubRun: " << evt.subRun()
				 << ", Event: " << evt.event()
				 << " is NOT VALID";
    throw cet::exception("raw NOT VALID");
    return;
  }

  auto trigs = DAQToOffline::PennFragmentToExternalTrigger(*fragment, fPTBMap);

  unsigned int total_Hits = trigs.size();

  for(unsigned int i = 0; i < total_Hits; i++)
  {
    int auxdetid = trigs.at(i).GetTrigID();

    if(auxdetid<=u_TSU)
    {
      fHist1->Fill(auxdetid);
    }
    else if(auxdetid>=l_BSU && auxdetid<=u_BSU)
    {
      fHist2->Fill(auxdetid);
    }
    else
    {
      fHist3->Fill(auxdetid);
    }
  }
  event_Count+=1;
}

void Muoncounter::beginJob()
{
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> tfs;
  fHist1 = tfs->make<TH1D>("h1", "h1", u_TSU-l_TSU+1, l_TSU, u_TSU+1); 
  fHist2 = tfs->make<TH1D>("h2", "h2", u_BSU-l_BSU, l_BSU, u_BSU); 
  fHist3 = tfs->make<TH1D>("h3", "h3", u_Trig-l_Trig, l_Trig, u_Trig);

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

  fHist1->Sumw2();
  fHist2->Sumw2();
  fHist3->Sumw2();
  fHist1->Scale(1/total_Time);
  fHist2->Scale(1/total_Time);
  fHist3->Scale(1/total_Time);

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

  for(int i = 0; i < 10; i++)
  {
    TString label = Form("WU%i", i);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 10; i < 20; i++)
  {
    TString label = Form("EL%i", i);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 20; i < 24; i++)
  {
    TString label = Form("EL Special%i", i);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 24; i < 30; i++)
  {
    TString label = Form("NU%i", i);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 30; i < 36; i++)
  {
    TString label = Form("SL%i", i);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 36; i < 42; i++)
  {
    TString label = Form("NL%i", i);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  for(int i = 42; i < 48; i++)
  {
    TString label = Form("SU%i", i);
    fHist1->GetXaxis()->SetBinLabel(i+1, label);
  }
  fHist1->GetXaxis()->SetLabelSize(0.018);

  for(int i = 48; i < 59; i++)
  {
    TString label = Form("Tel%i", i);
    fHist2->GetXaxis()->SetBinLabel(i-47, label);
  }
  for(int i = 59; i < 65; i++)
  {
    TString label = Form("CL%i", i);
    fHist2->GetXaxis()->SetBinLabel(i-47, label);
  }
  for(int i = 65; i < 74; i++)
  {
    TString label = Form("CU%i", i);
    fHist2->GetXaxis()->SetBinLabel(i-47, label);
  }
  for(int i = 74; i < 98; i++)
  {
    TString label = Form("Tel%i", i);
    fHist2->GetXaxis()->SetBinLabel(i-47, label);
  }
  fHist2->GetXaxis()->SetLabelSize(0.025);
  /*
  fHist3->GetXaxis()->SetBinLabel(, "Trigger 1");
  fHist3->GetXaxis()->SetBinLabel(, "Trigger 2");
  fHist3->GetXaxis()->SetBinLabel(, "Trigger 3");*/


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
