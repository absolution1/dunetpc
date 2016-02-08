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
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/TimeService.h"
#include "Simulation/AuxDetSimChannel.h"
#include "RawData/ExternalTrigger.h"
#include "dune/daqinput35t/PennToOffline.h"
#include "dune/daqinput35t/utilities/UnpackFragment.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TString.h"
#include "TText.h"

//PTB online to offline converter alg


const int kMaxHits       = 10000; //maximum number of hits
const int kMaxAuxDets = 100;
//const unsigned short kMaxTkIDs = 100;
namespace MyMuoncounter {

class Muoncounter : public art::EDAnalyzer {
public:
  explicit Muoncounter(fhicl::ParameterSet const & p);
  virtual ~Muoncounter();

   // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();
    void endJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
  //    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

private:

  void ResetVars();
  
  int total_Hits;
  int run;
  int subrun;
  int event;
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
};

Muoncounter::Muoncounter(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
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
  fFragType = p.get<std::string>("FragType");
  fRawDataLabel = p.get<std::string>("DataLabel");

  return;
}

Muoncounter::~Muoncounter()
{
  // Clean up dynamic memory and other resources here.
}

void Muoncounter::analyze(const art::Event& evt)
{
  // Implementation of required member function here.
  ResetVars();

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::TimeService> timeserv;

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

  auto trigs = DAQToOffline::PennFragmentToExternalTrigger(*fragment);

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

  TString fHist1_Title = Form("TSU Frequency, Run %i", run);
  TString fHist2_Title = Form("BSU Frequency, Run %i", run);
  TString fHist3_Title = Form("Trigger Frequency, Run %i", run);
  TString fHist1_Name  = Form("%i_TSUs", run);
  TString fHist2_Name  = Form("%i_BSUs", run);
  TString fHist3_Name  = Form("%i_Triggers", run);

  fHist1->SetTitle(fHist1_Title);
  fHist1->SetName(fHist1_Name);
  fHist2->SetTitle(fHist2_Title);
  fHist2->SetName(fHist2_Name);
  fHist3->SetTitle(fHist3_Title);
  fHist3->SetName(fHist3_Name);

  fHist1->GetXaxis()->SetTitle("Counter Number");
  fHist1->GetYaxis()->SetTitle("Frequency, [Hz]. (No Hits/Run Time)");
  fHist1->GetYaxis()->SetTitleOffset(1.3);

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
    TString label = Form("Telescope%i", i);
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
    TString label = Form("Telescope%i", i);
    fHist2->GetXaxis()->SetBinLabel(i-47, label);
  }
  fHist2->GetXaxis()->SetLabelSize(0.025);
  /*
  fHist3->GetXaxis()->SetBinLabel(, "Trigger 1");
  fHist3->GetXaxis()->SetBinLabel(, "Trigger 2");
  fHist3->GetXaxis()->SetBinLabel(, "Trigger 3");*/
}

void Muoncounter::ResetVars(){

  run = -99999;
  subrun = -99999;
  event = -99999;
}

DEFINE_ART_MODULE(Muoncounter)

  }  // namespace Muoncounter

#endif // Muoncounter_Module
