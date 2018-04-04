////////////////////////////////////////////////////////////////////////
// Class:       Muoncountertimes
// Module Type: analyzer
// File:        Muoncountertimes_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
////////////////////////////////////////////////////////////////////////
#ifndef Muoncountertimes_Module
#define Muoncountertimes_Module

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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "dune/daqinput35t/PennToOffline.h"
#include "dune/daqinput35t/utilities/UnpackFragment.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TString.h"

//PTB online to offline converter alg


//const int kMaxHits       = 10000; //maximum number of hits // unused
//const int kMaxAuxDets = 100; // unused
//const unsigned short kMaxTkIDs = 100;
namespace MyMuoncountertimes {

class Muoncountertimes : public art::EDAnalyzer {
public:
  explicit Muoncountertimes(fhicl::ParameterSet const & p);
  virtual ~Muoncountertimes();

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
  
  //int total_Hits; // unused
  int run;
  int subrun;
  int event;
  int event_Count = 0;
  //double event_Length = 10e-3; // unused

  std::vector<TH1D*> fHist_Vec;

  int         fLower_Block1;
  int         fUpper_Block1;
  int         fLower_Block2;
  int         fUpper_Block2;
  int         fCoincidence_Window;
  double      fCombinedTimeDelay;
  std::string fCounterModuleLabel;
  std::string fFragType;
  std::string fRawDataLabel;
  std::string fPTBMapFile;
  std::string fPTBMapDir;

  std::map<int,int> fPTBMap;
};

Muoncountertimes::Muoncountertimes(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
{
   // Read in the parameters from the .fcl file.
    this->reconfigure(pset);
}
  
void Muoncountertimes::reconfigure(fhicl::ParameterSet const& p)
{   // Read parameters from the .fcl file. The names in the arguments
  // to p.get<TYPE> must match names in the .fcl file.

  fLower_Block1       = p.get<int>("LowerBlock1");
  fUpper_Block1       = p.get<int>("UpperBlock1");
  fLower_Block2       = p.get<int>("LowerBlock2");
  fUpper_Block2       = p.get<int>("UpperBlock2");
  fCoincidence_Window = p.get<int>("CoincidenceWindow");
  fCounterModuleLabel = p.get< std::string >("CounterModuleLabel");
  //    fLArG4ModuleLabel(p.get< std::string >("LArGeantModuleLabel", "largeant"));
  fCombinedTimeDelay = p.get< double >("CombinedTimeDelay");
  fFragType = p.get<std::string>("FragType");
  fRawDataLabel = p.get<std::string>("DataLabel");
  fPTBMapFile   = p.get<std::string>("PTBMapFile");
  fPTBMapDir    = p.get<std::string>("PTBMapDir");

  DAQToOffline::BuildPTBChannelMap(fPTBMapDir, fPTBMapFile, fPTBMap);

  return;
}

Muoncountertimes::~Muoncountertimes()
{
  // Clean up dynamic memory and other resources here.
}

void Muoncountertimes::analyze(const art::Event& evt)
{
  // Implementation of required member function here.
  ResetVars();

  art::ServiceHandle<geo::Geometry> geom;

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
    mf::LogWarning("MuonCounterTimes") << "WARNING: Raw PTB data not found in event " << evt.event();
    return;
  }

 // Check that the data is valid
  if(!fragment.isValid()){
    mf::LogError("MuonCounterTimes") << "Run: " << evt.run()
				 << ", SubRun: " << evt.subRun()
				 << ", Event: " << evt.event()
				 << " is NOT VALID";
    throw cet::exception("raw NOT VALID");
    return;
  }

  std::vector<std::vector<int>>       coincidences(2);
  std::vector<std::vector<int>>       trig_ID_Vec(2);
  std::vector<std::vector<long long>> trig_Times_Vec(2); 
  lbne::PennMicroSlice::Payload_Timestamp *FirstPTBTimestamp = nullptr;
  auto trigs = DAQToOffline::PennFragmentToExternalTrigger(*fragment, fPTBMap, FirstPTBTimestamp);
  
  unsigned int total_Hits = trigs.size();

  for(unsigned int i = 0; i < total_Hits; i++)
  {
    int       auxdetid  = trigs.at(i).GetTrigID();
    long long trig_Time = trigs.at(i).GetTrigTime();   

    if(auxdetid>=fLower_Block1 && auxdetid<=fUpper_Block1)
    {
      trig_ID_Vec.at(0).push_back(auxdetid); 
      trig_Times_Vec.at(0).push_back(trig_Time);
    }
    else if(auxdetid>=fLower_Block2 && auxdetid<=fUpper_Block2)
    {
      trig_ID_Vec.at(1).push_back(auxdetid); 
      trig_Times_Vec.at(1).push_back(trig_Time);
    }
  }

  for(unsigned int i = 0; i < trig_ID_Vec.at(0).size(); i++)
  {
    for(unsigned int j = 0; j < trig_ID_Vec.at(1).size(); j++)
    {
      if(j<=i)
      {
        if(std::abs(trig_Times_Vec.at(0).at(i)-trig_Times_Vec.at(1).at(j))<fCoincidence_Window)
        {
          coincidences.at(0).push_back(trig_ID_Vec.at(0).at(i));
          coincidences.at(1).push_back(trig_ID_Vec.at(1).at(j));
        }
      }
    }
  }

  for(unsigned int i = 0; i < coincidences.at(0).size(); i++)
  {
    fHist_Vec.at(coincidences.at(0).at(i)-fLower_Block1)->Fill(coincidences.at(1).at(i));
  }
  event_Count+=1;
}

void Muoncountertimes::beginJob()
{
  // Implementation of optional member function here.

  fHist_Vec.resize(fUpper_Block1-fLower_Block1+1);

  art::ServiceHandle<art::TFileService> tfs;

  for(unsigned int i = 0; i < fHist_Vec.size(); i++)
  {
    TString key     = Form("h%i", i);
    TH1D *tempHist  = tfs->make<TH1D>(key, key, fUpper_Block2-fLower_Block2+1, fLower_Block2, fUpper_Block2+1); 
    fHist_Vec.at(i) = tempHist; 
  }
}

void Muoncountertimes::endJob()
{
  for(unsigned int i = 0; i < fHist_Vec.size(); i++)
  {
    TString fHist_Title1 = Form("Coincidence Plots, run %i, ", run); 
    TString fHist_Title2 = Form("counter %i.", i+fLower_Block1);
    TString fHist_Title  = fHist_Title1 + fHist_Title2;

    TString fHist_Name1 = Form("Run%i", run);
    TString fHist_Name2 = Form("Counter%i", i+fLower_Block1);
    TString fHist_Name  = fHist_Name1 + fHist_Name2;
 
    fHist_Vec.at(i)->GetXaxis()->SetTitle("Coincidence Partner");
    fHist_Vec.at(i)->GetYaxis()->SetTitle("Frequency");
    fHist_Vec.at(i)->GetYaxis()->SetTitleOffset(1.);

    fHist_Vec.at(i)->SetTitle(fHist_Title);
    fHist_Vec.at(i)->SetName(fHist_Name);

    fHist_Vec.at(i)->Sumw2();
  }
}

void Muoncountertimes::ResetVars(){

  run = -99999;
  subrun = -99999;
  event = -99999;
}

DEFINE_ART_MODULE(Muoncountertimes)

  }  // namespace Muoncountertimes

#endif // Muoncountertimes_Module
