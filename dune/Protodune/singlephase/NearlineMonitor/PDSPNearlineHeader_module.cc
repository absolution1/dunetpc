// Header module file. For now saves the run and subrun ids of the input file into a TTree

#ifndef PDSPNearlineHeader_module
#define PDSPNearlineHeader_module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

// C++ Includes
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace PDSPNearlineheader_module{
  
  class PDSPNearlineHeaderModule : public art::EDAnalyzer{
  public:
    
    explicit PDSPNearlineHeaderModule(fhicl::ParameterSet const& pset);
    virtual ~PDSPNearlineHeaderModule();
    
    void beginJob();
    void analyze(const art::Event& evt);
    
  private:
    
    int fRun;
    int fSubRun;
    uint64_t fTimeStamp;
    
    TTree *fPDSPNearlineHeaderTree;
    
    std::vector<TString> runset;
    
  };
  

  //-----------------------------------------------------------------------
  PDSPNearlineHeaderModule::PDSPNearlineHeaderModule(fhicl::ParameterSet const& pset) : EDAnalyzer(pset){
    
    art::ServiceHandle<art::TFileService> tfs;
    
    // Define output tree
    fPDSPNearlineHeaderTree = tfs->make<TTree>("PDSPNearlineHeader", "PDSP Nearline header tree");
    fPDSPNearlineHeaderTree->Branch("fRun",       &fRun,       "fRun/I");
    fPDSPNearlineHeaderTree->Branch("fSubRun",    &fSubRun,    "fSubRun/I");
    fPDSPNearlineHeaderTree->Branch("fTimeStamp", &fTimeStamp, "fTimeStamp/l");
  }
  
  //-----------------------------------------------------------------------
  PDSPNearlineHeaderModule::~PDSPNearlineHeaderModule(){}
  
  //-----------------------------------------------------------------------
  void PDSPNearlineHeaderModule::beginJob() {}
  
  //-----------------------------------------------------------------------
  void PDSPNearlineHeaderModule::analyze(const art::Event& evt){
    
    fRun = evt.run();
    fSubRun = evt.subRun();
    fTimeStamp = evt.time().value();
    
    TString tempstring = Form("%i-%i",fRun,fSubRun);
    
    int num_items = std::count(runset.begin(), runset.end(), tempstring);
    if(num_items == 0){ // run-subrun not in vector - add it and save
      runset.push_back(tempstring);
      fPDSPNearlineHeaderTree->Fill();
    }
  } 
} // namespace

DEFINE_ART_MODULE(PDSPNearlineheader_module::PDSPNearlineHeaderModule)

#endif
