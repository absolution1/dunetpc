// Header module file. For now saves the run and subrun ids of the input file into a TTree

#ifndef Header_module
#define Header_module

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

namespace header_module{
  
  class HeaderModule : public art::EDAnalyzer{
  public:
    
    explicit HeaderModule(fhicl::ParameterSet const& pset);
    virtual ~HeaderModule();
    
    void beginJob();
    void analyze(const art::Event& evt);
    
  private:
    
    int fRun;
    int fSubRun;
    
    TTree *fHeaderTree;
    
    std::vector<TString> runset;
    
  };
  

  //-----------------------------------------------------------------------
  HeaderModule::HeaderModule(fhicl::ParameterSet const& pset) : EDAnalyzer(pset){
    
    art::ServiceHandle<art::TFileService> tfs;
    
    // Define output tree
    fHeaderTree = tfs->make<TTree>("Header", "header tree");
    fHeaderTree->Branch("fRun",    &fRun,    "fRun/I");
    fHeaderTree->Branch("fSubRun", &fSubRun, "fSubRun/I");
  }
  
  //-----------------------------------------------------------------------
  HeaderModule::~HeaderModule(){}
  
  //-----------------------------------------------------------------------
  void HeaderModule::beginJob() {}
  
  //-----------------------------------------------------------------------
  void HeaderModule::analyze(const art::Event& evt){
    
    fRun = evt.run();
    fSubRun = evt.subRun();
    
    TString tempstring = Form("%i-%i",fRun,fSubRun);
    
    int num_items = std::count(runset.begin(), runset.end(), tempstring);
    if(num_items == 0){ // run-subrun not in vector - add it and save in output
      runset.push_back(tempstring);
      fHeaderTree->Fill();
    }
  } 
} // namespace

DEFINE_ART_MODULE(header_module::HeaderModule)

#endif
