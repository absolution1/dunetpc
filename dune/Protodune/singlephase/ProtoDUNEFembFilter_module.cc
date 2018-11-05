// A module to filter out events with out full set of fembs (on beam side?)
// owen.goodwin@postgrad.manchester.ac.uk

#include <iostream>
#include <utility>
#include <set>

#include "TH1.h"
#include "TFile.h"

#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Services/Optional/TFileService.h"

#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"
#include "lardataobj/RawData/RDTimeStamp.h"





namespace filt{

  class ProtoDUNEFembFilter : public art::EDFilter {
  public:
    explicit ProtoDUNEFembFilter(fhicl::ParameterSet const & pset);
    virtual ~ProtoDUNEFembFilter() {};
    virtual bool filter(art::Event& e);
    void beginJob();

  private:

    unsigned int fLogLevel;
    
    TH1D* fSelectedEvents;
    TH1D* fTotalEvents;

  };

  ProtoDUNEFembFilter::ProtoDUNEFembFilter::ProtoDUNEFembFilter(fhicl::ParameterSet const & pset) {

    fLogLevel = pset.get<unsigned int>("LogLevel");

    if ( fLogLevel >= 1 ) {
      std::cout  << "                LogLevel: " << fLogLevel << std::endl;

      }
    
  }

  void ProtoDUNEFembFilter::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;
      fSelectedEvents = tfs->make<TH1D>("fSelectedEvents", "Number of Selected Events", 3, 0, 3); //counts the number of selected events 
      fTotalEvents = tfs->make<TH1D>("fTotalEvents", "Total Events", 3, 0, 3); //counts the initial number of events in the unfiltered root input file
    
  }

  bool ProtoDUNEFembFilter::filter(art::Event & evt) {
    std::vector<int> BeamsideAPAs;
    // Add some elements to myIntVector
    BeamsideAPAs.push_back(0);
    BeamsideAPAs.push_back(2);
    BeamsideAPAs.push_back(4);

    const std::string myname = "ProtoDUNEFembFilter::filter: ";

    bool keep = true;
    // Helper utility functions
    protoana::ProtoDUNEDataUtils dataUtil;
    //ProtoDUNEDataUtils dataUtil;

    fTotalEvents->Fill(1); //count total events
    for (auto APA = BeamsideAPAs.begin(); APA != BeamsideAPAs.end(); ++APA){ //loop through beam side APAs
      //std::cout<<"APA:"<<*APA<<std::endl;
      //std::cout<<dataUtil.GetNActiveFembsForAPA(evt, *APA)<<std::endl;
      if (dataUtil.GetNActiveFembsForAPA(evt, *APA)!=20){ //check if APA has all 20 fembs active

        if (fLogLevel >=2) std::cout<<"Missing FEMBs on APA: "<<*APA<<std::endl; 
        keep=false; //if not remove event
      }
      
    
    
      
    }
    if ( fLogLevel >=2 ) std::cout << myname << (keep ? "Keep" : "Reject") << "ing event." << std::endl;
    if (keep==true) fSelectedEvents->Fill(1); //count total events
    
    return keep;

  }

  DEFINE_ART_MODULE(ProtoDUNEFembFilter)

}