#ifndef disambigcheck_h
#define disambigcheck_h
////////////////////////////////////////////////////////////////////////
//
// Disambigcheck class
// 
// tjyang@fnal.gov
//  
// Based on Tom Junk's idea of 3-view matching
// HitFinder for 35t, input is undisambiguated hits, out is disambiguated hits
// Start from code written by talion@gmail.com
//
//
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <stdint.h>
#include <string>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"   
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDProducer.h" 
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft Includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "DisambigAlg35t.h"


// ROOT Includes 
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TVector2.h"
#include "TVector3.h"

namespace disambigcheck{
  class disambigcheck : public art::EDAnalyzer
  {
    
  public:
    
    explicit disambigcheck(fhicl::ParameterSet const& pset); 
    
    void analyze(const art::Event& evt); 
    void beginJob(); 
    void beginRun(const art::Run& run);
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& pset);                
    
     virtual ~disambigcheck();
    
  private:
    
   art::ServiceHandle<geo::Geometry> fGeom;
    
    std::string fChanHitDisambig;
    std::string fChanHitCheater;
    
    TH1D* fCorrect;
    TH1D* fMissed;
    TH1D* fIncorrect;
    
  protected: 
    
    
  }; 
  
  //-------------------------------------------------
  //-------------------------------------------------
  disambigcheck::disambigcheck(fhicl::ParameterSet const& pset) 
    : EDAnalyzer (pset)  
  {
    this->reconfigure(pset);
    //  produces< std::vector<recob::Hit> >();
  }
  
  disambigcheck::~disambigcheck() {}
  //-------------------------------------------------
  //-------------------------------------------------
  void disambigcheck::reconfigure(fhicl::ParameterSet const& pset)
  {
    
    fChanHitDisambig =  pset.get< std::string >("ChanHitDisambig");
    fChanHitCheater =  pset.get< std::string >("ChanHitCheater");
    
  }  
  
  //-------------------------------------------------
  //-------------------------------------------------
  void disambigcheck::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    
    fCorrect = tfs->make<TH1D>("correct","Correct disambiguated",     100, 0, 1.01);
    fMissed = tfs->make<TH1D>("missed","Missed hits",     100, 0, 1.01);
    fIncorrect = tfs->make<TH1D>("incorrect","Incorrect disambiguated",     100, 0, 1.01);
  }
  
  //-------------------------------------------------
  //-------------------------------------------------
  void disambigcheck::beginRun(const art::Run& run){}
  void disambigcheck::endJob()
  {
    
  }
  

  //-------------------------------------------------
  void disambigcheck::analyze(const art::Event& evt)
  {
    
    std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>);
    art::Handle< std::vector<recob::Hit> > ChannelHitsDisambig;
    art::Handle< std::vector<recob::Hit> > ChannelHitsCheater;
    evt.getByLabel(fChanHitDisambig, ChannelHitsDisambig);
    evt.getByLabel(fChanHitCheater, ChannelHitsCheater);
    
    // Make unambiguous collection hits
    std::vector< art::Ptr<recob::Hit> >  ChHitsDisambig;
    std::vector< art::Ptr<recob::Hit> >  ChHitsCheater;
    art::fill_ptr_vector(ChHitsDisambig, ChannelHitsDisambig);
    art::fill_ptr_vector(ChHitsCheater, ChannelHitsCheater);

    //   int cheathits = ChHitsCheater.size();
    int correcthits = 0;
    int incorrecthits = 0; 
    int missinghits = 0;
    //  int disambighits = ChHitsDisambig.size();
      
   

    for( size_t h = 0; h < ChHitsCheater.size(); h++ ){
      uint32_t cheatchannel = ChHitsCheater[h]->Channel();
      double cheatpeaktime = ChHitsCheater[h]->PeakTime();
      uint32_t cheatwire = ChHitsCheater[h]->WireID().Wire;
  
      bool found = false;

      if(ChHitsCheater[h]->View() == geo::kZ ) continue;

      for( size_t w = 0; w < ChHitsDisambig.size(); w++ ) {
	uint32_t disambigchannel = ChHitsDisambig[w]->Channel();
	double disambigpeaktime = ChHitsDisambig[w]->PeakTime();
	uint32_t disambigwire = ChHitsDisambig[w]->WireID().Wire;

	
	if((cheatchannel == disambigchannel)&&(TMath::Abs(cheatpeaktime - disambigpeaktime)< 0.1)) {
            
            found=true;
	    if(cheatwire==disambigwire){
	      correcthits++;
	    
	    }
	    else{
	      incorrecthits++;
	    }
	    
	    
	  }	
	
	
      }
      
      if(!found){
	missinghits++;
      }
    }    
    
      double dcorrecthits = correcthits;
      double dincorrecthits = incorrecthits;
      double dmissinghits = missinghits;
      //   double ddisambighits = disambighits;
      
      double dtotalhits = (dcorrecthits + dincorrecthits + dmissinghits);    

      std::cout << "Total hits " << dtotalhits << ", correct hits " << dcorrecthits << ", incorrect hits " << dincorrecthits << ", missing hits " << dmissinghits << std::endl;
      
      //double disambighitsFraction = ddisambighits / dtotalhits;
      double correcthitsFraction = dcorrecthits / dtotalhits;
      double incorrecthitsFraction = dincorrecthits / dtotalhits; 
      double missinghitsFraction = dmissinghits / dtotalhits; 
      


    fCorrect->Fill(correcthitsFraction);
    fMissed->Fill(missinghitsFraction);
    fIncorrect->Fill(incorrecthitsFraction);
    
  }
  
  DEFINE_ART_MODULE(disambigcheck)
 
}
  
  
  
#endif
