//Use this module to create 2D histos  raw event display for protoDUNE detector
//Dec. 2016: taken from RawEVD35t module.


#ifndef RawEventDisplay_module
#define RawEventDisplay_module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Data type includes
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"



// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes.
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<int> >+;
#endif

namespace raw_event_display{

  class RawEventDisplay : public art::EDAnalyzer{
  public:
 
    explicit RawEventDisplay(fhicl::ParameterSet const& pset);
    virtual ~RawEventDisplay();

    void beginJob();
    void beginRun(const art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt); 

  private:

    // Parameters in .fcl file
    std::string fRawDigitLabel;
    std::string fTPCInput;
    std::string fTPCInstance;

    // Branch variables for tree
    unsigned int fEvent;
    unsigned int fRun;
    unsigned int fSubRun;


    // TPC
    unsigned int fNUCh;
    unsigned int fNVCh;
    unsigned int fNZCh;

    // find channel boundaries for each view
    unsigned int fUChanMin;
    unsigned int fUChanMax;
    unsigned int fVChanMin;
    unsigned int fVChanMax;
    unsigned int fZChanMin;
    unsigned int fZChanMax;
    unsigned int fNticks;

    unsigned int fNofAPA;
    unsigned int fChansPerAPA;


    //unsigned int fMinT, fMaxT, fMaxTimeRange; // unused

    std::vector<TH2S*> fTimeChanU;
    std::vector<TH2S*> fTimeChanV;
    std::vector<TH2S*> fTimeChanZ;




    // define nADC counts for uncompressed vs compressed
    unsigned int nADC_uncompPed;


    geo::GeometryCore const * fGeom = &*(art::ServiceHandle<geo::Geometry>());


 }; // RawEventDisplay

  //-----------------------------------------------------------------------

  RawEventDisplay::RawEventDisplay(fhicl::ParameterSet const& parameterSet): EDAnalyzer(parameterSet) {
    this->reconfigure(parameterSet);

  }

  //-----------------------------------------------------------------------


  RawEventDisplay::~RawEventDisplay() {}
   

  //-----------------------------------------------------------------------
  void RawEventDisplay::beginJob() {
    // place to define the histograms

    art::ServiceHandle<art::TFileService> tfs;
    //Histogram names and titles                                                                                                                                                         
    std::stringstream  name, title;

    unsigned int UChMin;
    unsigned int UChMax;
    unsigned int VChMin;
    unsigned int VChMax;
    unsigned int ZChMin;
    unsigned int ZChMax;
    TH2S* TempHisto;


    // Accquiring geometry data
    fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2;
    fChansPerAPA = fGeom->Nchannels()/fNofAPA;

    // taken from dune35t module a way to organise the channel mapping:
    // loop through channels in the first APA to find the channel boundaries for each view
    // will adjust for desired APA after
    fUChanMin = 0;
    fZChanMax = fChansPerAPA - 1;
    for ( unsigned int c = fUChanMin + 1; c < fZChanMax; c++ ){
      if ( fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU ){
        fVChanMin = c;
        fUChanMax = c - 1;
      }
      if ( fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV ){
        fZChanMin = c;
        fVChanMax = c-1;
      }
    }
    

    fNUCh=fUChanMax-fUChanMin+1;
    fNVCh=fVChanMax-fVChanMin+1;
    fNZCh=fZChanMax-fZChanMin+1;

    
    unsigned int minT = 0;
    unsigned int maxT = 0;
    minT = 0;
    maxT = fNticks;
    unsigned int binT = (maxT-minT);

    for(unsigned int i=0;i<fNofAPA;i++){
      UChMin=fUChanMin + i*fChansPerAPA;
      UChMax=fUChanMax + i*fChansPerAPA;
      VChMin=fVChanMin + i*fChansPerAPA;
      VChMax=fVChanMax + i*fChansPerAPA;
      ZChMin=fZChanMin + i*fChansPerAPA;
      ZChMax=fZChanMax + i*fChansPerAPA;

      // construct the histograms; TH2 constructors: ("Name", "Title", NxBin, xMin, xMax, NyBin, yMax, yMin)
      name.str("");
      name << "fTimeChanU";
      name <<  i;
      title.str("");
      title << "Time vs Channel(Plane U, APA";
      title << i<<")";
      TempHisto = tfs->make<TH2S>(name.str().c_str(),title.str().c_str(), UChMax - UChMin + 1, UChMin, UChMax, binT, minT, maxT);
      fTimeChanU.push_back(TempHisto);

      name.str("");
      name << "fTimeChanV";
      name << i;
      title.str("");
      title << "Time vs Channel(Plane V, APA";
      title << i<<")";
      TempHisto = tfs->make<TH2S>(name.str().c_str(),title.str().c_str(), VChMax - VChMin + 1, VChMin, VChMax, binT, minT, maxT);
      fTimeChanV.push_back(TempHisto);

      name.str("");
      name << "fTimeChanZ";
      name << i;
      title.str("");
      title << "Time vs Channel(Plane Z, APA";
      title <<i<<")";
      TempHisto = tfs->make<TH2S>(name.str().c_str(),title.str().c_str(), ZChMax - ZChMin + 1, ZChMin, ZChMax, binT, minT, maxT);
      fTimeChanZ.push_back(TempHisto);


      fTimeChanU[i]->SetStats(0);
      fTimeChanV[i]->SetStats(0);    
      fTimeChanZ[i]->SetStats(0);    


      fTimeChanU[i]->GetXaxis()->SetTitle("Channel"); fTimeChanU[i]->GetYaxis()->SetTitle("TDC");
      fTimeChanV[i]->GetXaxis()->SetTitle("Channel"); fTimeChanV[i]->GetYaxis()->SetTitle("TDC");
      fTimeChanZ[i]->GetXaxis()->SetTitle("Channel"); fTimeChanZ[i]->GetYaxis()->SetTitle("TDC");
    }


  }

  //-----------------------------------------------------------------------

  void RawEventDisplay::beginRun(const art::Run& run) {
    // place to read databases or run independent info
  }

  //-----------------------------------------------------------------------

  void RawEventDisplay::reconfigure(fhicl::ParameterSet const& p){

    // reconfigure without recompiling
    // read in the parameters from the .fcl file
    // allows for interactive changes of the parameter values

    fTPCInput       = p.get< std::string >("TPCInputModule");
    fTPCInstance    = p.get< std::string >("TPCInstanceName");
    auto const *fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fNticks         = fDetProp->NumberTimeSamples();
    return;
  }



  //-----------------------------------------------------------------------

  void RawEventDisplay::analyze(const art::Event& event) {

    // called once per event

    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    std::cout << "EventNumber = " << fEvent << std::endl;

    // Get the objects holding raw information: RawDigit for TPC data
    art::Handle< std::vector<raw::RawDigit> > RawTPC;
    event.getByLabel(fTPCInput, fTPCInstance, RawTPC);


    // Fill pointer vectors - more useful form for the raw data
    // a more usable form
    std::vector< art::Ptr<raw::RawDigit> > RawDigits;
    art::fill_ptr_vector(RawDigits, RawTPC);



    // Loop over all RawDigits (entire channels)                                                                                                        
    for(auto const & dptr : RawDigits) {
      const raw::RawDigit & digit = *dptr;
      
      // Get the channel number for this digit
      uint32_t chan = digit.Channel();
      // number of samples in uncompressed ADC
      int nSamples = digit.Samples();
      unsigned int apa = std::floor( chan/fChansPerAPA );	  
      int pedestal = (int)digit.GetPedestal();
      
      std::vector<short> uncompressed(nSamples);
      // with pedestal	  
      raw::Uncompress(digit.ADCs(), uncompressed, pedestal, digit.Compression());
      // subtract pedestals
      std::vector<short> uncompPed(nSamples);
      for (int i=0; i<nSamples; i++) uncompPed.at(i)=uncompressed.at(i)-pedestal;
      
      // number of ADC uncompressed without pedestal
      nADC_uncompPed=uncompPed.size();	  
      

      //Induction Plane	  
      if( fGeom->View(chan) == geo::kU){	
				for(unsigned int l=0;l<nADC_uncompPed;l++) {
	  			if(uncompPed.at(l)!=0){
	    			fTimeChanU[apa]->Fill(chan,l, uncompPed.at(l));
	  			}
				}
      }// end of U View
      
      //Induction Plane	  
      if( fGeom->View(chan) == geo::kV){	
				for(unsigned int l=0;l<nADC_uncompPed;l++) {
	  			if(uncompPed.at(l)!=0){
	    			fTimeChanV[apa]->Fill(chan,l, uncompPed.at(l));
	  			}
				}
      }// end of V View

      if ( fGeom->View(chan) == geo::kZ){
				for(unsigned int l=0;l<nADC_uncompPed;l++) {
	  			if(uncompPed.at(l)!=0){
	    			fTimeChanZ[apa]->Fill(chan,l, uncompPed.at(l));
	  			}
				}	
      }

            
    } // RawDigits   
      
    return;
  }
  
}
DEFINE_ART_MODULE(raw_event_display::RawEventDisplay)
  


#endif // RawEventDisplay

