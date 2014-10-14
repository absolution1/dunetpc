//Use this module to create a Tree for raw data display for 35t detector
//Mar. 20, 2014, Seongtae Park

#ifndef RawEVD35tTree_module
#define RawEVD35t_module

// LArSoft includes
#include "Utilities/DetectorProperties.h"
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RawData/raw.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
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

namespace AnalysisExample{

  class RawEVD35tTree : public art::EDAnalyzer{
  public:
 
    explicit RawEVD35tTree(fhicl::ParameterSet const& pset);
    virtual ~RawEVD35tTree();

    void beginJob();

    void beginRun(const art::Run& run);

    void reconfigure(fhicl::ParameterSet const& pset);
 
    void initOutput();

    void analyze(const art::Event& evt); 

  private:

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;

	 //Branch variables for tree
    unsigned int fEvent;
    unsigned int fRun;
    unsigned int fSubRun;

    unsigned int fNUCh;
    unsigned int fNVCh;
    unsigned int fNZ0Ch;
    unsigned int fNZ1Ch;
    // find channel boundaries for each view
    unsigned int fUChanMin;
    unsigned int fUChanMax;
    unsigned int fVChanMin;
    unsigned int fVChanMax;
    unsigned int fZ0ChanMin;
    unsigned int fZ0ChanMax;
    unsigned int fZ1ChanMin;
    unsigned int fZ1ChanMax;
    unsigned int fNticks;
    unsigned int fNofAPA;
    unsigned int fChansPerAPA;
    std::vector<unsigned int> fChan;
    std::vector<unsigned int> fAPA;
    std::vector<unsigned int> fPlane;
    std::vector<std::vector<int> > fADC;
    std::vector<std::vector<int> > fTDC;

    art::ServiceHandle<geo::Geometry> fGeom;
    art::ServiceHandle<util::DetectorProperties> fDetProp;

    TTree *tRD;
    
 }; // class RawEVD35tTree

  //-----------------------------------------------------------------------

  RawEVD35tTree::RawEVD35tTree(fhicl::ParameterSet const& parameterSet)
	: EDAnalyzer(parameterSet)
	{
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------

  // read in the parameters from the .fcl file
  void RawEVD35tTree::reconfigure(fhicl::ParameterSet const& p){
    fRawDigitLabel  =  p.get< std::string >("RawDigitLabel");
    fNticks         = fDetProp->NumberTimeSamples();
    return;
  }


  //-----------------------------------------------------------------------

  RawEVD35tTree::~RawEVD35tTree(){
  }
   
  //-----------------------------------------------------------------------

  void RawEVD35tTree::beginJob(){
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Accquiring geometry data
	 fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2;
    fChansPerAPA = fGeom->Nchannels()/fNofAPA;

	//Defining a Tree
	tRD = tfs->make<TTree>("RawData","Raw Data Display");
	tRD->Branch("Run",&fRun,"Run/i");
	tRD->Branch("SubRun",&fSubRun,"SubRun/i");
	tRD->Branch("Event",&fEvent,"Event/i");
	tRD->Branch("ADC",&fADC);
	tRD->Branch("TDC",&fTDC);
	tRD->Branch("Nticks",&fNticks,"Nticks/i");
	tRD->Branch("NofUChan",&fNUCh,"NofUChan/i");
	tRD->Branch("NofVChan",&fNVCh,"NofVChan/i");
	tRD->Branch("NofZ0Chan",&fNZ0Ch,"NofZ0Chan/i");
	tRD->Branch("NofZ1Chan",&fNZ1Ch,"NofZ1Chan/i");
	tRD->Branch("Chan",&fChan);
	tRD->Branch("APA",&fAPA);
	tRD->Branch("Plane",&fPlane);

    // loop through channels in the first APA to find the channel boundaries for each view
    // will adjust for desired APA after
    fUChanMin = 0;
	 fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2;
    fChansPerAPA = fGeom->Nchannels()/fNofAPA;
    fZ1ChanMax = fChansPerAPA - 1;
    for ( unsigned int c = fUChanMin + 1; c < fZ1ChanMax; c++ ){
      if ( fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU ){
        fVChanMin = c;
        fUChanMax = c - 1;
      }
      if ( fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV ){
        fZ0ChanMin = c;
        fVChanMax = c-1;
      }
      if ( fGeom->View(c) == geo::kZ && fGeom->ChannelToWire(c)[0].TPC == fGeom->ChannelToWire(c-1)[0].TPC + 1 ){
        fZ1ChanMin = c;
        fZ0ChanMax = c-1;
      }
    }

	fNUCh=fUChanMax-fUChanMin+1;
	fNVCh=fVChanMax-fVChanMin+1;
	fNZ0Ch=fZ0ChanMax-fZ0ChanMin+1;
	fNZ1Ch=fZ1ChanMax-fZ1ChanMin+1;

//	ofstream outfile;
//	outfile.open("msglog.txt");
//outfile<<fNUCh<<"  "<<fNVCh<<"  "<<fNZ0Ch<<"  "<<fNZ1Ch<<std::endl;
  }

  //-----------------------------------------------------------------------

  void RawEVD35tTree::beginRun(const art::Run& run){

  }


  //-----------------------------------------------------------------------

  void RawEVD35tTree::analyze( const art::Event& event ){

    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
//	 unsigned int tpcid, cryoid;

    fADC.clear();
    fTDC.clear();
    fChan.clear();
	 fAPA.clear();
	 fPlane.clear();

    // get the objects holding all of the raw data information
    art::Handle< std::vector<raw::RawDigit> > Raw;
    event.getByLabel(fRawDigitLabel, Raw);

    // put it in a more easily usable form
    std::vector< art::Ptr<raw::RawDigit> >  Digits;
    art::fill_ptr_vector(Digits, Raw);

    //loop through all RawDigits (over entire channels)
    for(size_t d = 0; d < Digits.size(); d++){
    art::Ptr<raw::RawDigit> digit;
	 digit=Digits.at(d);

    // get the channel number for this digit
    uint32_t chan = digit->Channel();
	 int nSamples = digit->Samples();
    unsigned int apa = std::floor( chan/fChansPerAPA );
//	 tpcid=fGeom->ChannelToWire(chan)[0].TPC;
//	 cryoid=fGeom->ChannelToWire(chan)[0].Cryostat;
	 int Plane=0;
	 if( fGeom->View(chan) == geo::kU) Plane=0;
	 if( fGeom->View(chan) == geo::kV) Plane=1;
	 if ( fGeom->View(chan) == geo::kZ && fGeom->ChannelToWire(chan)[0].TPC % 2 == 0) Plane=2;
	 if ( fGeom->View(chan) == geo::kZ && fGeom->ChannelToWire(chan)[0].TPC % 2 == 1) Plane=3;

    std::vector<short> uncompressed(digit->Samples());
    raw::Uncompress(digit->fADC, uncompressed, digit->Compression());

    std::vector<int> tmpADC;
    std::vector<int> tmpTDC;
    for (int i=0; i<nSamples; i++) {
		short adc = uncompressed[i];
		if(adc!=0) {
			tmpADC.push_back(int(adc));
			tmpTDC.push_back(int(i));				
			}	
		}
    fADC.push_back(tmpADC);
    fTDC.push_back(tmpTDC);
    fChan.push_back(chan);
    fAPA.push_back(apa);
    fPlane.push_back(Plane);
    } // end RawDigit loop

   tRD->Fill();
   return;
  }

  DEFINE_ART_MODULE(RawEVD35tTree)

} // namespace AnalysisExample

#endif // RawEVD35tTree_module

