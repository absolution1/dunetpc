//Use this module to create histograms for 10kt FD Dual-phase TPC

#ifndef RawEVDDP_Module
#define RawEVDDP_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/raw.h"

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

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

namespace AnalysisExample{

  class RawEVDDP : public art::EDAnalyzer{
  public:
 
    explicit RawEVDDP(fhicl::ParameterSet const& pset);
    virtual ~RawEVDDP();

    void beginJob();

    void beginRun(const art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt); 
    
    
  private:

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;
    //unsigned int fEvent; // unused

    // number of time bins for histogram ranges
    unsigned int fNticks;   
    unsigned int fChPerView; // total number of channels in given view

    // find channel boundaries for each view
    //unsigned int fV0ChanMin; // unused
    //unsigned int fV0ChanMax; // unused
    //unsigned int fV1ChanMin; // unused
    //unsigned int fV1ChanMax; // unused

    //unsigned int fNofAPA; // unused
    //unsigned int fChansPerAPA; // unused
    //unsigned int fUWireMax; // unused
    //unsigned int fVWireMax; // unused
    //unsigned int fZ0WireMax; // unused
    //unsigned int fZ1WireMax; // unused

    

    // art handles
    art::ServiceHandle<geo::Geometry> fGeom; 
    detinfo::DetectorProperties const* fDetProp = nullptr;

    std::vector<TH2I*> fTimeChanU; //data in view 0
    std::vector<TH2I*> fTimeChanV; //data in view 1

    std::vector<TH2I*> fTimeChanThumbU;
    std::vector<TH2I*> fTimeChanThumbV;

    std::vector<TH1I*> fADCMaxDistU;  //max adc per channel per view
    std::vector<TH1I*> fADCMaxDistV;  //max adc per channel per view
    
 }; // class RawEVDDP

  //-----------------------------------------------------------------------

  // read in the parameters from the .fcl file
  RawEVDDP::RawEVDDP(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    this->reconfigure(parameterSet);
  }


  //-----------------------------------------------------------------------

  void RawEVDDP::reconfigure(fhicl::ParameterSet const& p){
    fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fRawDigitLabel  =  p.get< std::string >("RawDigitLabel");
    fNticks         = fDetProp->NumberTimeSamples();
    return;
  }


  //-----------------------------------------------------------------------

  RawEVDDP::~RawEVDDP(){
  }
   
  //-----------------------------------------------------------------------

  void RawEVDDP::beginJob(){
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
    
    //Histogram names and titles
    std::stringstream  name, title;
    
    TH2I* TempHisto2;
    TH1I* TempHisto1;

    //ofstream outfile;
    //outfile.open("msglog.txt");

    
    fChPerView = fGeom->Nchannels()/(fGeom->NTPC() * fGeom->Ncryostats() * fGeom->Nviews());
    mf::LogInfo("RawEVDP")<< "Number of channels per view is "<<fChPerView;
    
    if(fGeom->Nviews() != 2 )
      throw cet::exception("RawEVDDP") << "For DUNE DP expected to have only 2 views ";
    

    //
    unsigned int minT = 0;
    unsigned int maxT = fNticks;
    unsigned int binT = (maxT-minT);

    float minADC = 0.0;    
    float maxADC = 4096;
    int nBins    = (int)((maxADC - minADC)/2.0);

    for(unsigned int i=0;i<fGeom->NTPC();i++)
      {
	
	name.str("");
	name << "fTimeChanU";
	name <<  i;
	title.str("");
	title << "Time vs Channel(Plane U, CRM";
	title << i<<")";
	TempHisto2 = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), fChPerView, 0, fChPerView, binT, minT, maxT);
	fTimeChanU.push_back(TempHisto2);
	
	name.str("");
	name << "fTimeChanThumbU";
	name <<  i;
	title.str("");
	title << "Time vs Channel(Plane U, CRM";
	title << i<<")";
	TempHisto2 = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 32, 0, fChPerView, 32, minT, maxT);
	fTimeChanThumbU.push_back(TempHisto2);

	name.str("");
	name << "fADCMaxU";
	name <<  i;
	title.str("");
	title << "Max ADC per channel (Plane U, CRM";
	title << i<<")";
	TempHisto1 = tfs->make<TH1I>(name.str().c_str(),title.str().c_str(), nBins, minADC, maxADC);
	fADCMaxDistU.push_back(TempHisto1);
	
	//
	name.str("");
	name << "fTimeChanV";
	name <<  i;
	title.str("");
	title << "Time vs Channel(Plane V, CRM";
	title << i<<")";
	TempHisto2 = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), fChPerView, 0, fChPerView, binT, minT, maxT);
	fTimeChanV.push_back(TempHisto2);
	
	name.str("");
	name << "fTimeChanThumbV";
	name <<  i;
	title.str("");
	title << "Time vs Channel(Plane V, CRM";
	title << i<<")";
	TempHisto2 = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 32, 0, fChPerView, 32, minT, maxT);
	fTimeChanThumbV.push_back(TempHisto2);
	    
	name.str("");
	name << "fADCMaxV";
	name <<  i;
	title.str("");
	title << "Max ADC per channel (Plane V, CRM";
	title << i<<")";
	TempHisto1 = tfs->make<TH1I>(name.str().c_str(),title.str().c_str(), nBins, minADC, maxADC);
	fADCMaxDistV.push_back(TempHisto1);
      }
    
  }
  //-----------------------------------------------------------------------

  void RawEVDDP::beginRun(const art::Run& /*run*/){

  }


  //-----------------------------------------------------------------------

  void RawEVDDP::analyze( const art::Event& event )
  {
    
    unsigned int tpcid;//, cryoid;
    //std::stringstream  thumbnameZ0, thumbnameZ1;

    // get the objects holding all of the raw data information
    art::Handle< std::vector<raw::RawDigit> > Raw;
    event.getByLabel(fRawDigitLabel, Raw);

    // put it in a more easily usable form
    std::vector< art::Ptr<raw::RawDigit> >  Digits;
    art::fill_ptr_vector(Digits, Raw);

    //loop through all RawDigits (over entire channels)
    
    for(size_t d = 0; d < Digits.size(); d++)
      {
	art::Ptr<raw::RawDigit> digit;
	digit=Digits.at(d);
	
	// get the channel number for this digit
	uint32_t chan    = digit->Channel();
	tpcid  = fGeom->ChannelToWire(chan)[0].TPC;
	//cryoid = fGeom->ChannelToWire(chan)[0].Cryostat;
	
	// ok this loop does not really work for more than 1 cryostat
	
	std::vector<short> uncompressed(digit->Samples());
	raw::Uncompress(digit->ADCs(), uncompressed, digit->Compression());
	
	short maxadc = 0;
	if( fGeom->View(chan) == geo::kU )
	  {
	    for(unsigned int l=0;l<uncompressed.size();l++) 
	      {
		short adcval = uncompressed.at(l);
		if(uncompressed.at(l) > maxadc)
		  maxadc = adcval;
		
		fTimeChanU[tpcid]->Fill(chan - 2*tpcid*fChPerView, l, adcval);
		if(adcval>0) 
		  fTimeChanThumbU[tpcid]->Fill(chan - 2*tpcid*fChPerView,l, adcval);
	      }
	    if(maxadc>0) fADCMaxDistU[tpcid]->Fill(maxadc);
	  }
	else if( fGeom->View(chan) == geo::kV )
	  {
	    for(unsigned int l=0;l<uncompressed.size();l++) 
	      {
		short adcval = uncompressed.at(l);
		if(uncompressed.at(l) > maxadc)
		  maxadc = adcval;
		
		fTimeChanV[tpcid]->Fill(chan - 2*tpcid*fChPerView-fChPerView, l, adcval);
		if(adcval>0) 
		  fTimeChanThumbV[tpcid]->Fill(chan - 2*tpcid*fChPerView-fChPerView,l, adcval);
	      }
	    if(maxadc > 0) fADCMaxDistV[tpcid]->Fill(maxadc);
	  }
	
      } // end RawDigit loop

    return;
  }

  DEFINE_ART_MODULE(RawEVDDP)

} // namespace AnalysisExample

#endif // RawEVDDP_Module
