//Use this module to create histograms for 35t FD TPCs
//The output file can be used to display time vs. channel hit information 
//for U,V,Z planes using 'rawEVD35t.c' root script.
//Sep. 25, 2013, Seongtae Park
//Nov. 4, 2013, Thumbnail histograms are included in the root file

#ifndef RawEVD35t_Module
#define RawEVD35t_Module

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

namespace AnalysisExample{

  class RawEVD35t : public art::EDAnalyzer{
  public:
 
    explicit RawEVD35t(fhicl::ParameterSet const& pset);
    virtual ~RawEVD35t();

    void beginJob();

    void beginRun(const art::Run& run);

    void reconfigure(fhicl::ParameterSet const& pset);
 
    void analyze(const art::Event& evt); 

	unsigned int getAPAindex(unsigned int apa) {
	 for (int k=0; k<4; k++){if(apa==UVPlane[k]) {return k; }}
	 return -1;
	}

	unsigned int getTPCindex(unsigned int tpc) {
	 for (int k=0; k<8; k++){if(tpc==ZPlane[k]) {return k; }}
	 return -1;
	}

  private:

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;
    //unsigned int fEvent; // unused

    // number of time bins for histogram ranges
    unsigned int fNticks;
    bool fUncompressWithPed;

    // find channel boundaries for each view
    unsigned int fUChanMin;
    unsigned int fUChanMax;
    unsigned int fVChanMin;
    unsigned int fVChanMax;
    unsigned int fZ0ChanMin;
    unsigned int fZ0ChanMax;
    unsigned int fZ1ChanMin;
    unsigned int fZ1ChanMax;

    unsigned int fNofAPA;
    unsigned int fChansPerAPA;
    //unsigned int fUWireMax; // unused
    //unsigned int fVWireMax; // unused
    //unsigned int fZ0WireMax; // unused
    //unsigned int fZ1WireMax; // unused

    //unsigned int fMinT, fMaxT, fMaxTimeRange; // unused


    art::ServiceHandle<geo::Geometry> fGeom;

    std::vector<TH2I*> fTimeChanU;
    std::vector<TH2I*> fTimeChanV;
    std::vector<TH2I*> fTimeChanZ0;
    std::vector<TH2I*> fTimeChanZ1;

    std::vector<TH2I*> fTimeChanThumbU;
    std::vector<TH2I*> fTimeChanThumbV;
    std::vector<TH2I*> fTimeChanThumbZ0;
    std::vector<TH2I*> fTimeChanThumbZ1;

    TH2I* fChargeSumU;
    TH2I* fChargeSumV;
    TH2I* fChargeSumZ;

	 unsigned int UVPlane[4]={3,2,1,0};
	 unsigned int ZPlane[8]={7,6,5,4,3,2,1,0};

 }; // class RawEVD35t

  //-----------------------------------------------------------------------

  // read in the parameters from the .fcl file
  RawEVD35t::RawEVD35t(fhicl::ParameterSet const& parameterSet)
	: EDAnalyzer(parameterSet)
	{
    this->reconfigure(parameterSet);
  }


  //-----------------------------------------------------------------------

  void RawEVD35t::reconfigure(fhicl::ParameterSet const& p){
    fRawDigitLabel  =  p.get< std::string >("RawDigitLabel");
    // auto const* fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    //    fNticks         = fDetProp->NumberTimeSamples();
    fNticks  =  (unsigned int) p.get< int >("TicksToDraw");
    fUncompressWithPed  = p.get< bool         >("UncompressWithPed", true);
    return;
  }


  //-----------------------------------------------------------------------

  RawEVD35t::~RawEVD35t(){
  }
   
  //-----------------------------------------------------------------------

  void RawEVD35t::beginJob(){
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
       art::ServiceHandle<art::TFileService> tfs;

	//Histogram names and titles
	std::stringstream  name, title;

    unsigned int UChMin;
    unsigned int UChMax;
    unsigned int VChMin;
    unsigned int VChMax;
    unsigned int Z0ChMin;
    unsigned int Z0ChMax;
    unsigned int Z1ChMin;
    unsigned int Z1ChMax;
    TH2I* TempHisto;

    std::ofstream outfile;
    outfile.open("msglog.txt");

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

outfile<<fChansPerAPA<<"  "<<fGeom->Ncryostats()<<"  "<<fNofAPA<<std::endl;

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
    Z0ChMin=fZ0ChanMin + i*fChansPerAPA;
    Z0ChMax=fZ0ChanMax + i*fChansPerAPA;
    Z1ChMin=fZ1ChanMin + i*fChansPerAPA;
    Z1ChMax=fZ1ChanMax + i*fChansPerAPA;

    // construct the histograms; TH2 constructors: ("Name", "Title", NxBin, xMin, xMax, NyBin, yMax, yMin)
    name.str("");
    name << "fTimeChanU";
    name <<  i;
    title.str("");
    title << "Time vs Channel(Plane U, APA";
    title << i<<")";
	 TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), UChMax - UChMin + 1, UChMin, UChMax, binT, minT, maxT);
	 fTimeChanU.push_back(TempHisto);

    name.str("");
    name << "fTimeChanThumbU";
    name <<  i;
    title.str("");
    title << "Time vs Channel(Plane U, APA";
    title << i<<")";
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 138, UChMin, UChMax, 138, minT, maxT);
	 fTimeChanThumbU.push_back(TempHisto);

    name.str("");
    name << "fTimeChanV";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane V, APA";
    title << i<<")";
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), VChMax - VChMin + 1, VChMin, VChMax, binT, minT, maxT);
	 fTimeChanV.push_back(TempHisto);

    name.str("");
    name << "fTimeChanThumbV";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane V, APA";
    title << i<<")";
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 138, VChMin, VChMax, 138, minT, maxT);
	 fTimeChanThumbV.push_back(TempHisto);

    name.str("");
    name << "fTimeChanZ0";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane Z0, APA";
    title <<i<<")";
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), Z0ChMax - Z0ChMin + 1, Z0ChMin, Z0ChMax, binT, minT, maxT);
	 fTimeChanZ0.push_back(TempHisto);

    name.str("");
    name << "fTimeChanThumbZ0";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane Z0, APA";
    title <<i<<")";
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 111, Z0ChMin, Z0ChMax, 111, minT, maxT);
	 fTimeChanThumbZ0.push_back(TempHisto);

    name.str("");
    name << "fTimeChanZ1";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane Z1, APA";
    title << i<<")";
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), Z1ChMax - Z1ChMin + 1, Z1ChMin, Z1ChMax, binT, minT, maxT);
	 fTimeChanZ1.push_back(TempHisto);

    name.str("");
    name << "fTimeChanThumbZ1";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane Z1, APA";
    title << i<<")";
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 111, Z1ChMin, Z1ChMax, 111, minT, maxT);
	 fTimeChanThumbZ1.push_back(TempHisto);

    fTimeChanU[i]->SetStats(0);    fTimeChanV[i]->SetStats(0);    
	 fTimeChanZ0[i]->SetStats(0);    fTimeChanZ1[i]->SetStats(0);
    fTimeChanThumbU[i]->SetStats(0);    fTimeChanThumbV[i]->SetStats(0);    
	 fTimeChanThumbZ0[i]->SetStats(0);    fTimeChanThumbZ1[i]->SetStats(0);

    fTimeChanU[i]->GetXaxis()->SetTitle("Channel"); fTimeChanU[i]->GetYaxis()->SetTitle("TDC");
    fTimeChanV[i]->GetXaxis()->SetTitle("Channel"); fTimeChanV[i]->GetYaxis()->SetTitle("TDC");
    fTimeChanZ0[i]->GetXaxis()->SetTitle("Channel"); fTimeChanZ0[i]->GetYaxis()->SetTitle("TDC");
    fTimeChanZ1[i]->GetXaxis()->SetTitle("Channel"); fTimeChanZ1[i]->GetYaxis()->SetTitle("TDC");
	}

    fChargeSumU= tfs->make<TH2I>("hChargeSumU","Charge Sum on U-planes", 2,0.5,2.5,2,0.5,2.5 );
    fChargeSumV= tfs->make<TH2I>("hChargeSumV","Charge Sum on V-planes", 2,0.5,2.5,2,0.5,2.5 );
    fChargeSumZ= tfs->make<TH2I>("hChargeSumZ","Charge Sum on Z-planes", 4,0.5,4.5,4,0.5,4.5);
    fChargeSumU->SetStats(0);    fChargeSumV->SetStats(0);    fChargeSumZ->SetStats(0);  

  }

  //-----------------------------------------------------------------------

  void RawEVD35t::beginRun(const art::Run& /*run*/){

  }


  //-----------------------------------------------------------------------

  void RawEVD35t::analyze( const art::Event& event ){

    unsigned int tpcid, cryoid;
    std::stringstream  thumbnameZ0, thumbnameZ1;
    
    // get the objects holding all of the raw data information
    art::Handle< std::vector<raw::RawDigit> > Raw;
    std::cout << "raw digit label check: " << fRawDigitLabel << std::endl;
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
      unsigned int apa = std::floor( chan/fChansPerAPA );
      tpcid=fGeom->ChannelToWire(chan)[0].TPC;
      cryoid=fGeom->ChannelToWire(chan)[0].Cryostat;
      
      int nSamples = digit->Samples();
      std::vector<short> uncompressed(nSamples);
      //    raw::Uncompress(digit->ADCs(), uncompressed, digit->Compression());
      int pedestal = (int)digit->GetPedestal();
      //      std::cout << "channel " << chan << " pedestal " << pedestal << std::endl;
      // uncompress the data
      if (fUncompressWithPed){
	raw::Uncompress(digit->ADCs(), uncompressed, pedestal, digit->Compression());
      }
      else{
	raw::Uncompress(digit->ADCs(), uncompressed, digit->Compression());
      }
      
      // subtract pedestals
	std::vector<short> ladc(nSamples);
	      for (int i=0; i<nSamples; i++) ladc[i]=uncompressed[i]-pedestal;
	// for (int i=0; i<nSamples; i++) {ladc[i]=uncompressed[i]-pedestal;
	//   if (i<10) std::cout << uncompressed[i] << " " << ladc[i] << std::endl;
	      //	}
      if( fGeom->View(chan) == geo::kU ){
	for(unsigned int l=0;l<ladc.size();l++) {
	  if(ladc.at(l)!=0){
	    fTimeChanU[apa]->Fill(chan,l, ladc.at(l));
	    if(ladc.at(l)>0) fTimeChanThumbU[apa]->Fill(chan,l, ladc.at(l));
	    fChargeSumU->Fill(1+getAPAindex(apa)%2, 2-(getAPAindex(apa)/2),std::abs(ladc.at(l))/2);
	  }
	}
      }
      if( fGeom->View(chan) == geo::kV ){
	for(unsigned int l=0;l<ladc.size();l++) {
	     if(ladc.at(l)!=0){
	       fTimeChanV[apa]->Fill(chan,l, ladc.at(l));
	       if(ladc.at(l)>0) fTimeChanThumbV[apa]->Fill(chan,l, ladc.at(l));
	       fChargeSumV->Fill(1+getAPAindex(apa)%2, 2-(getAPAindex(apa)/2),std::abs(ladc.at(l))/2);
	     }
	   }
      }
      if ( fGeom->View(chan) == geo::kZ && fGeom->ChannelToWire(chan)[0].TPC % 2 == 0 ){
	for(unsigned int l=0;l<ladc.size();l++) {
	  if(ladc.at(l)!=0){
	    fTimeChanZ0[apa]->Fill(chan,l, ladc.at(l));
	       if(ladc.at(l)>0) fTimeChanThumbZ0[apa]->Fill(chan,l, ladc.at(l));
	       fChargeSumZ->Fill(1+getTPCindex(cryoid*4+tpcid)%4, 4-(getTPCindex(cryoid*4+tpcid)/4),ladc.at(l));
	  }
	}
      }
      if ( fGeom->View(chan) == geo::kZ && fGeom->ChannelToWire(chan)[0].TPC % 2 == 1 ){
	for(unsigned int l=0;l<ladc.size();l++) {
	  if(ladc.at(l)!=0){
	    fTimeChanZ1[apa]->Fill(chan,l, ladc.at(l));
	    if(ladc.at(l)>0) fTimeChanThumbZ1[apa]->Fill(chan,l, ladc.at(l));
	    fChargeSumZ->Fill(1+getTPCindex(cryoid*4+tpcid)%4, 4-(getTPCindex(cryoid*4+tpcid)/4),ladc.at(l));
	  }
	}
      }
      
    } // end RawDigit loop

    return;
  }

  DEFINE_ART_MODULE(RawEVD35t)

} // namespace AnalysisExample

#endif // RawEVD35t_Module
