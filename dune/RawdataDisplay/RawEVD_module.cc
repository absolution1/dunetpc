//Use this module to create histograms for 10kt FD TPCs
//The output file can be used to display time vs. channel hit information 
//for U,V,Z planes using 'rawEVD.c' root script.
//Sep. 25, 2013, Seongtae Park
//Mov. 4, 2013, Thumbnail histograms are included in the root file

#ifndef RawEVD_Module
#define RawEVD_Module

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

  class RawEVD : public art::EDAnalyzer{
  public:
 
    explicit RawEVD(fhicl::ParameterSet const& pset);
    virtual ~RawEVD();

    void beginJob();

    void beginRun(const art::Run& run);

    void reconfigure(fhicl::ParameterSet const& pset);
 
    void analyze(const art::Event& evt); 

	unsigned int getAPAindex(unsigned int apa) {
	 for (int k=0; k<120; k++){if(apa==UVPlane[k]) {return k; }}
	 return -1;
	}

	unsigned int getTPCindex(unsigned int tpc) {
	 for (int k=0; k<240; k++){if(tpc==ZPlane[k]) {return k; }}
	 return -1;
	}

  private:

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;
    //unsigned int fEvent; // unused

    // number of time bins for histogram ranges
    unsigned int fNticks;

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

	unsigned int UVPlane[120]={
	119,118,117,59,58,57,
	113,112,111,53,52,51,
	107,106,105,47,46,45,
	101,100,99,41,40,39,
	95,94,93,35,34,33,
	89,88,87,29,28,27,
	83,82,81,23,22,21,
	77,76,75,17,16,15,
	71,70,69,11,10,9,
	65,64,63,5,4,3,
	116,115,114,56,55,54,
	110,109,108,50,49,48,
	104,103,102,44,43,42,
	98,97,96,38,37,36,
	92,91,90,32,31,30,
	86,85,84,26,25,24,
	80,79,78,20,19,18,
	74,73,72,14,13,12,
	68,67,66,8,7,6,
	62,61,60,2,1,0};

	unsigned int ZPlane[240]={
	239,238,237,236,235,234,119,118,117,116,115,114,
	227,226,225,224,223,222,107,106,105,104,103,102,
	215,214,213,212,211,210,95,94,93,92,91,90,
	203,202,201,200,199,198,83,82,81,80,79,78,
	191,190,189,188,187,186,71,70,69,68,67,66,
	179,178,177,176,175,174,59,58,57,56,55,54,
	167,166,165,164,163,162,47,46,45,44,43,42,
	155,154,153,152,151,150,35,34,33,32,31,30,
	143,142,141,140,139,138,23,22,21,20,19,18,
	131,130,129,128,127,126,11,10,9,8,7,6,
	233,232,231,230,229,228,113,112,111,110,109,108,
	221,220,219,218,217,216,101,100,99,98,97,96,
	209,208,207,206,205,204,89,88,87,86,85,84,
	197,196,195,194,193,192,77,76,75,74,73,72,
	185,184,183,182,181,180,65,64,63,62,61,60,
	173,172,171,170,169,168,53,52,51,50,49,48,
	161,160,159,158,157,156,41,40,39,38,37,36,
	149,148,147,146,145,144,29,28,27,26,25,24,
	137,136,135,134,133,132,17,16,15,14,13,12,
	125,124,123,122,121,120,5,4,3,2,1,0
	};

 }; // class RawEVD

  //-----------------------------------------------------------------------

  // read in the parameters from the .fcl file
  RawEVD::RawEVD(fhicl::ParameterSet const& parameterSet)
	: EDAnalyzer(parameterSet)
	{
    this->reconfigure(parameterSet);
  }


  //-----------------------------------------------------------------------

  void RawEVD::reconfigure(fhicl::ParameterSet const& p){
    fRawDigitLabel  =  p.get< std::string >("RawDigitLabel");
    auto const *fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fNticks         = fDetProp->NumberTimeSamples();
    return;
  }


  //-----------------------------------------------------------------------

  RawEVD::~RawEVD(){
  }
   
  //-----------------------------------------------------------------------

  void RawEVD::beginJob(){
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
	 fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2; //NTPC:120 for each Cryostat
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
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 32, UChMin, UChMax, 32, minT, maxT);
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
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 32, VChMin, VChMax, 32, minT, maxT);
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
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 32, Z0ChMin, Z0ChMax, 32, minT, maxT);
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
    TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), 32, Z1ChMin, Z1ChMax, 32, minT, maxT);
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

    fChargeSumU= tfs->make<TH2I>("hChargeSumU","Charge Sum on U-planes", 6,0.5,6.5,20,0.5,20.5 );
    fChargeSumV= tfs->make<TH2I>("hChargeSumV","Charge Sum on V-planes", 6,0.5,6.5,20,0.5,20.5 );
    fChargeSumZ= tfs->make<TH2I>("hChargeSumZ","Charge Sum on Z-planes", 12,0.5,12.5,20,0.5,20.5);
    fChargeSumU->SetStats(0);    fChargeSumV->SetStats(0);    fChargeSumZ->SetStats(0);  

  }

  //-----------------------------------------------------------------------

  void RawEVD::beginRun(const art::Run& /*run*/){

  }


  //-----------------------------------------------------------------------

  void RawEVD::analyze( const art::Event& event ){

	unsigned int tpcid, cryoid;
	std::stringstream  thumbnameZ0, thumbnameZ1;

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
    unsigned int apa = std::floor( chan/fChansPerAPA );
	 tpcid=fGeom->ChannelToWire(chan)[0].TPC;
	 cryoid=fGeom->ChannelToWire(chan)[0].Cryostat;

    std::vector<short> uncompressed(digit->Samples());
    raw::Uncompress(digit->ADCs(), uncompressed, digit->Compression());

	if( fGeom->View(chan) == geo::kU ){
		for(unsigned int l=0;l<uncompressed.size();l++) {
			if(uncompressed.at(l)!=0){
				fTimeChanU[apa]->Fill(chan,l, uncompressed.at(l));
				if(uncompressed.at(l)>0) fTimeChanThumbU[apa]->Fill(chan,l, uncompressed.at(l));
				fChargeSumU->Fill(1+getAPAindex(apa)%6, 20-(getAPAindex(apa)/6),std::abs(uncompressed.at(l))/2);
				}
			}
		}
	if( fGeom->View(chan) == geo::kV ){
		for(unsigned int l=0;l<uncompressed.size();l++) {
			if(uncompressed.at(l)!=0){
				fTimeChanV[apa]->Fill(chan,l, uncompressed.at(l));
				if(uncompressed.at(l)>0) fTimeChanThumbV[apa]->Fill(chan,l, uncompressed.at(l));
				fChargeSumV->Fill(1+getAPAindex(apa)%6, 20-(getAPAindex(apa)/6),std::abs(uncompressed.at(l))/2);
				}
			}
		}
	if ( fGeom->View(chan) == geo::kZ && fGeom->ChannelToWire(chan)[0].TPC % 2 == 0 ){
		for(unsigned int l=0;l<uncompressed.size();l++) {
			if(uncompressed.at(l)!=0){
				fTimeChanZ0[apa]->Fill(chan,l, uncompressed.at(l));
				if(uncompressed.at(l)>0) fTimeChanThumbZ0[apa]->Fill(chan,l, uncompressed.at(l));
				fChargeSumZ->Fill(1+getTPCindex(cryoid*120+tpcid)%12, 20-(getTPCindex(cryoid*120+tpcid)/12),uncompressed.at(l));
				}
			}
		}
	if ( fGeom->View(chan) == geo::kZ && fGeom->ChannelToWire(chan)[0].TPC % 2 == 1 ){
		for(unsigned int l=0;l<uncompressed.size();l++) {
			if(uncompressed.at(l)!=0){
				fTimeChanZ1[apa]->Fill(chan,l, uncompressed.at(l));
				if(uncompressed.at(l)>0) fTimeChanThumbZ1[apa]->Fill(chan,l, uncompressed.at(l));
				fChargeSumZ->Fill(1+getTPCindex(cryoid*120+tpcid)%12, 20-(getTPCindex(cryoid*120+tpcid)/12),uncompressed.at(l));
				}
			}
		}

    } // end RawDigit loop

    return;
  }

  DEFINE_ART_MODULE(RawEVD)

} // namespace AnalysisExample

#endif // RawEVD_Module
