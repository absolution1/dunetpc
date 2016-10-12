// Module to create 1D channel occupancy histograms for the DUNE 35t detector
// totally copied from Bill's AnalysisExample_module and Seongtae's Rawevd35t_module
// Matthew Worcester (BNL) 3/7/14

// Added code to check the muon 'track' (as defined by its initial momentum 
// vector, not reconstruction) to see if it hits any of the muon counters.
// MW 5/13/14

#ifndef MuonOccupancy_Module
#define MuonOccupancy_Module

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
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
#include <ctime>

// DUNE code includes
#include "dune/Geometry/MuonCounter35Alg.h"

namespace AnalysisExample{

  class MuonOccupancy : public art::EDAnalyzer{
  
  public:
 
    // constructor and destructor
    explicit MuonOccupancy(fhicl::ParameterSet const& pset);

    virtual ~MuonOccupancy();

    // called once at the beginning of the job
    void beginJob();

    // called once at the beginning of each run
    void beginRun(const art::Run& run);

    // reads in .fcl file parameters, can be called any time
    void reconfigure(fhicl::ParameterSet const& pset);

    // called once per event
    void analyze(const art::Event& evt); 

  private:

    // the parameters we'll read from the .fcl file
    std::string fSimulationProducerLabel;
    std::string fRawDigitLabel;
    int fSelectedPDG; // PDG code for primary particle 

    // log file
    std::ofstream outfile; 

    // define the histograms and ntuples

    // histograms and ntuple for simulated particles
    TH1D* fPositionHist;
    TH1D* fMomentumHist;
    TH1D* fTrackLengthHist;
    TH1D* hit_counter_occupancy;
    TH1D* trigger_occupancy;

    TTree* fSimulationNtuple;

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fPDG;
    int fTrackID;

    double fStartXYZT[4];
    double fEndXYZT[4];
    double fStartPE[4];
    double fEndPE[4];

    // variables and histograms for occupancy plots
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

    std::vector<TH1I*> fChanU;
    std::vector<TH1I*> fChanV;
    std::vector<TH1I*> fChanZ0;
    std::vector<TH1I*> fChanZ1;

    // muon counter geometry
    geo::MuonCounter35Alg *muon_counter;
    int counters_loaded = -1;
    std::vector< std::vector<double> > countergeometry;

    // pointers to services
    art::ServiceHandle<geo::Geometry> fGeom;

 }; // class MuonOccupancy

  //-----------------------------------------------------------------------
  // constructor

  MuonOccupancy::MuonOccupancy(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet){

    // read in the parameters from the .fcl file
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  // destructor
  
  MuonOccupancy::~MuonOccupancy(){}

  //-----------------------------------------------------------------------
  // read in the parameters from the .fcl file

  void MuonOccupancy::reconfigure(fhicl::ParameterSet const& p){

    fSimulationProducerLabel = p.get< std::string >("SimulationLabel");
    fRawDigitLabel = p.get< std::string >("RawDigitLabel");

    fSelectedPDG = p.get< int >("PDGcode");

    return;
  }

  //-----------------------------------------------------------------------
  // executes once at the beginning of the job

  void MuonOccupancy::beginJob(){
 
    // get local time
    time_t rawtime;
    struct tm * timeinfo;

    time (&rawtime);
    timeinfo = localtime (&rawtime);

    // open a basic log file, will overwrite a pre-existing one
    outfile.open("muonoccupancy.log");

    outfile << "MuonOccupancy_module log file, " << asctime(timeinfo) << std::endl;

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // get the detector size to define the histogram ranges
    outfile << "Number of cryostats: " << fGeom->Ncryostats() << std::endl;

    double xrange = 0.; // physical drift axis of cryostat
    double yrange = 0.; // physical height of cryostat
    double zrange = 0.; // physical width of cryostat

    // loop over cryostats
    for(unsigned int c = 0; c < fGeom->Ncryostats(); ++c){

      outfile << "   Cryostat " << c << " [x,y,z] (cm) = ";

      xrange = fGeom->CryostatLength(c);
      outfile << xrange << ", ";

      yrange = fGeom->CryostatHalfHeight(c);
      outfile << yrange << ", ";

      zrange = fGeom->CryostatHalfWidth(c);
      outfile << zrange << std::endl;

    }

    outfile << std::endl;

    // simulated particle histograms
    fPositionHist    = tfs->make<TH1D>("position",";initial position (cm);", int(xrange), 0., xrange*2);
    fMomentumHist    = tfs->make<TH1D>("momentum",";initial momentum (GeV);", 100, 0., 500.);
    fTrackLengthHist = tfs->make<TH1D>("length",  ";particle track length (cm);", 100, 0., 1000.);
    hit_counter_occupancy = tfs->make<TH1D>("Hit counter occupancy",";detector number;",104,0,104);
    trigger_occupancy = tfs->make<TH1D>("Trigger occupancy",";trigger number;",4,1,5);

    // Define our n-tuples, which are limited forms of ROOT
    // TTrees. Start with the TTree itself.
    fSimulationNtuple = tfs->make<TTree>("MuonOccupancySimulation","MuonOccupancySimulation");

    // Define the branches (columns) of our simulation n-tuple. When
    // we write a variable, we give the address of the variable to TTree::Branch.
    fSimulationNtuple->Branch("Event",       &fEvent,          "Event/I");
    fSimulationNtuple->Branch("SubRun",      &fSubRun,         "SubRun/I");
    fSimulationNtuple->Branch("Run",         &fRun,            "Run/I");
    fSimulationNtuple->Branch("TrackID",     &fTrackID,        "TrackID/I");

    // When we write arrays, we give the address of the array to
    // TTree::Branch; in C++ this is simply the array name.
    fSimulationNtuple->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/D");
    fSimulationNtuple->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/D");
    fSimulationNtuple->Branch("StartPE",     fStartPE,         "StartPE[4]/D");
    fSimulationNtuple->Branch("EndPE",       fEndPE,           "EndPE[4]/D");

    // occupancy histogram names and titles
    std::stringstream name, title;

    unsigned int UChMin;
    unsigned int UChMax;
    unsigned int VChMin;
    unsigned int VChMax;
    unsigned int Z0ChMin;
    unsigned int Z0ChMax;
    unsigned int Z1ChMin;
    unsigned int Z1ChMax;

    TH1I* TempHisto;

    // loop through the channels and get the U, V, and Z channel numbers

    // number of APA, assumes that 2 TPC share one APA
    fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2;
    outfile << "Number of APA: " << fNofAPA << std::endl;

    // channels per APA, assumes all channels even split among APA
    fChansPerAPA = fGeom->Nchannels()/fNofAPA;
    outfile << "Number of channels: " << fGeom->Nchannels() << ", " << fChansPerAPA << " per APA" << std::endl;

    // starting and ending points for the loop
    // assumes channel list starts with U and ends with Z
    fUChanMin = 0;
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

    outfile << "U channel number minimum and maximum: " << fUChanMin << "," << fUChanMax 
	    << " (" << fUChanMax-fUChanMin+1 << " total channels)" << std::endl;
    outfile << "V channel number minimum and maximum: " << fVChanMin << "," << fVChanMax
	    << " (" << fVChanMax-fVChanMin+1 << " total channels)" << std::endl;
    outfile << "Z0 channel number minimum and maximum: " << fZ0ChanMin << "," << fZ0ChanMax
	    << " (" << fZ0ChanMax-fZ0ChanMin+1 << " total channels)" << std::endl;
    outfile << "Z1 channel number minimum and maximum: " << fZ1ChanMin << "," << fZ1ChanMax
	    << " (" << fZ1ChanMax-fZ1ChanMin+1 << " total channels)" << std::endl;

    // create the 1D occupancy histograms for each APA
    for(unsigned int i=0;i<fNofAPA;i++){

      UChMin=fUChanMin + i*fChansPerAPA;
      UChMax=fUChanMax + i*fChansPerAPA;
      VChMin=fVChanMin + i*fChansPerAPA;
      VChMax=fVChanMax + i*fChansPerAPA;
      Z0ChMin=fZ0ChanMin + i*fChansPerAPA;
      Z0ChMax=fZ0ChanMax + i*fChansPerAPA;
      Z1ChMin=fZ1ChanMin + i*fChansPerAPA;
      Z1ChMax=fZ1ChanMax + i*fChansPerAPA;

      // construct the histograms; TH1 constructors: ("Name", "Title", NxBin, xMin, xMax)
      name.str("");
      name << "fChanU";
      name <<  i;
      title.str("");
      title << "Hit Channel (Plane U, APA";
      title << i<<")";
      TempHisto = tfs->make<TH1I>(name.str().c_str(),title.str().c_str(), UChMax - UChMin + 2, UChMin-1, UChMax);
      fChanU.push_back(TempHisto);

      name.str("");
      name << "fChanV";
      name << i;
      title.str("");
      title << "Hit Channel (Plane V, APA";
      title << i<<")";
      TempHisto = tfs->make<TH1I>(name.str().c_str(),title.str().c_str(), VChMax - VChMin + 2, VChMin-1, VChMax);
      fChanV.push_back(TempHisto);
      
      name.str("");
      name << "fChanZ0";
      name << i;
      title.str("");
      title << "Hit Channel (Plane Z0, APA";
      title <<i<<")";
      TempHisto = tfs->make<TH1I>(name.str().c_str(),title.str().c_str(), Z0ChMax - Z0ChMin + 2, Z0ChMin-1, Z0ChMax);
      fChanZ0.push_back(TempHisto);

      name.str("");
      name << "fChanZ1";
      name << i;
      title.str("");
      title << "Hit Channel (Plane Z1, APA";
      title << i<<")";
      TempHisto = tfs->make<TH1I>(name.str().c_str(),title.str().c_str(), Z1ChMax - Z1ChMin + 2, Z1ChMin-1, Z1ChMax);
      fChanZ1.push_back(TempHisto);

      fChanU[i]->SetStats(0);    fChanV[i]->SetStats(0);    
      fChanZ0[i]->SetStats(0);    fChanZ1[i]->SetStats(0);

      fChanU[i]->GetXaxis()->SetTitle("Channel number"); fChanU[i]->GetYaxis()->SetTitle("hits");
      fChanV[i]->GetXaxis()->SetTitle("Channel number"); fChanV[i]->GetYaxis()->SetTitle("hits");
      fChanZ0[i]->GetXaxis()->SetTitle("Channel number"); fChanZ0[i]->GetYaxis()->SetTitle("hits");
      fChanZ1[i]->GetXaxis()->SetTitle("Channel number"); fChanZ1[i]->GetYaxis()->SetTitle("hits");

    } // end loop over APA

    outfile << std::endl;

    // load the muon counter positions from a text file

    char counterfile[] = "../Geometry/muoncounters.txt";

    counters_loaded = muon_counter->loadMuonCounterGeometry(counterfile,countergeometry);

    if(!counters_loaded){

      outfile << "ERROR: muon counter geometry failed to load." << std::endl;

    }

    outfile << std::endl;

  }

  //-----------------------------------------------------------------------
  // executes once at the beginning of each run

  void MuonOccupancy::beginRun(const art::Run& /*run*/){
    
  }

  //-----------------------------------------------------------------------
  // executes once per event

  void MuonOccupancy::analyze( const art::Event& event ){

    // start with the basics
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    outfile << "Event " << fEvent << ", run " << fRun << ", subrun " << fSubRun << std::endl;

    // then get the simulated particle information

    // define the handle as an MCParticle vector and fill it with events from the simulation
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    event.getByLabel(fSimulationProducerLabel, particleHandle);

    // define a sorted map in which to put the particles
    std::map< int, const simb::MCParticle* > particleMap;

    // loop over all the particles, find the primary muon, and get the initial/final positions
    for ( auto const& particle : (*particleHandle) )
      {
	// For the methods you can call to get particle information,
	// see $NUTOOLS_INC/SimulationBase/MCParticle.h.
	fTrackID = particle.TrackId();

	// Add the address of the MCParticle to the map, with the track ID as the key.
	particleMap[fTrackID] = &particle; 

	// PDG code of every particle in the event.
	fPDG = particle.PdgCode();

	// locate the primary muon
	if ( particle.Process() == "primary"  &&  fPDG == fSelectedPDG )
	  {
	    outfile << "Primary particle PDG: " << fPDG << std::endl;

	    // A particle has a trajectory, consisting of a set of
	    // 4-positions and 4-mommenta.
	    size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

	    // For trajectories, as for vectors and arrays, the
	    // first point is #0, not #1.
	    int last = numberTrajectoryPoints - 1;
	    const TLorentzVector& positionStart = particle.Position(0);
	    const TLorentzVector& positionEnd   = particle.Position(last);
	    const TLorentzVector& momentumStart = particle.Momentum(0);
	    const TLorentzVector& momentumEnd   = particle.Momentum(last);

	    // move the initial position (for muon counter studies)

	    double x_increment = 0.; // 349.666
	    double z_increment = 0.; // 231.065

	    TVector3 trackStart(positionStart.X()+x_increment,positionStart.Y(),positionStart.Z()+z_increment);

	    //std::cout << trackStart[0] << " " << trackStart[1] << " " << trackStart[2] << std::endl;

	    // fill the histogram of the starting position.
	    fPositionHist->Fill( positionStart.P() );

	    // fill the histogram of the starting momentum.
	    fMomentumHist->Fill( momentumStart.P() );

	    // Fill arrays with the 4-values. (Don't be fooled by
	    // the name of the method; it just puts the numbers from
	    // the 4-vector into the array.)
	    positionStart.GetXYZT( fStartXYZT );
	    positionEnd.GetXYZT( fEndXYZT );
	    momentumStart.GetXYZT( fStartPE );
	    momentumEnd.GetXYZT( fEndPE );

	    outfile << "Initial position: " << fStartXYZT[0] << ", " << fStartXYZT[1] << ", " 
		    << fStartXYZT[2] << " (cm)" << std::endl;
	    outfile << "Initial time: " << fStartXYZT[3] << " (nsec)" << std::endl;

	    outfile << "Initial momentum: " << fStartPE[0] << ", " << fStartPE[1] << ", " 
		    << fStartPE[2] << ", " << fStartPE[3] << " (GeV)" << std::endl;

	    outfile << "Final position: " << fEndXYZT[0] << ", " << fEndXYZT[1] << ", " 
		    << fEndXYZT[2] << " (cm)" << std::endl;
	    outfile << "Final time: " << fEndXYZT[3] << " (nsec)" << std::endl; 

	    outfile << "Final momentum: " << fEndPE[0] << ", " << fEndPE[1] << ", " 
		    << fEndPE[2] << ", " << fEndPE[3] << " (GeV)" << std::endl;

	    // Use a polar-coordinate view of the 4-vectors to
	    // get the track length.
	    double trackLength = ( positionEnd - positionStart ).Rho();

	    // Make a histogram of the track length.
	    fTrackLengthHist->Fill( trackLength );

	    outfile << "Track length: " << trackLength << std::endl;

	    // Check for track intersections with the muon counters.

	    // make sure the muon counter geometry was loaded
	    if(counters_loaded){

	      unsigned int counters_hit = 0;
	      std::vector< std::vector<double> > hitcounters;
    
	      counters_hit = muon_counter->testTrackInAllCounters(fTrackID,
								  trackStart, momentumStart.Vect(), 
								  countergeometry, hitcounters);

	      if(counters_hit != hitcounters.size()){

		outfile << "ERROR: size of hit counters vector is not the same as number of hit counters." 
			<<std::endl;

	      }

	      // condition flags for each layer
	      bool Layer_1_2 = false;
	      bool Layer_3_4_5 = false;
	      bool Layer_E = false;
	      bool Layer_W = false;
	      bool Layer_N_U = false;
	      bool Layer_N_L = false;
	      bool Layer_S_U = false;
	      bool Layer_S_L = false;

	      int Trigger= 0;

	      // loop over the hit counters
	      for(unsigned int hc=0; hc<hitcounters.size(); hc++){

		hit_counter_occupancy->Fill(hitcounters[hc][0]);

		outfile << "Hit counter ID " << hitcounters[hc][0] << ", flag " 
			<< hitcounters[hc][1] << ", track ID " << hitcounters[hc][2] 
			<< ", intersection point: ";

		// loop over the rest of the data for each hit counter
		for(unsigned int nd=3; nd<hitcounters[hc].size(); nd++){

		  outfile << hitcounters[hc][nd] << " ";

		}

		outfile << std::endl;

		// check which layer is hit
		if( 40 <= hitcounters[hc][0] && hitcounters[hc][0] <= 61){
		  Layer_1_2 = true;
		}
		if( hitcounters[hc][0] > 61){
		  Layer_3_4_5 = true;
		}
		if (14 <= hitcounters[hc][0] && hitcounters[hc][0] <=19){
		  Layer_N_U = true;
		}
		if ( 34 <= hitcounters[hc][0] && hitcounters[hc][0] <= 39){
		  Layer_S_L = true;
		}
		if (8 <= hitcounters[hc][0] && hitcounters[hc][0] <= 13){
		  Layer_N_L = true;
		}
		if (28 <= hitcounters[hc][0] && hitcounters[hc][0] <= 33){
		  Layer_S_U = true;
		}
		if (hitcounters[hc][0] <= 7){
		  Layer_E = true;
		}
		if (20 <= hitcounters[hc][0] && hitcounters[hc][0] <= 27){
		  Layer_W = true;
		}
	      }
      
	      // check for a satisfied trigger condition
	      if (Layer_1_2 && Layer_3_4_5){
		Trigger = 1;
	      }
	      if (Layer_N_U && Layer_S_L){
		Trigger = 2;
	      }
	      if (Layer_N_L && Layer_S_U){
		Trigger = 3;
	      }
	      if (Layer_E && Layer_W){
		Trigger = 4;
	      }
	      trigger_occupancy->Fill(Trigger);

	      if(Trigger != 0)
		outfile << "Muon trigger satisfied: " << Trigger << std::endl;

	    } // check that the counter data is loaded

	  } // done with primary muon

      } // end loop over particles

    //outfile << std::endl;

    // now define a handle for the simulated channels
    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel(fSimulationProducerLabel, simChannelHandle);

    // loop over the SimChannel objects in the event.
    for ( auto const& channel : (*simChannelHandle) )
      {

	// Get the numeric ID associated with this channel.
	//auto channelNumber = channel.Channel();

	// count the time slices and energy deposits
	int ntimeslice = 0;
	int nenergydep = 0;

	// Each channel has a map inside it that connects
	// a time slice to energy deposits in the detector
	auto const& timeSlices = channel.TDCIDEMap();
		    
	// For every time slice in this channel:
	for ( auto const& timeSlice : timeSlices )
	  {

	    ntimeslice++;

	    // Each entry in a map is a pair<first,second>.
	    // For the timeSlices map, the 'first' is a time
	    // slice number, which we don't care about in this
	    // example. The 'second' is a vector of IDE
	    // objects.
	    auto const& energyDeposits = timeSlice.second;
			
	    // Loop over the energy deposits. The type of
	    // 'energyDeposit' will be sim::IDE, which is
	    // defined in SimChannel.h.
	    for ( auto const& energyDeposit : energyDeposits )
	      {

		nenergydep++;

		// Check if the track that deposited the
		// energy matches the track of the particle.
		if ( energyDeposit.trackID != fTrackID ){

		  //outfile << "Energy deposit track: " << energyDeposit.trackID
		  //	  << ", track ID: " << fTrackID << std::endl;
		}

	      } // end loop over energy deposits

	  } // end loop over time slices
	/*
	outfile << "simulated channel " << channelNumber << ", type " << fGeom->View(channelNumber)
		<< ", with " << ntimeslice << " time slices and " 
		<< nenergydep << " total energy deposits" << std::endl;
	*/
      } // For each SimChannel

    //outfile << std::endl;

    // fill the simulation ntuple
    fSimulationNtuple->Fill();
	    
    // now get the objects holding all of the raw data information
    art::Handle< std::vector<raw::RawDigit> > Raw;
    event.getByLabel(fRawDigitLabel, Raw);

    // put it in a more easily usable form
    std::vector< art::Ptr<raw::RawDigit> >  Digits;
    art::fill_ptr_vector(Digits, Raw);

    // loop through all RawDigits (over entire channels)
    for(size_t d = 0; d < Digits.size(); d++){

      //outfile << "start digit loop" << std::endl;

      art::Ptr<raw::RawDigit> digit;
      digit=Digits.at(d);

      //outfile << "got digit " << d << std::endl;

      // get the channel number for this digit
      uint32_t chan = digit->Channel();
      //outfile << "got channel " << chan << std::endl;

      //if(chan != 1061){

      unsigned int apa = std::floor( chan/fChansPerAPA );
      //outfile << "got apa " << apa << std::endl;

      // get a vector to hold the uncompressed raw information
      std::vector<short> uncompressed(digit->Samples());

      /*
      // check that the ADC data is not corrupt
      if(uncompressed[0] < 0 || uncompressed[1] < 0 || uncompressed[2] < 0)
	throw cet::exception("MuonOccupancy") << "Bad ADC header, channel " << chan << ", APA " << apa 
					      << ", ADC[0] = " << uncompressed[0] << 

      outfile << d << " (of " << Digits.size() << " digits), channel " << chan 
      	      << " from APA " << apa << std::endl;

      outfile << "compressed length " << uncompressed.size() << ", first/last memory address " 
	      << &uncompressed[0] << "/" << &uncompressed[uncompressed.size()-1] << std::endl;

      outfile << "ADC length " << digit->ADCs().size() << ", compression type " << digit->Compression() << std::endl;

      for(unsigned int i = 0; i < uncompressed.size(); ++i)
	outfile << "ADC vector at " << i << " = " << uncompressed[i] << std::endl;
      */

      raw::Uncompress(digit->ADCs(), uncompressed, digit->Compression());

      //outfile << "uncompressed length " << uncompressed.size() << ", first/last memory address " 
      //	      << &uncompressed[0] << "/" << &uncompressed[uncompressed.size()-1] << std::endl;

      bool hit = false;
      int nsamples = 0;
      
      // check for uncompressed data
      for(unsigned int l=0;l<uncompressed.size();l++) {

	if(uncompressed.at(l)!=0){

	  hit = true;
	  nsamples++;

	} // end check for uncompressed data

      } // end loop over uncompressed digits

      // fill the hits
      if(hit){

	/*
	outfile << "raw hit channel: " << chan << ", samples " << nsamples << ", type " 
		<< fGeom->View(chan) 
		<< " from APA " << apa << ", C2W " << fGeom->ChannelToWire(chan)[0].TPC << std::endl;
	*/

	if( fGeom->View(chan) == geo::kU ){
	  fChanU[apa]->Fill(chan);
	}
	if( fGeom->View(chan) == geo::kV ){
	  fChanV[apa]->Fill(chan);
	}
	if ( fGeom->View(chan) == geo::kZ && fGeom->ChannelToWire(chan)[0].TPC % 2 == 0 ){
	  fChanZ0[apa]->Fill(chan);
	}
	if ( fGeom->View(chan) == geo::kZ && fGeom->ChannelToWire(chan)[0].TPC % 2 == 1 ){
	  fChanZ1[apa]->Fill(chan);
	  //outfile << "filled" << std::endl;
	}

	//outfile << "done filling" << std::endl;

      } // end check for hit

      //outfile << "done looping over digits" << std::endl;

      //}

    } // end RawDigit loop

    outfile << std::endl << "End event " << fEvent << std::endl << std::endl;

    return;
  }

  DEFINE_ART_MODULE(MuonOccupancy)

} // namespace AnalysisExample

#endif // MuonOccupancy_Module
