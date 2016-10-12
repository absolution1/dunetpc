#ifndef DbigEVD4_Module
#define DbigEVD4_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
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

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

namespace AnalysisExample{

  class DbigEVD4apaFD : public art::EDAnalyzer{
  public:
 
    explicit DbigEVD4apaFD(fhicl::ParameterSet const& pset);
    virtual ~DbigEVD4apaFD();

    void beginJob();

    void beginRun(const art::Run& run);

    void reconfigure(fhicl::ParameterSet const& pset);
 
    void analyze(const art::Event& evt); 

  private:

    // the parameters we'll read from the .fcl
    std::string fChanHitLabel;
    std::string fWidHitLabel;
    std::string fHitCheatLabel;
    std::string fClusLabel;
    unsigned int fAPA;
    unsigned int fEvent;


    // diagnostics on other apas to help find other activity
    std::map< unsigned int, unsigned int > fAPAToNumHits;

    // pointers to the time-channel/wire histograms
    std::vector<std::vector<TH2I*> > fAmbigWireU0;
    std::vector<std::vector<TH2I*> > fAmbigWireV0;
    std::vector<std::vector<TH2I*> > fAmbigWireU1;
    std::vector<std::vector<TH2I*> > fAmbigWireV1;
    std::vector<std::vector<TH2I*> > fDisambigWireU0;
    std::vector<std::vector<TH2I*> > fDisambigWireV0;
    std::vector<std::vector<TH2I*> > fDisambigWireU1;
    std::vector<std::vector<TH2I*> > fDisambigWireV1;
    std::vector<std::vector<TH2I*> > fCheatWireU0;
    std::vector<std::vector<TH2I*> > fCheatWireV0;
    std::vector<std::vector<TH2I*> > fCheatWireU1;
    std::vector<std::vector<TH2I*> > fCheatWireV1;
    std::vector<std::vector<TH2I*> > fTimeChanZ0; // fcl option to time reverse for display clarity
    std::vector<std::vector<TH2I*> > fTimeChanZ1;



    // show color-coded channel clusters
    // index in the vector corresponds with root Color_t 0-9
    //std::vector<TH2I*> fClusChanU;
    //std::vector<TH2I*> fClusChanV;
    //std::vector<TH2I*> fClusChanZ0;
    //std::vector<TH2I*> fClusChanZ1;
    // ^^^ removed for now

    // number of time bins for histograms
    unsigned int fNticks;

    // find channel boundaries for each view
    unsigned int fUChMin;
    unsigned int fUChMax;
    unsigned int fVChMin;
    unsigned int fVChMax;
    unsigned int fZ0ChMin;
    unsigned int fZ0ChMax;
    unsigned int fZ1ChMin;
    unsigned int fZ1ChMax;

    unsigned int fnAPAs;
    unsigned int fnEvts;
    unsigned int fChansPerAPA;

    unsigned int fUWireMax;
    unsigned int fVWireMax;
    unsigned int fZ0WireMax;
    unsigned int fZ1WireMax;

    unsigned int fMinT, fMaxT, fMaxTimeRange;
    bool fTimeReverseZ0 = true;

    bool fGotFirstEvtNum = false;
    int  fFirstEvt;

    art::ServiceHandle<geo::Geometry> fGeom;

  }; // class DbigEVD4apaFD

  //-----------------------------------------------------------------------

  // read in the parameters from the .fcl file
  DbigEVD4apaFD::DbigEVD4apaFD(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer( parameterSet )
    , fEvent(-1)
  {
    this->reconfigure(parameterSet);
  }


  //-----------------------------------------------------------------------

  void DbigEVD4apaFD::reconfigure(fhicl::ParameterSet const& p){
    fChanHitLabel   =  p.get< std::string >("ChanHitLabel");
    fWidHitLabel    =  p.get< std::string >("WidHitLabel");
    fHitCheatLabel  =  p.get< std::string >("HitCheatLabel");
    fClusLabel      =  p.get< std::string >("ClusLabel");
    fnEvts          =  p.get< unsigned int >("nEvents");
 
    fTimeReverseZ0  =  p.get< bool         >("TimeReverseZ0"); 

    auto const *fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fNticks         = fDetProp->NumberTimeSamples();
    fnAPAs          = fGeom->NTPC()/2;
    return;
  }


  //-----------------------------------------------------------------------

  DbigEVD4apaFD::~DbigEVD4apaFD(){
  }
   
  //-----------------------------------------------------------------------

  void DbigEVD4apaFD::beginJob(){

    // loop through channels in the first APA to find the channel boundaries for each view
    // will adjust for desired APA after
    fUChMin = 0;
    fChansPerAPA = ( fGeom->Nchannels() )/ (fGeom->NTPC()/2 );
    fZ1ChMax = ( fGeom->Nchannels() )/ (fGeom->NTPC()/2 ) - 1;
    for ( unsigned int c = fUChMin + 1; c < fZ1ChMax; c++ ){
      if ( fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU ){
        fVChMin = c;
        fUChMax = c - 1;
      }
      if ( fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV ){
        fZ0ChMin = c;
        fVChMax = c-1;
      }
      if ( fGeom->View(c) == geo::kZ && fGeom->ChannelToWire(c)[0].TPC == fGeom->ChannelToWire(c-1)[0].TPC + 1 ){
        fZ1ChMin = c;
        fZ0ChMax = c-1;
      }
    }

    //fUWireMax = fGeom->Nwires(0, std::floor( fAPA/2 ), 0);
    //fVWireMax = fGeom->Nwires(1, std::floor( fAPA/2 ), 0);
    //fZ0WireMax = fGeom->Nwires(2, std::floor( fAPA/2 ), 0);
    fUWireMax = fGeom->Nwires(0, 0, 0);
    fVWireMax = fGeom->Nwires(1, 0, 0);
    fZ0WireMax = fGeom->Nwires(2, 0, 0);
    fZ1WireMax = fZ0WireMax;


    unsigned int minT = 0;
    unsigned int maxT = fNticks;

    float ThumbFactor = 0.85;
    unsigned int binT = (maxT-minT)*ThumbFactor;

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;

    int                 minZ0T        = minT;
    int                 maxZ0T        = maxT;
    // makes for more logical display in the root display script
    if(fTimeReverseZ0){ minZ0T        = -fNticks;       
                        maxZ0T        = 0;  }


    std::stringstream hname;
    fAmbigWireU0.resize(fnEvts);       fAmbigWireU1.resize(fnEvts);
    fCheatWireU0.resize(fnEvts);       fCheatWireU1.resize(fnEvts);
    fDisambigWireU0.resize(fnEvts);    fDisambigWireU1.resize(fnEvts);
    fAmbigWireV0.resize(fnEvts);       fAmbigWireV1.resize(fnEvts);
    fCheatWireV0.resize(fnEvts);       fCheatWireV1.resize(fnEvts);
    fDisambigWireV0.resize(fnEvts);    fDisambigWireV1.resize(fnEvts);
    fTimeChanZ0.resize(fnEvts);        fTimeChanZ1.resize(fnEvts);

    // art starts events at number 1
    for(unsigned int evt=1; evt<=fnEvts; evt++){
      // index with evt-1
      fAmbigWireU0[evt-1].resize(fnAPAs);       fAmbigWireU1[evt-1].resize(fnAPAs);
      fCheatWireU0[evt-1].resize(fnAPAs);       fCheatWireU1[evt-1].resize(fnAPAs);
      fDisambigWireU0[evt-1].resize(fnAPAs);    fDisambigWireU1[evt-1].resize(fnAPAs);
      fAmbigWireV0[evt-1].resize(fnAPAs);       fAmbigWireV1[evt-1].resize(fnAPAs);
      fCheatWireV0[evt-1].resize(fnAPAs);       fCheatWireV1[evt-1].resize(fnAPAs);
      fDisambigWireV0[evt-1].resize(fnAPAs);    fDisambigWireV1[evt-1].resize(fnAPAs);
      fTimeChanZ0[evt-1].resize(fnAPAs);        fTimeChanZ1[evt-1].resize(fnAPAs);
      

      for(unsigned int apa=0; apa<fnAPAs;  apa++){
	//unsigned int fUChanMin = fUChMin  + apa*fChansPerAPA;
	//unsigned int fUChanMax = fUChMax  + apa*fChansPerAPA;
	//unsigned int fVChanMin = fVChMin  + apa*fChansPerAPA;
	//unsigned int fVChanMax = fVChMax  + apa*fChansPerAPA;
	unsigned int fZ0ChanMin = fZ0ChMin + apa*fChansPerAPA;
	unsigned int fZ0ChanMax = fZ0ChMax + apa*fChansPerAPA;
	unsigned int fZ1ChanMin = fZ1ChMin + apa*fChansPerAPA;
	unsigned int fZ1ChanMax = fZ1ChMax + apa*fChansPerAPA;
	
	hname.str(""); hname << "fAmbigWireU0_" << apa << "_" << evt;
	fAmbigWireU0[evt-1][apa]    = tfs->make<TH2I>(hname.str().c_str(), "U Ambiguous Hits (side 0)", fUWireMax, 0, fUWireMax, binT, minT, maxT);
	hname.str(""); hname << "fCheatWireU0_" << apa << "_" << evt;
	fCheatWireU0[evt-1][apa]    = tfs->make<TH2I>(hname.str().c_str(), "U Cheated Hits (side 0)", fUWireMax, 0, fUWireMax, binT, minT, maxT);
	hname.str(""); hname << "fDisambigWireU0_" << apa << "_" << evt;
	fDisambigWireU0[evt-1][apa] = tfs->make<TH2I>(hname.str().c_str(), "U Disambiguated Hits (side 0)", fUWireMax, 0, fUWireMax, binT, minT, maxT);
	
	hname.str(""); hname << "fAmbigWireV0_" << apa << "_" << evt;
	fAmbigWireV0[evt-1][apa]    = tfs->make<TH2I>(hname.str().c_str(), "V Ambiguous Hits (side 0)", fVWireMax, 0, fVWireMax, binT, minT, maxT);
	hname.str(""); hname << "fCheatWireV0_" << apa << "_" << evt;
	fCheatWireV0[evt-1][apa]    = tfs->make<TH2I>(hname.str().c_str(), "V Cheated Hits (side 0)", fVWireMax, 0, fVWireMax, binT, minT, maxT);
	hname.str(""); hname << "fDisambigWireV0_" << apa << "_" << evt;
	fDisambigWireV0[evt-1][apa] = tfs->make<TH2I>(hname.str().c_str(), "V Disambiguated Hits (side 0)", fVWireMax, 0, fVWireMax, binT, minT, maxT);
	
	hname.str(""); hname << "fAmbigWireU1_" << apa << "_" << evt;
	fAmbigWireU1[evt-1][apa]    = tfs->make<TH2I>(hname.str().c_str(), "U Ambiguous Hits (side 1)", fUWireMax, 0, fUWireMax, binT, minT, maxT);
	hname.str(""); hname << "fCheatWireU1_" << apa << "_" << evt;
	fCheatWireU1[evt-1][apa]    = tfs->make<TH2I>(hname.str().c_str(), "U Cheated Hits (side 1)", fUWireMax, 0, fUWireMax, binT, minT, maxT);
	hname.str(""); hname << "fDisambigWireU1_" << apa << "_" << evt;
	fDisambigWireU1[evt-1][apa] = tfs->make<TH2I>(hname.str().c_str(), "U Disambiguated Hits (side 1)", fUWireMax, 0, fUWireMax, binT, minT, maxT);
	
	hname.str(""); hname << "fAmbigWireV1_" << apa << "_" << evt;
	fAmbigWireV1[evt-1][apa]    = tfs->make<TH2I>(hname.str().c_str(), "V Ambiguous Hits (side 1)", fVWireMax, 0, fVWireMax, binT, minT, maxT);
	hname.str(""); hname << "fCheatWireV1_" << apa << "_" << evt;
	fCheatWireV1[evt-1][apa]    = tfs->make<TH2I>(hname.str().c_str(), "V Cheated Hits (side 1)", fVWireMax, 0, fVWireMax, binT, minT, maxT);
	hname.str(""); hname << "fDisambigWireV1_" << apa << "_" << evt;
	fDisambigWireV1[evt-1][apa] = tfs->make<TH2I>(hname.str().c_str(), "V Disambiguated Hits (side 1)", fVWireMax, 0, fVWireMax, binT, minT, maxT);
	
      
	fAmbigWireU0[evt-1][apa]->SetStats(0);       fAmbigWireV0[evt-1][apa]->SetStats(0);
	fAmbigWireU1[evt-1][apa]->SetStats(0);       fAmbigWireV1[evt-1][apa]->SetStats(0);
	fDisambigWireU0[evt-1][apa]->SetStats(0);    fDisambigWireV0[evt-1][apa]->SetStats(0);
	fDisambigWireU1[evt-1][apa]->SetStats(0);    fDisambigWireV1[evt-1][apa]->SetStats(0);
	fCheatWireU0[evt-1][apa]->SetStats(0);       fCheatWireV0[evt-1][apa]->SetStats(0);
	fCheatWireU1[evt-1][apa]->SetStats(0);       fCheatWireV1[evt-1][apa]->SetStats(0);
	
	fDisambigWireU0[evt-1][apa]->SetMarkerColor(2);   fDisambigWireV0[evt-1][apa]->SetMarkerColor(2);
	fDisambigWireU1[evt-1][apa]->SetMarkerColor(2);   fDisambigWireV1[evt-1][apa]->SetMarkerColor(2);

	hname.str(""); hname << "fTimeChanZ0_" << apa << "_" << evt;
	fTimeChanZ0[evt-1][apa] = tfs->make<TH2I>(  hname.str().c_str(), "", 
						    (fZ0ChanMax - fZ0ChanMin + 1)*ThumbFactor, 
						    fZ0ChanMin, fZ0ChanMax, 
						    binT, minZ0T, maxZ0T); // NOTE: Time reversed with a fcl trigger
	
	hname.str(""); hname << "fTimeChanZ1_" << apa << "_" << evt;
	fTimeChanZ1[evt-1][apa] = tfs->make<TH2I>( hname.str().c_str(), "", 
						   (fZ1ChanMax - fZ1ChanMin + 1)*ThumbFactor, 
						   fZ1ChanMin, fZ1ChanMax, binT, minT, maxT);

	fTimeChanZ0[evt-1][apa]->SetStats(0);    fTimeChanZ1[evt-1][apa]->SetStats(0);

      } // apas
    } // evts


  }

  //-----------------------------------------------------------------------

  void DbigEVD4apaFD::beginRun(const art::Run& run){

  }


  //-----------------------------------------------------------------------

  void DbigEVD4apaFD::analyze( const art::Event& event ){

    // Get the index for Hist vectors - not 1 if reading a file that 
    // started at a custom firstEvent in the fcl file
    if(!fGotFirstEvtNum){ fFirstEvt = event.id().event(); fGotFirstEvtNum = true; }
    // This works if all evts are contiguously numbered, which if produced by a filter
    // may not be the case. 
    //     unsigned int evtIndex = event.id().event() - fFirstEvt;
    unsigned int evtIndex = ++fEvent;

    art::Handle< std::vector<recob::Hit> > ChHits;
    event.getByLabel(fChanHitLabel, ChHits);
    art::Handle< std::vector<recob::Hit> > WidHits;
    event.getByLabel(fWidHitLabel, WidHits);
    art::Handle< std::vector<recob::Hit> > CheatedHits;
    event.getByLabel(fHitCheatLabel, CheatedHits);


    // plot hits on a given side of an APA vs all 
    // wireIDs corresponding to the it channel
    for( auto const& hit : (*ChHits) ){
      uint32_t     chan   = hit.Channel();
      unsigned int apa    = std::floor( chan/fChansPerAPA );
      double       startT = hit.PeakTimeMinusRMS();
      double       endT   = hit.PeakTimePlusRMS();
      double       charge = hit.PeakAmplitude();
      std::vector< geo::WireID > wids = fGeom->ChannelToWire(chan);

      //if( fAPAToNumHits.count(apa) == 0 ) fAPAToNumHits[apa] = 0;
      //fAPAToNumHits[apa]++;

      for(size_t i=0; i<wids.size(); i++){

	unsigned int wire = wids[i].Wire;

	if( wids[i].TPC % 2 == 0 ){
 
	  if( fGeom->View(chan) == geo::kU ){
	    for( double t = startT; t <= endT; t++ ){
	      fAmbigWireU0[evtIndex][apa]->Fill(wire, t, charge);
	    }
	  }
	  else if( fGeom->View(chan) == geo::kV ){
	    for( double t = startT; t <= endT; t++ ){
	      fAmbigWireV0[evtIndex][apa]->Fill(wire, t, charge);
	    }
	  }
	  else if ( fGeom->View(chan) == geo::kZ ){	
	    for( double t = startT; t <= endT; t++ ){
	      double Z0fillT = t;
	      if(fTimeReverseZ0) Z0fillT = -t; // time reverse for geometrically logical display
	      fTimeChanZ0[evtIndex][apa]->Fill(chan, Z0fillT, charge);
	    }
	  }
	  
	} else if( wids[i].TPC % 2 == 1 ){

	  if( fGeom->View(chan) == geo::kU ){
	    for( double t = startT; t <= endT; t++ ){
	      fAmbigWireU1[evtIndex][apa]->Fill(wire, t, charge);
	    }
	  }
	  else if( fGeom->View(chan) == geo::kV ){
	    for( double t = startT; t <= endT; t++ ){
	      fAmbigWireV1[evtIndex][apa]->Fill(wire, t, charge);
	    }
	  }
	  else if ( hit.View() == geo::kZ ){
	    for( double t = startT; t <= endT; t++ ){
	      fTimeChanZ1[evtIndex][apa]->Fill(chan, t, charge);
	    }
      }

	}

      } // end wid loop
    } // end channel hit per wireID loop



    // for now, let user know where other activity is from here

    //std::map<unsigned int, unsigned int>::iterator apa_it;
    //for( apa_it=fAPAToNumHits.begin(); apa_it!=fAPAToNumHits.end(); apa_it++ )
    //  mf::LogVerbatim("DbigEVD4apaFD") << "APA "  << apa_it->first 
    //			 << " has " << apa_it->second << " hits.";



    // plot hits on a given side of an APA vs. disambiguated wireID
    for( auto const& hit : (*WidHits) ){
      uint32_t     chan   = hit.Channel();
      unsigned int apa    = std::floor( chan/fChansPerAPA );
      unsigned int wire   = hit.WireID().Wire;
      double       startT = hit.PeakTimeMinusRMS();
      double       endT   = hit.PeakTimePlusRMS();
      double       charge = hit.PeakAmplitude();
 
      if(hit.WireID().TPC % 2 == 0){

	if( fGeom->View(chan) == geo::kU ){
	  for( double t = startT; t <= endT; t++ ){
	    fDisambigWireU0[evtIndex][apa]->Fill(wire, t, charge);
	  }
	}
	if( fGeom->View(chan) == geo::kV ){
	  for( double t = startT; t <= endT; t++ ){
	    fDisambigWireV0[evtIndex][apa]->Fill(wire, t, charge);
	  }
	}
	
      } else if(hit.WireID().TPC % 2 == 1){
	
	if( fGeom->View(chan) == geo::kU ){
	  for( double t = startT; t <= endT; t++ ){
	    fDisambigWireU1[evtIndex][apa]->Fill(wire, t, charge);
	  }
	}
	if( fGeom->View(chan) == geo::kV ){
	  for( double t = startT; t <= endT; t++ ){
	    fDisambigWireV1[evtIndex][apa]->Fill(wire, t, charge);
	  }
	}
	
      }

    } // end wireID hit per wireID loop


    // plot cheated hits on a given side of an APA vs. cheated wireID
    for( auto const& hit : (*CheatedHits) ){
      uint32_t     chan   = hit.Channel();
      unsigned int apa    = std::floor( chan/fChansPerAPA );
      unsigned int wire   = hit.WireID().Wire;
      double       startT = hit.PeakTimeMinusRMS();
      double       endT   = hit.PeakTimePlusRMS();
      double       charge = hit.PeakAmplitude();
 
      if(hit.WireID().TPC % 2 == 0){

	if( fGeom->View(chan) == geo::kU ){
	  for( double t = startT; t <= endT; t++ ){
	    fCheatWireU0[evtIndex][apa]->Fill(wire, t, charge);
	  }
	}
	if( fGeom->View(chan) == geo::kV ){
	  for( double t = startT; t <= endT; t++ ){
	    fCheatWireV0[evtIndex][apa]->Fill(wire, t, charge);
	  }
	}
	
      } else if(hit.WireID().TPC % 2 == 1){
	
	if( fGeom->View(chan) == geo::kU ){
	  for( double t = startT; t <= endT; t++ ){
	    fCheatWireU1[evtIndex][apa]->Fill(wire, t, charge);
	  }
	}
	if( fGeom->View(chan) == geo::kV ){
	  for( double t = startT; t <= endT; t++ ){
	    fCheatWireV1[evtIndex][apa]->Fill(wire, t, charge);
	  }
	}
	
      }

    } // end cheated hit loop
    

    return;
  }

  DEFINE_ART_MODULE(DbigEVD4apaFD)

} // namespace AnalysisExample

#endif // DbigEVD4_Module
