//////////////////////////////////////////////////////////////////////// 
// Class:       GapFilter                                             //
// Module Type: Analyzer                                              // 
// File:        GapFilter_module.cc                                   //
// Author:      Tristan Blackburn - t.blackburn@sussex.ac.uk          //
// Description: Module takes RawDigit input, produces hits, and uses  //
//              hit information to filter gap crossing events.        //
////////////////////////////////////////////////////////////////////////
#ifndef GapFilter_Module
#define GapFilter_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Principal/Handle.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h" 
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/Simulation/SimChannel.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>


const int kMaxHits       = 10000; //maximum number of hits

namespace MyGapFilter {

  class GapFilter : public art::EDAnalyzer {
  public:
    explicit GapFilter(fhicl::ParameterSet const & p);
    virtual ~GapFilter();

    void beginJob();

    void reconfigure(fhicl::ParameterSet const& pset);

    void analyze (const art::Event& evt); 

    TH2F* fPrimaryMom_vs_TrackLength;

  private:

    void ResetVars();
  
    TTree* fTree;

    //run information
    int    run;
    int    subrun;
    int    event;
    double evttime;
    double efield[3];
    int    t0;
    
    //Comparison
    int    nhits;
    int    hit_cryostat[kMaxHits];
    int    hit_tpc[kMaxHits];
    int    hit_plane[kMaxHits];
    int    hit_wire[kMaxHits];
    int    hit_channel[kMaxHits];
    double hit_peakT[kMaxHits];
    double hit_charge[kMaxHits];
    double hit_ph[kMaxHits];

    //Basic crossings
    std::vector<int> gap1;
    std::vector<int> gap2;
    std::vector<int> gap3;
    std::vector<int> gap4;
    std::vector<int> gap5;

    //Horizontal crossings
    std::vector<int> cross12;
    std::vector<int> cross34;

    //Diagonal crossings
    std::vector<int> cross14;
    std::vector<int> cross23;

    //Inclusive crossings
    std::vector<int> cross1and2or4;
    std::vector<int> cross2and1or3;
    std::vector<int> cross1or3and2or4;

    //Bools for deciding on gap crossers
    bool TPC1     = false; 
    bool TPC7     = false;

    bool TPC5GAP1 = false;
    bool GAP1     = false;

    bool TPC5GAP2 = false;
    bool GAP2     = false;

    bool TPC3GAP3 = false;
    bool GAP3     = false;

    bool TPC3GAP4 = false;
    bool GAP4     = false;

    bool TPC3     = false;
    bool TPC5     = false;
    bool GAP5     = false;
    
    //FCL labels
    std::string fHitsModuleLabel;
    std::string fSimulationProducerLabel;

  };


  GapFilter::GapFilter(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    // Read in the parameters from the .fcl file.
    this->reconfigure(pset);
  }

  void GapFilter::reconfigure(fhicl::ParameterSet const& p)
  {
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.

    fHitsModuleLabel         = p.get< std::string >("HitsModuleLabel");
    fSimulationProducerLabel = p.get< std::string >("SimulationProducerLabel");
    return;
  }

  GapFilter::~GapFilter()
  {
    // Clean up dynamic memory and other resources here.
  }

  void GapFilter::analyze(const art::Event& evt)
  {
    // Implementation of required member function here.
    ResetVars();

    art::ServiceHandle<geo::Geometry> geom;
    
    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();

    //Reset vectors per event
    gap1.clear();
    gap2.clear();
    gap3.clear();
    gap4.clear();
    gap5.clear();
    cross12.clear();
    cross34.clear();
    cross14.clear();
    cross23.clear();
    cross1and2or4.clear();
    cross2and1or3.clear();
    cross1or3and2or4.clear();

    //Art handles

    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    evt.getByLabel(fSimulationProducerLabel, particleHandle);
    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    evt.getByLabel(fSimulationProducerLabel, simChannelHandle);
    std::map< int, const simb::MCParticle* > particleMap;

    // Get hit list
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
    nhits = hitlist.size();
  
    //loop to determine hits and assign bool metrics for deciding on crossers
    for (size_t k = 0; k < hitlist.size(); ++k){	  
      unsigned int channel = hitlist[k]->Channel();
      geo::WireID wireid   = hitlist[k]->WireID();
      hit_cryostat[k]      = wireid.Cryostat;
      hit_tpc[k]           = wireid.TPC;
      hit_plane[k]         = wireid.Plane;
      hit_wire[k]          = wireid.Wire;
      hit_channel[k]       = channel;
      hit_peakT[k]         = hitlist[k]->PeakTime();
      hit_charge[k]        = hitlist[k]->Integral();
      hit_ph[k]            = hitlist[k]->PeakAmplitude();
   
      if (hit_channel[k] == 511  && hit_charge[k] != 0) TPC1                           = true;
      if (hit_channel[k] == 1424 && hit_charge[k] != 0) TPC5GAP1                       = true;
      if (hit_channel[k] == 1535 && hit_charge[k] != 0) TPC5GAP2                       = true;
      if (hit_channel[k] == 912  && hit_charge[k] != 0) TPC3GAP3                       = true;
      if (hit_channel[k] == 1023 && hit_charge[k] != 0) TPC3GAP4                       = true;
      if (hit_channel[k] == 1936 && hit_charge[k] != 0) TPC7                           = true;
      if (hit_channel[k] >= 912  && hit_channel[k] <= 1023 && hit_charge[k] != 0) TPC3 = true;  
      if (hit_channel[k] >= 1424 && hit_channel[k] <= 1535 && hit_charge[k] != 0) TPC5 = true;  
    }  

    //Use bools to decide whether an event crosses a gap
    
    if (TPC1 == true && TPC5GAP1 == true) {
      GAP1 = true;
      gap1.push_back(event);
    }
    else GAP1 = false;
	
    if (TPC5GAP2 == true && TPC7 == true) {
      GAP2 = true;
      gap2.push_back(event);
    }
    else GAP2 = false;
	
    if (TPC1 == true && TPC3GAP3 == true) {
      GAP3 = true;
      gap3.push_back(event);
    }
    else GAP3 = false;	

    if (TPC3GAP4 == true && TPC7 == true) {
      GAP4 = true;
      gap4.push_back(event);
    }
    else GAP4 = false;

    if (TPC3 == true && TPC5 == true) {
      GAP5 = true;
      gap5.push_back(event);
    }
    else GAP5 = false;
    
    //Use gap crossers to determine multiple gap crossing events
    if (GAP1 == true && GAP2 == true) cross12.push_back(event);
    if (GAP3 == true && GAP4 == true) cross34.push_back(event);
    if (GAP1 == true && GAP4 == true) cross14.push_back(event);
    if (GAP2 == true && GAP3 == true) cross23.push_back(event);
    if (GAP1 == true && (GAP2 == true || GAP4 == true)) cross1and2or4.push_back(event);
    if (GAP2 == true && (GAP1 == true || GAP3 == true)) cross2and1or3.push_back(event);
    if ((GAP1 == true || GAP3 == true) && (GAP2 == true || GAP4 == true)) cross1or3and2or4.push_back(event);

    //Reset bools for next event

    TPC1     = false; 
    TPC7     = false;
    TPC5GAP1 = false;
    GAP1     = false;    
    TPC5GAP2 = false;
    GAP2     = false;    
    TPC3GAP3 = false;
    GAP3     = false; 
    TPC3GAP4 = false;
    GAP4     = false;
    TPC3     = false;
    TPC5     = false;
    GAP5     = false;
  
    fTree->Fill();
  }
  
  void GapFilter::beginJob()
  {    
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
  
    fTree = tfs->make<TTree>("hitdumper","analysis tree");
    fTree->Branch("run",&run,"run/I");
    fTree->Branch("subrun",&subrun,"subrun/I");
    fTree->Branch("event",&event,"event/I");
    fTree->Branch("evttime",&evttime,"evttime/D");
    fTree->Branch("efield",efield,"efield[3]/D");
    fTree->Branch("t0",&t0,"t0/I");

    fTree->Branch("nhits",&nhits,"nhits/I");
    fTree->Branch("hit_cryostat",hit_cryostat,"hit_cryostat[nhits]/I");
    fTree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits]/I");
    fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I");
    fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I");
    fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I");
    fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/D");
    fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/D");
    fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/D");

    fTree->Branch("gap1", &gap1);
    fTree->Branch("gap2", &gap2);
    fTree->Branch("gap3", &gap3);
    fTree->Branch("gap4", &gap4);
    fTree->Branch("gap5", &gap5);

    fTree->Branch("cross12", &cross12);
    fTree->Branch("cross34", &cross34);

    fTree->Branch("cross14", &cross14);
    fTree->Branch("cross23", &cross23);

    fTree->Branch("cross1and2or4", &cross1and2or4);
    fTree->Branch("cross2and1or3", &cross2and1or3);
    fTree->Branch("cross1or3and2or4", &cross1or3and2or4);
  }

  void GapFilter::ResetVars(){

    run = -99999;
    subrun = -99999;
    event = -99999;
    evttime = -99999;
    for (int i = 0; i<3; ++i){
      efield[i] = -99999;
    }
    t0 = -99999;
    nhits = -99999;
    for (int i = 0; i<kMaxHits; ++i){
      hit_cryostat[i] = -99999;
      hit_tpc[i] = -99999;
      hit_plane[i] = -99999;
      hit_wire[i] = -99999;
      hit_channel[i] = -99999;
      hit_peakT[i] = -99999;
      hit_charge[i] = -99999;
      hit_ph[i] = -99999;
    }
  }
  DEFINE_ART_MODULE(GapFilter)
}  // namespace GapFilter

#endif // GapFilter_Module
