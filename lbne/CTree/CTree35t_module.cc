// Read data from MC raw files and convert it into ROOT tree
// Chao Zhang (chao@bnl.gov) 2/24/2014

#ifndef CTree35t_Module
#define CTree35t_Module
// test

// LArSoft includes
#include "Utilities/DetectorProperties.h"
#include "Utilities/GeometryUtilities.h"
// #include "Utilities/LArProperties.h"

#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
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
#include "art/Persistency/Common/PtrVector.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"

// ROOT includes.
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

// C++ Includes
#include <map>
#include <vector>
// #include <algorithm>
#include <fstream>
#include <iostream>
// #include <string>
// #include <sstream>
// #include <cmath>

// #ifdef __MAKECINT__
// #pragma link C++ class vector<vector<int> >+;
// #pragma link C++ class vector<vector<float> >+;
// #endif

#define MAX_TPC 8
#define MAX_PLANE 3
// #define MAX_CHANNEL 1992
#define MAX_CHANNEL 2048
#define MAX_TRACKS 2000
#define MAX_HITS 20000

using namespace std;

namespace LBNE{

class CTree35t : public art::EDAnalyzer {
public:

    explicit CTree35t(fhicl::ParameterSet const& pset);
    virtual ~CTree35t();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);
    void initOutput();
    void saveChannelWireMap();
    void saveWireGeometry(int plane, int tpc);
    void printGeometry();
    void processRaw(const art::Event& evt);
    void processCalib(const art::Event& evt);
    void processMC(const art::Event& evt);
    void processHits(const art::Event& evt);
    void printEvent();
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;
    std::string fCalibLabel;
    std::string fHitsModuleLabel;
    std::string fOutFileName;
    bool fSaveChannelWireMap;

    art::ServiceHandle<geo::Geometry> fGeom;
    // art::ServiceHandle<util::LArProperties> larp;

    // art::ServiceHandle<art::TFileService> fTfs;
    TFile *fOutFile;
    TTree *fGeoTree;
    TTree *fEventTree;

    // Geometry Tree Leafs 
    int fNcryostats;
    int fNTPC;
    float fTPC_x[MAX_TPC];  // TPC length in x
    float fTPC_y[MAX_TPC];  // TPC length in y
    float fTPC_z[MAX_TPC];  // TPC length in z
    int fNplanes;
    int fPlane_type[MAX_PLANE];  // plane type: 0 == induction, 1 == collection
    int fPlane_view[MAX_PLANE];  // wire orientation: 0 == U, 1 == V, 2 == X
    double fPlane_wirepitch[MAX_PLANE];  // wire pitch of each plane
    double fPlane_wireangle[MAX_PLANE];  // wire angle (to vertical) of each plane
    int fPlane_wires[MAX_PLANE];  // number of wires in each plane
    int fNchannels;


    // Event Tree Leafs
    int fEvent;
    int fRun;
    int fSubRun;

    int fRaw_Nhit;
    int fRaw_channelId[MAX_CHANNEL];  // hit channel id; size == raw_Nhit
    int fRaw_charge[MAX_CHANNEL];  // hit channel charge (simple alg); size == raw_Nhit
    int fRaw_time[MAX_CHANNEL];  // hit channel time (simple alg); size == raw_Nhit
    std::vector<std::vector<int> > fRaw_wfADC;
    std::vector<std::vector<int> > fRaw_wfTDC;

    int fCalib_Nhit;
    int fCalib_channelId[MAX_CHANNEL];  // hit channel id; size == fCalib_Nhit
    // FIXEME:: cannot save e.g std::vector<std::vector<float> > in ttree
    std::vector<std::vector<int> > fCalib_wfADC;  
    std::vector<std::vector<int> > fCalib_wfTDC;


    int fMC_Ntrack;  // number of tracks in MC
    int fMC_id[MAX_TRACKS];  // track id; size == mc_Ntrack
    int fMC_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
    int fMC_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
    float fMC_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
    float fMC_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
    float fMC_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
    float fMC_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
    std::vector<std::vector<int> > fMC_daughters;  // daughters id of this track; vector


    int    no_hits;                  //number of hits
    int    hit_channel[MAX_HITS];    //channel ID
    float  hit_peakT[MAX_HITS];      //peak time
    float  hit_charge[MAX_HITS];     //charge (area)


	 // unsigned int UVPlane[4]={3,2,1,0};
	 // unsigned int ZPlane[8]={7,6,5,4,3,2,1,0};

 }; // class CTree35t


//-----------------------------------------------------------------------
CTree35t::CTree35t(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
    initOutput();
}


//-----------------------------------------------------------------------
CTree35t::~CTree35t()
{
}


//-----------------------------------------------------------------------
void CTree35t::reconfigure(fhicl::ParameterSet const& p){
    fRawDigitLabel = p.get< std::string >("RawDigitLabel");
    fCalibLabel = p.get< std::string >("CalibLabel");
    fHitsModuleLabel = p.get< std::string >("HitsModuleLabel");
    fOutFileName = p.get< std::string >("outFile");
    fSaveChannelWireMap = p.get< bool >("saveChannelWireMap");
}


//-----------------------------------------------------------------------
void CTree35t::initOutput()
{
    TDirectory* tmpDir = gDirectory;

    fOutFile = new TFile(fOutFileName.c_str(), "recreate");

    // init Detector Geometry TTree
    TDirectory* subDir = fOutFile->mkdir("Detector");
    subDir->cd();
    fGeoTree = new TTree("Geometry", "Detector Geometry");
    fGeoTree->Branch("Ncryostats", &fNcryostats, "Ncryostats/I");  // number of cryostats, == 1
    
    fGeoTree->Branch("NTPC"      , &fNTPC);  // number of (virtual) Time Projection Chambers, == 8, NAPA = NTPC/2
    fGeoTree->Branch("TPC_x"     , &fTPC_x, "TPC_x[NTPC]/F"); // TPC length in x
    fGeoTree->Branch("TPC_y"     , &fTPC_y, "TPC_y[NTPC]/F"); // TPC length in y
    fGeoTree->Branch("TPC_z"     , &fTPC_z, "TPC_z[NTPC]/F"); // TPC length in z
    
    fGeoTree->Branch("Nplanes"     , &fNplanes);  // number of wire planes in each TPC, == 3
    fGeoTree->Branch("plane_type"  , &fPlane_type, "plane_type[Nplanes]/I"); // plane type: 0 == induction, 1 == collection
    fGeoTree->Branch("plane_view"  , &fPlane_view, "plane_view[Nplanes]/I"); // wire orientation: 0 == U, 1 == V, 2 == X
    fGeoTree->Branch("plane_wirepitch"  , &fPlane_wirepitch, "plane_wirepitch[Nplanes]/D"); // wire pitch of each plane
    fGeoTree->Branch("plane_wireangle"  , &fPlane_wireangle, "plane_wireangle[Nplanes]/D"); // wire pitch of each plane
    fGeoTree->Branch("plane_wires" , &fPlane_wires, "Plane_wires[Nplanes]/I"); // number of wires in each plane

    fGeoTree->Branch("Nchannels" , &fNchannels);  // number of total channels


    // init Event TTree
    TDirectory* subDir2 = fOutFile->mkdir("Event");
    subDir2->cd();
    fEventTree = new TTree("Sim", "Event Tree from Simulation");
    fEventTree->Branch("eventNo", &fEvent);
    fEventTree->Branch("runNo", &fRun);
    fEventTree->Branch("subRunNo", &fSubRun);

    fEventTree->Branch("raw_Nhit", &fRaw_Nhit);  // number of hit channels above threshold
    fEventTree->Branch("raw_channelId" , &fRaw_channelId, "raw_channelId[raw_Nhit]/I"); // hit channel id; size == raw_Nhit
    fEventTree->Branch("raw_charge" , &fRaw_charge, "raw_charge[raw_Nhit]/I"); // hit channel id; size == raw_Nhit
    fEventTree->Branch("raw_time" , &fRaw_time, "raw_time[raw_Nhit]/I"); // hit channel id; size == raw_Nhit
    fEventTree->Branch("raw_wfADC", &fRaw_wfADC);  // raw waveform adc of each channel
    fEventTree->Branch("raw_wfTDC", &fRaw_wfTDC);  // raw waveform tdc of each channel

    fEventTree->Branch("calib_Nhit", &fCalib_Nhit);  // number of hit channels above threshold
    fEventTree->Branch("calib_channelId" , &fCalib_channelId, "calib_channelId[calib_Nhit]/I"); // hit channel id; size == calib_Nhit
    fEventTree->Branch("calib_wfADC", &fCalib_wfADC);  // calib waveform adc of each channel
    fEventTree->Branch("calib_wfTDC", &fCalib_wfTDC);  // calib waveform tdc of each channel


    fEventTree->Branch("mc_Ntrack", &fMC_Ntrack);  // number of tracks in MC
    fEventTree->Branch("mc_id", &fMC_id, "mc_id[mc_Ntrack]/I");  // track id; size == mc_Ntrack
    fEventTree->Branch("mc_pdg", &fMC_pdg, "mc_id[mc_Ntrack]/I");  // track particle pdg; size == mc_Ntrack
    fEventTree->Branch("mc_mother", &fMC_mother, "mc_mother[mc_Ntrack]/I");  // mother id of this track; size == mc_Ntrack
    fEventTree->Branch("mc_daughters", &fMC_daughters);  // daughters id of this track; vector
    fEventTree->Branch("mc_startXYZT", &fMC_startXYZT, "mc_startXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
    fEventTree->Branch("mc_endXYZT", &fMC_endXYZT, "mc_endXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
    fEventTree->Branch("mc_startMomentum", &fMC_startMomentum, "mc_startMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
    fEventTree->Branch("mc_endMomentum", &fMC_endMomentum, "mc_endMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack


    fEventTree->Branch("no_hits", &no_hits);  //number of hits
    fEventTree->Branch("hit_channel", &hit_channel, "hit_channel[no_hits]/I");  // channel ID
    fEventTree->Branch("hit_peakT", &hit_peakT, "hit_peakT[no_hits]/F");  // peak time
    fEventTree->Branch("hit_charge", &hit_charge, "hit_charge[no_hits]/F");  // charge (area)

    gDirectory = tmpDir;

}

//-----------------------------------------------------------------------
void CTree35t::beginJob()
{
    fNcryostats = fGeom->Ncryostats();  // 1 for 35t

    // 8 TPC for 35t
    // TPC (0,1), (6,7): large APA
    // TPC (2,3), (4,5): small APA
    // 0,2,4,6: short drift; 1,3,5,7: long drift
    fNTPC = fGeom->NTPC();  
    for (int i=0; i<fNTPC; i++) {
        fTPC_x[i] = fGeom->DetHalfWidth(i)*2;
        fTPC_y[i] = fGeom->DetHalfHeight(i)*2;
        fTPC_z[i] = fGeom->DetLength(i);
    }

    // 3 planes for 35t. geo::kCollection plane: 1; induction plane: 0
    // plane 0( type: 0, Nwires: 357)
    // plane 1( type: 0, Nwires: 343)
    // plane 2( type: 1, Nwires: 111)
    fNplanes = fGeom->Nplanes();
    for (int i=0; i<fNplanes; i++) {
        fPlane_type[i] = fGeom->Plane(i).SignalType();
        fPlane_view[i] = fGeom->Plane(i).View();
        // fPlane_wirepitch[i] = fGeom->WirePitch(fPlane_view[i]);  // this doesn't seem to return the correct value!
        fPlane_wirepitch[i] = fGeom->WirePitch(0, 1, fPlane_view[i], 1, 0);  // this doesn't seem to return the correct value!
        fPlane_wireangle[i] = fGeom->WireAngleToVertical(fGeom->Plane(i).View());
        fPlane_wires[i] = fGeom->Nwires(i);
    }
    
    fNchannels = fGeom->Nchannels();

    printGeometry();

    // for (int i=0; i<fNplanes; i++) {
    //   int nw = fGeom->Nwires(i);
    //   for (int j=0; j<nw; j++) {
    //     cout << "channel (" << i << ", " << j << "): " << fGeom->PlaneWireToChannel(i, j) << endl;
    //   }
    
    // }

    // Write fGeoTree to Disk (once)
    fGeoTree->Fill();
    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Detector");
    fGeoTree->Write();
    gDirectory = tmpDir;
    

    // Save Channel Map to text file.
    if (fSaveChannelWireMap) {
        saveChannelWireMap();
    }
    saveWireGeometry(1, 1);
    saveWireGeometry(1, 3);
    saveWireGeometry(1, 5);
    saveWireGeometry(1, 7);

}


//-----------------------------------------------------------------------
void CTree35t::saveChannelWireMap()
{
    ofstream mapfile;
    mapfile.open("ChannelWireMap.txt");
    mapfile << "# total channels: " << fNchannels << "\n\n";
    for (int i=0; i<fNchannels; i++) {
        std::vector<geo::WireID> wireids = fGeom->ChannelToWire(i);
        mapfile << "Channel " << i << "\n";
        mapfile <<  "Nwires " << wireids.size() << "\n";

        mapfile << "TPC ";
        for (auto const& wid : wireids) {
            mapfile << wid.TPC << " ";
        }
        mapfile << "\n";

        mapfile << "Plane ";
        for (auto const& wid : wireids) {
            mapfile << wid.Plane << " ";
        }
        mapfile << "\n";

        // redundant: Plane defines View
        // mapfile << "View ";
        // for (auto const& wid : wireids) {
        //   mapfile << fGeom->Plane(wid.Plane, wid.TPC).View() << " ";
        // }
        // mapfile << "\n";

        mapfile << "Wire ";
        for (auto const& wid : wireids) {
            mapfile << wid.Wire << " ";
        }
        mapfile << "\n" << endl;
    }

    mapfile.close();
}


//-----------------------------------------------------------------------
void CTree35t::saveWireGeometry(int plane, int tpc)
{
    int cstat = 0;
    int Nwires = fGeom->Nwires(plane, tpc);
    double xyzStart[3];
    double xyzEnd[3];
    for (int wire=0; wire<Nwires; wire++) {
        fGeom->WireEndPoints(cstat, tpc, plane, wire, xyzStart, xyzEnd);
        cout << plane << "\t" << wire << "\t";
        for (int i=0; i<3; i++) {
            cout << xyzStart[i] << "\t";
        }
        for (int i=0; i<3; i++) {
            cout << xyzEnd[i] << "\t";
        }
        cout << endl;
    }

    // cout << " Temperature: " << larp->Temperature() << endl;
    // cout << " E field: " << larp->Efield() << endl;
    // cout << " Drift Velocity: " << larp->DriftVelocity(larp->Efield(), larp->Temperature()) << endl;

}

//-----------------------------------------------------------------------
void CTree35t::endJob()
{
    // Write fGeoTree to fEventTree
    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Event");
    fEventTree->Write();
    gDirectory = tmpDir;
    fOutFile->Close();
}


//-----------------------------------------------------------------------
void CTree35t::printGeometry()
{

    cout << "fNcryostats: " << fNcryostats << endl;
    cout << "fNTPC: " << fNTPC << endl;
    for (int i=0; i<fNTPC; i++) {
        cout << "\tTPC " << i << ": " << fTPC_x[i] << ", " << fTPC_y[i] << ", " << fTPC_z[i] << endl;
    }
    cout << "fNplanes: " << fNplanes << endl;
    for (int i=0; i<fNplanes; i++) {
        cout 
            << "\tplane " << i 
            << "( type: " << fPlane_type[i]
            << ", view: " << fPlane_view[i]
            << ", wirepitch: " << fPlane_wirepitch[i]
            << ", wire angle: " << fPlane_wireangle[i]
            << ", wires: " << fPlane_wires[i]
            << ")" << endl;
    }
    cout << "fNchannels: " << fNchannels << endl;
    cout << "fNOpDet: " << fGeom->NOpDet() << endl;
    cout << "fAuxDetectors: " << fGeom->NAuxDets() << endl;
    cout << endl;
}


//-----------------------------------------------------------------------
void CTree35t::beginRun(const art::Run& /*run*/)
{
  mf::LogInfo("CTree35t")  << "begin run";
}


//-----------------------------------------------------------------------
void CTree35t::analyze( const art::Event& event )
{
    reset();

    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    processMC(event);
    processRaw(event);
    processCalib(event);
    processHits(event);
    printEvent();
    fEventTree->Fill();
}


//-----------------------------------------------------------------------
void CTree35t::reset()
{
    fRaw_wfADC.clear();
    fRaw_wfTDC.clear();
    fCalib_wfADC.clear();
    fCalib_wfTDC.clear();
    for (int i=0; i<MAX_CHANNEL; i++) {
        fRaw_channelId[i] = 0;
        fRaw_charge[i] = 0;
        fRaw_time[i] = 0;
        fCalib_channelId[i] = 0;
    }

    for (int i=0; i<MAX_TRACKS; i++) {
        fMC_id[i] = 0;
        fMC_pdg[i] = 0;
        fMC_mother[i] = 0;
        for (int j=0; j<4; j++) {
            fMC_startXYZT[i][j]      = 0;
            fMC_endXYZT[i][j]        = 0;
            fMC_startMomentum[i][j] = 0;
            fMC_endMomentum[i][j]   = 0;
        }
    }
    fMC_daughters.clear();

    for (int i=0; i<MAX_HITS; i++) {
        hit_channel[i] = 0;
        hit_peakT[i] = 0;
        hit_charge[i] = 0;
    }

}


//-----------------------------------------------------------------------
void CTree35t::processRaw( const art::Event& event )
{
    //    unsigned int tpcNo, cryoNo;
    // get the objects holding all of the raw data information
    art::Handle< std::vector<raw::RawDigit> > rawdigit;
    event.getByLabel(fRawDigitLabel, rawdigit);

    // put it in a more easily usable form
    std::vector< art::Ptr<raw::RawDigit> >  rawhits;
    art::fill_ptr_vector(rawhits, rawdigit);

    // rawhits size should == Nchannels == 1992; (no hit channel has a flat 0-waveform)
    // cout << "\n Raw Hits size: " << rawhits.size() << endl;

    //loop through all RawDigits (over entire channels)
    fRaw_Nhit = 0;
    for (auto const& hit : rawhits) {      
        int chanId = hit->Channel();
        int nSamples = hit->Samples();
        int pedstal = hit->GetPedestal(); // should be 0 as of 2/25

        std::vector<short> uncompressed(nSamples);
        raw::Uncompress(hit->fADC, uncompressed, hit->Compression());
        // uncompressed size is 3200 samples per waveform
        short thresh = pedstal + 1; // threshold set to 1 adc;
        bool isHit = false;
        for (auto const& adc : uncompressed) {
            if (adc > thresh) {
                isHit = true;
                break;
            }
        }
        if (!isHit) continue; // skip empty channels

        int id = fRaw_Nhit;
        fRaw_channelId[id] = chanId;

        vector<int> wfADC;
        vector<int> wfTDC;
        int nSavedSamples = 0;
        bool hasTime = false;
        for (int i=0; i<nSamples; i++) {
            short adc = uncompressed[i];
            if (adc != pedstal) {
                // cout << i << "," << adc << " | ";
                nSavedSamples++;
                wfADC.push_back(int(adc));
                wfTDC.push_back(i);
                fRaw_charge[id] += adc;
                if (!hasTime) {
                    fRaw_time[id] = i;
                    hasTime = true;
                }
            }
        }

        fRaw_wfADC.push_back(wfADC);
        fRaw_wfTDC.push_back(wfTDC);
        fRaw_Nhit++;
      // std::vector<geo::WireID> wireids = fGeom->ChannelToWire(chanNo);
      // cout 
      //   << "\n channelID: " << fRaw_channelId[id] 
      //   << "\n charge: " << fRaw_charge[id] 
      //   << "\n time: " << fRaw_time[id] 
      //   << "\n nSamples: " << nSamples
      //   << "\n pedestal: " << pedstal
      //   << "\n nSavedSamples: " << nSavedSamples
      //   << endl;

    }

}

void CTree35t::processCalib( const art::Event& event )
{

    art::Handle< std::vector<recob::Wire> > wires_handle;
    if (! event.getByLabel(fCalibLabel, wires_handle)) return;

    // put it in a more easily usable form
    std::vector< art::Ptr<recob::Wire> >  wires;
    art::fill_ptr_vector(wires, wires_handle);

    // wires size should == Nchannels == 1992; (no hit channel has a flat 0-waveform)
    // cout << "\n wires size: " << wires.size() << endl;

    fCalib_Nhit = 0;
    for (auto const& wire : wires) {      
        std::vector<float> calibwf = wire->Signal(); 
        int chanId = wire->Channel();
        int nSamples = calibwf.size();
        int pedstal = 0;

        short thresh = pedstal + 1; // threshold set to 1 adc;
        bool isHit = false;
        for (auto const& adc : calibwf) {
            if (adc > thresh) {
                isHit = true;
                break;
            }
        }
        if (!isHit) continue; // skip empty channels

        int id = fCalib_Nhit;
        fCalib_channelId[id] = chanId;

        vector<int> wfADC;
        vector<int> wfTDC;
        int nSavedSamples = 0;
        for (int i=0; i<nSamples; i++) {
            int adc = int(calibwf[i]);
            if (adc != pedstal) {
                // cout << i << "," << adc << " | ";
                nSavedSamples++;
                wfADC.push_back(adc);
                wfTDC.push_back(i);
            }
        }
        // cout << endl;

        fCalib_wfADC.push_back(wfADC);
        fCalib_wfTDC.push_back(wfTDC);
        fCalib_Nhit++;
    }

}

void CTree35t::processMC( const art::Event& event )
{
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    event.getByLabel("largeant", particleHandle);

    // put it in a more easily usable form
    std::vector< art::Ptr<simb::MCParticle> > particles;
    art::fill_ptr_vector(particles, particleHandle);

    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel("largeant", simChannelHandle);    

    fMC_Ntrack = particles.size();

    int i=0; // track index
    for (auto const& particle: particles ) {
        fMC_id[i] = particle->TrackId();
        fMC_pdg[i] = particle->PdgCode();
        fMC_mother[i] = particle->Mother();
        int Ndaughters = particle->NumberDaughters();
        vector<int> daughters;
        for (int i=0; i<Ndaughters; i++) {
            daughters.push_back(particle->Daughter(i));
        }
        fMC_daughters.push_back(daughters);
        size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
        int last = numberTrajectoryPoints - 1;
        const TLorentzVector& positionStart = particle->Position(0);
        const TLorentzVector& positionEnd   = particle->Position(last);
        const TLorentzVector& momentumStart = particle->Momentum(0);
        const TLorentzVector& momentumEnd   = particle->Momentum(last);
        positionStart.GetXYZT(fMC_startXYZT[i]);
        positionEnd.GetXYZT(fMC_endXYZT[i]);
        momentumStart.GetXYZT(fMC_startMomentum[i]);
        momentumEnd.GetXYZT(fMC_endMomentum[i]);
        i++;
    } // particle loop done 
}


void CTree35t::processHits( const art::Event& event )
{
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    if (! event.getByLabel(fHitsModuleLabel, hitListHandle)) return;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    art::fill_ptr_vector(hitlist, hitListHandle);

    no_hits = hitlist.size();
    if (no_hits>MAX_HITS) {
        mf::LogError("CTree35t") << "Event has " << no_hits
        << " hits, MAX ALLOWED: " << MAX_HITS;
    }
    for (int i=0; i<no_hits; i++) {
        art::Ptr<recob::Hit> hit = hitlist[i];
        hit_channel[i] = hit->Channel();
        hit_charge[i] = hit->Charge();
        hit_peakT[i] = hit->PeakTime();        
    }
}


void CTree35t::printEvent()
{
    cout << " Event: " << fEvent << endl;
    cout << "   Run: " << fRun << endl;
    cout << "SubRun: " << fSubRun << endl;

    cout << " Hit Channels: " << fRaw_Nhit << endl;
    cout << "      Ntracks:" << fMC_Ntrack;
    for (int i=0; i<fMC_Ntrack; i++) {
        cout << "\n              id: " << fMC_id[i];
        cout << "\n             pdg: " << fMC_pdg[i];
        cout << "\n          mother: " << fMC_mother[i];
        cout << "\n      Ndaughters: " << fMC_daughters.at(i).size();
        cout << "\n      start XYZT: (" << fMC_startXYZT[i][0] << ", " << fMC_startXYZT[i][1] << ", " << fMC_startXYZT[i][2] << ", " << fMC_startXYZT[i][3] << ")";
        cout << "\n        end XYZT: (" << fMC_endXYZT[i][0] << ", " << fMC_endXYZT[i][1] << ", " << fMC_endXYZT[i][2] << ", " << fMC_endXYZT[i][3] << ")";
        cout << "\n  start momentum: (" << fMC_startMomentum[i][0] << ", " << fMC_startMomentum[i][1] << ", " << fMC_startMomentum[i][2] << ", " << fMC_startMomentum[i][3] << ")";
        cout << "\n    end momentum: (" << fMC_endMomentum[i][0] << ", " << fMC_endMomentum[i][1] << ", " << fMC_endMomentum[i][2] << ", " << fMC_endMomentum[i][3] << ")";

        cout << endl;
    }

    cout << "Number of Hits Found: " << no_hits << endl;
    // for (int i=0; i<no_hits; i++) {
    //     cout << hit_channel[i] << ", ";
    //     cout << hit_charge[i] << ", ";
    //     cout << hit_peakT[i] << ", ";
    //     cout << endl;
    // }

}

DEFINE_ART_MODULE(CTree35t)

} // namespace RawDataReader

#endif // CTree35t_Module
