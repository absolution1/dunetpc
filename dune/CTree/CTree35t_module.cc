// Read data from MC raw files and convert it into ROOT tree
// Chao Zhang (chao@bnl.gov) 2/24/2014

#ifndef CTree35t_Module
#define CTree35t_Module
// test

// LArSoft includes
//#include "lardata/Utilities/DetectorProperties.h"
#include "lardata/Utilities/GeometryUtilities.h"
// #include "Utilities/LArProperties.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Simulation/SimListUtils.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larana/OpticalDetector/OpDetResponseInterface.h"
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
// lbne-artdaq includes
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh" 
#include "lbne-raw-data/Overlays/SSPFragment.hh" 
#include "lbne-raw-data/Overlays/anlTypes.hh"
#include "artdaq-core/Data/Fragment.hh" 
#include "../daqinput35t/tpcFragmentToRawDigits.h"
#include "../daqinput35t/PennToOffline.h"
#include "../daqinput35t/SSPReformatterAlgs.h"
#include "../daqinput35t/utilities/UnpackFragment.h"
// ROOT includes.
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TGeoTube.h"
#include "TGeoNode.h"
#include "TH1D.h"
// C++ Includes
#include <map>
#include <vector>
// #include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
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
#define MAX_OPDET 8
#define MAX_CLUSTER 10
#define MAX_OPWAVEFORMS 1000

using namespace std;

namespace DUNE{

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
    void processRecoTracks(const art::Event& event);
    void processOpDet(const art::Event& event);
    void processTiming(const art::Event& event);
    void printEvent();
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;
    std::string fCalibLabel;
    std::string fHitsModuleLabel;
    std::string fOutFileName;
    std::string fTrackModuleLabel;
    std::string fClusterModuleLabel; 
    std::string fOpDetInputModule;
    std::string fInstanceName;
    std::string fOpHitModule;
    bool fSaveChannelWireMap;
    bool fSaveChannelWireGeo;
  
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
    int fNOpDets;
    float fOpDetPositions_Y[MAX_OPDET];
    float fOpDetPositions_Z[MAX_OPDET];
    float fOpDetHalfWidths[MAX_OPDET];
    float fOpDetHalfHeights[MAX_OPDET];

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
    TObjArray *fMC_trackMomentum;
    TObjArray *fMC_trackPosition;

    int    no_hits;                  //number of hits
    int    hit_channel[MAX_HITS];    //channel ID
    float  hit_peakT[MAX_HITS];      //peak time
    float  hit_charge[MAX_HITS];     //charge (area)
    float  hit_wireID[MAX_HITS];
    float  hit_plane[MAX_HITS]; 
    float  hit_tpc[MAX_HITS]; 

    int    nthits;
    int    thit_channel[MAX_HITS];
    float  thit_peakT[MAX_HITS];
    float  thit_charge[MAX_HITS];
    int    thit_wireID[MAX_HITS];
    int    thit_plane[MAX_HITS];
    int    thit_tpc[MAX_HITS];

    int    nclhits;
    int    chit_cryostat[MAX_HITS];
    int    chit_tpc[MAX_HITS];
    int    chit_plane[MAX_HITS];
    int    chit_charge[MAX_HITS];
    float  chit_peakT[MAX_HITS];
    int    chit_wire[MAX_HITS];
    int    chit_channel[MAX_HITS];
    int    chit_cluster[MAX_HITS];

    int reco_nTrack;    // number of reco tracks
    TObjArray *fReco_trackPosition;
  
  // unsigned int UVPlane[4]={3,2,1,0};
  // unsigned int ZPlane[8]={7,6,5,4,3,2,1,0};

  //Photon detector related
    TTree * fThePhotonTreeAll;
    TTree * fThePhotonTreeDetected;
    TTree * fTheOpDetTree;
    TTree * fTheEventTree;
    std::string fInputModule;
    bool fMakeDetectedPhotonsTree;
    bool fMakeAllPhotonsTree;
    bool fMakeOpDetsTree;
    bool fMakeOpDetEventsTree;
    bool fUncompressWithPed;
    bool fProcessMCtruth;
    bool fProcessCalib;
    bool fProcessHits;
    bool fProcessReco;
    bool fProcessOpDet;
    //float fQE; // unused
    //float fWavelengthCutLow; // unused
    //float fWavelengthCutHigh; // unused
    Float_t fWavelength;
    Float_t fTime;
    Int_t fCount;
    Int_t fCountOpDetAll[MAX_OPDET];
    Int_t fCountOpDetDetected[MAX_OPDET];
    Int_t fCountEventAll;
    Int_t fCountEventDetected;
    Int_t fOpChannel;
  // new photon detector display
  TObjArray *averageWaveforms;
  std::map<int, int> waveformCount; 
  std::vector<std::vector<int> > OpChannelToOpDet;
  std::vector<std::vector<int> > timestamp;
  //std::map< int, int   > OpHitCount; 
  //std::map< int, double   > FirstHitTimePerChannel;
  /*
  TH1D*       fPedestalMeanPerChannel;
  TH1D*       fPedestalSigmaPerChannel;
  TH1D*       fIntegratedSignalMeanPerChannel;
  TH1D*       fFractionSamplesNearMaximumPerChannel;
  TH1D*       fNumberOfWaveformsProcessedPerChannel;
  TH1D*       fFirstOpHitTimeMean;
  //TH1D*       fFirstOpHitTimeSigma;
  TH1D*       fSecondOpHitTimeMean;
  //TH1D*       fSecondOpHitTimeSigma;
  TH1D*       fFirstSecondDiffOpHitTimeMean;
  //TH1D*       fFirstSecondDiffOpHitTimeSigma; 
  TH1D*       fNumberOfOpHitsPerChannelPerEvent;
  */

  // timing
    std::string fRCERawDataLabel, fRCEFragType;
    std::string fSSPFragType, fSSPRawDataLabel;
    double fSampleFreq;
    double RCETimeBegin;
    double RCETimeEnd;
    double SSPTimeBegin;
    double SSPTimeEnd;

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
    fClusterModuleLabel = p.get< std::string >("ClusterModuleLabel");
    fTrackModuleLabel = p.get< std::string >("TrackModuleLabel");
    fOutFileName = p.get< std::string >("outFile");
    fSaveChannelWireMap = p.get< bool >("saveChannelWireMap");
    fSaveChannelWireGeo = p.get< bool >("saveChannelWireGeo");
    fInputModule = p.get<std::string>("InputModule");
    fMakeAllPhotonsTree = p.get<bool>("MakeAllPhotonsTree");
    fMakeDetectedPhotonsTree = p.get<bool>("MakeDetectedPhotonsTree");
    fMakeOpDetsTree = p.get<bool>("MakeOpDetsTree");
    fMakeOpDetEventsTree = p.get<bool>("MakeOpDetEventsTree");
    fUncompressWithPed  = p.get< bool         >("UncompressWithPed", true);
    fProcessMCtruth = p.get< bool         >("ProcessMCtruth", true);
    fProcessCalib = p.get< bool         >("ProcessCalib", true);
    fProcessHits = p.get< bool         >("ProcessHits", true);
    fProcessReco = p.get< bool         >("ProcessReco", true);
    fProcessOpDet = p.get< bool         >("ProcessOpDet", true);
    fOpDetInputModule = p.get<std::string >("OpDetInputModule");
    fInstanceName = p.get<std::string >("InstanceName");
    fOpHitModule = p.get<std::string >("OpHitModule");
    fRCEFragType = p.get<std::string >("RCEFragType");
    fRCERawDataLabel = p.get<std::string >("RCERawDataLabel");
    fSSPFragType = p.get<std::string >("SSPFragType");
    fSSPRawDataLabel = p.get<std::string >("SSPRawDataLabel");
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

    fGeoTree->Branch("NOpDets", &fNOpDets, "NOpDets/I");
    fGeoTree->Branch("OpDetPositions_Y", &fOpDetPositions_Y, "OpDetPositions_Y[NOpDets]/F");
    fGeoTree->Branch("OpDetPositions_Z", &fOpDetPositions_Z, "OpDetPositions_Z[NOpDets]/F");
    fGeoTree->Branch("OpDetHalfWidths", &fOpDetHalfWidths, "OpDetHalfWidths[NOpDets]/F");
    fGeoTree->Branch("OpDetHalfHeights", &fOpDetHalfHeights, "OpDetHalfHeights[NOpDets]/F");

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

    fMC_trackPosition = new TObjArray();
    fMC_trackMomentum= new TObjArray();
    fMC_trackPosition->SetOwner(kTRUE);
    fMC_trackMomentum->SetOwner(kTRUE);

    fEventTree->Branch("mc_trackPosition", &fMC_trackPosition);
    fEventTree->Branch("mc_trackMomentum", &fMC_trackMomentum);

    fEventTree->Branch("no_hits", &no_hits);  //number of hits
    fEventTree->Branch("hit_channel", &hit_channel, "hit_channel[no_hits]/I");  // channel ID
    fEventTree->Branch("hit_peakT", &hit_peakT, "hit_peakT[no_hits]/F");  // peak time
    fEventTree->Branch("hit_charge", &hit_charge, "hit_charge[no_hits]/F");  // charge (area)
    fEventTree->Branch("hit_wireID", &hit_wireID, "hit_wireID[no_hits]/F");
    fEventTree->Branch("hit_plane", &hit_plane, "hit_plane[no_hits]/F");
    fEventTree->Branch("hit_tpc", &hit_tpc, "hit_tpc[no_hits]/F");

    fEventTree->Branch("nthits", &nthits);
    fEventTree->Branch("thit_channel", &thit_channel, "thit_channel[nthits]/I");
    fEventTree->Branch("thit_peakT", &thit_peakT, "thit_peakT[nthits]/F");
    fEventTree->Branch("thit_charge", &thit_charge, "thit_charge[nthits]/F");
    fEventTree->Branch("thit_wireID", &thit_wireID, "thit_wireID[nthits]/I");
    fEventTree->Branch("thit_plane", &thit_plane, "thit_plane[nthits]/I");
    fEventTree->Branch("thit_tpc", &thit_tpc, "thit_tpc[nthits]/I");

    fEventTree->Branch("nclhits", &nclhits, "nclhits/I");
    fEventTree->Branch("chit_cryostat", &chit_cryostat, "chit_cryostat[nclhits]/I");
    fEventTree->Branch("chit_tpc", &chit_tpc, "chit_tpc[nclhits]/I");
    fEventTree->Branch("chit_plane", &chit_plane, "chit_plane[nclhits]/I");
    fEventTree->Branch("chit_charge", &chit_charge, "chit_charge[nclhits]/I");
    fEventTree->Branch("chit_peakT", &chit_peakT, "chit_peakT[nclhits]/F");
    fEventTree->Branch("chit_wire", &chit_wire, "chit_wire[nclhits]/I");
    fEventTree->Branch("chit_channel", &chit_channel, "chit_channel[nclhits]/I");
    fEventTree->Branch("chit_cluster", &chit_cluster, "chit_cluster[nclhits]/I");

    fEventTree->Branch("reco_nTrack", &reco_nTrack); 
    fReco_trackPosition = new TObjArray();
    fReco_trackPosition->SetOwner(kTRUE);
    fEventTree->Branch("reco_trackPosition", &fReco_trackPosition);

    averageWaveforms = new TObjArray();
    averageWaveforms->SetOwner(kTRUE);
    fEventTree->Branch("averageWaveforms", &averageWaveforms);
    fEventTree->Branch("waveformCount", &waveformCount);
    fEventTree->Branch("OpChannelToOpDet", &OpChannelToOpDet);
    fEventTree->Branch("timestamp", &timestamp);

    fEventTree->Branch("sampleFreq", &fSampleFreq, "sampleFreq/D");
    fEventTree->Branch("RCETimeBegin", &RCETimeBegin, "RCETimeBegin/D");
    fEventTree->Branch("RCETimeEnd", &RCETimeEnd, "RCETimeEnd/D");
    fEventTree->Branch("SSPTimeBegin", &SSPTimeBegin, "SSPTimeBegin[1000]/D");
    fEventTree->Branch("SSPTimeEnd", &SSPTimeEnd, "SSPTimeEnd[1000]/D");
    gDirectory = tmpDir;

    // init Photon TTree
    TDirectory *subDir3 = fOutFile->mkdir("OpDet");
    subDir3->cd();
    if(fMakeAllPhotonsTree){
      fThePhotonTreeAll = new TTree("AllPhotons","AllPhotons");
      fThePhotonTreeAll->Branch("eventNo", &fEvent);
      fThePhotonTreeAll->Branch("runNo", &fRun);
      fThePhotonTreeAll->Branch("subRunNo", &fSubRun);
      fThePhotonTreeAll->Branch("Wavelength", &fWavelength, "Wavelength/F");
      fThePhotonTreeAll->Branch("OpChannel", &fOpChannel, "OpChannel/I");
      fThePhotonTreeAll->Branch("Time", &fTime, "Time/F");
    }
    if(fMakeDetectedPhotonsTree){
      fThePhotonTreeDetected = new TTree("DetectedPhotons", "DetectedPhotons");
      fThePhotonTreeDetected->Branch("eventNo", &fEvent);
      fThePhotonTreeDetected->Branch("runNo", &fRun);
      fThePhotonTreeDetected->Branch("subRunNo", &fSubRun);
      fThePhotonTreeDetected->Branch("Wavelength", &fWavelength, "Wavelength/F");
      fThePhotonTreeDetected->Branch("OpChannel", &fOpChannel, "OpChannel/I");
      fThePhotonTreeDetected->Branch("Time", &fTime, "Time/F");
    }
    if(fMakeOpDetsTree){
      fTheOpDetTree = new TTree("OpDets", "OpDets");
      fTheOpDetTree->Branch("eventNo", &fEvent);
      fTheOpDetTree->Branch("runNo", &fRun);
      fTheOpDetTree->Branch("subRunNo", &fSubRun);
      fTheOpDetTree->Branch("NOpDets", &fNOpDets, "NOpDets/I");
      fTheOpDetTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
      fTheOpDetTree->Branch("CountOpDetAll", &fCountOpDetAll, "CountOpDetAll[NOpDets]/I");
      fTheOpDetTree->Branch("CountOpDetDetected", &fCountOpDetDetected, "CountOpDetDetected[NOpDets]/I");
    }
    if(fMakeOpDetEventsTree){
      fTheEventTree = new TTree("OpDetEvents", "OpDetEvents");
      fTheEventTree->Branch("eventNo", &fEvent);
      fTheEventTree->Branch("runNo", &fRun);
      fTheEventTree->Branch("subRunNo", &fSubRun);
      fTheEventTree->Branch("CountAll", &fCountEventAll, "CountAll/I");
      fTheEventTree->Branch("CountDetected", &fCountEventDetected, "CountDetected/I");
    }


}

//-----------------------------------------------------------------------
void CTree35t::beginJob()
{
  //art::ServiceHandle< util::TimeService> timeService;
    auto const* timeService = lar::providerFrom<detinfo::DetectorClocksService>();
    fSampleFreq = timeService->OpticalClock().Frequency();

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
        fPlane_type[i] = fGeom->SignalType(geo::PlaneID(0, 0, i));
        fPlane_view[i] = fGeom->Plane(i).View();
        // fPlane_wirepitch[i] = fGeom->WirePitch(fPlane_view[i]);  // this doesn't seem to return the correct value!
        fPlane_wirepitch[i] = fGeom->WirePitch(fPlane_view[i], 1, 0);  // this doesn't seem to return the correct value!
        fPlane_wireangle[i] = fGeom->WireAngleToVertical(fGeom->Plane(i).View());
        fPlane_wires[i] = fGeom->Nwires(i);
    }
    
    fNchannels = fGeom->Nchannels();

    // photon detector
    fNOpDets = fGeom->NOpDets();
    OpChannelToOpDet.resize(MAX_OPDET);
    timestamp.resize(MAX_OPDET);
    double xyz[3];
    double tmp;
    for (int i=0; i<fNOpDets; i++) {
      const geo::OpDetGeo& fOpDetNode = fGeom->Cryostat(0).OpDet(i);
      fOpDetNode.GetCenter(xyz,0.);
      fOpDetPositions_Y[i] = (float)xyz[1];
      fOpDetPositions_Z[i] = (float)xyz[2];
      TGeoNode *fOpNode = (TGeoNode*)fOpDetNode.Node();
      TGeoTube *fOpTube = (TGeoTube*)fOpNode->GetVolume()->GetShape();
      tmp = fOpTube->GetDZ();//GetRmax();//fOpDetNode.RMax();
      fOpDetHalfWidths[i] = (float)tmp;
      tmp = fOpTube->GetDY();//fOpDetNode.HalfL();
      fOpDetHalfHeights[i] = (float)tmp;
    }

    printGeometry();

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
    // saveWireGeometry(1, 1);
    // saveWireGeometry(1, 3);
    // saveWireGeometry(1, 5);
    // saveWireGeometry(1, 7);

    /*
    std::ofstream wireGeoFile;
    wireGeoFile.open("WireGeometry.txt");
    for (unsigned int plane=0; plane<(unsigned int)fNplanes; plane++) {
      wireGeoFile << "***************** PLANE " << plane <<" ********************\n";
      for(unsigned int tpc=0; tpc<(unsigned int)fNTPC; tpc++) {
	wireGeoFile << "----------- TPC "<< tpc << " ------------\n";
	wireGeoFile << "Wire#    WireID    ChannelID     Start        End\n";
	int Nwires = fGeom->Nwires(plane, tpc);
	double xyzStart[3];
	double xyzEnd[3];
	for(int wire=0; wire<Nwires; wire++) {
	  fGeom->WireEndPoints(0, tpc, plane, wire, xyzStart, xyzEnd);
	  uint32_t channelid = fGeom->PlaneWireToChannel(plane, wire, tpc, 0);
	  int wireid = wire;
	  wireGeoFile << "                                " << xyzStart[0] << " " << xyzEnd[0] <<"\n";
	  wireGeoFile << wire << "         " << wireid << "         ";
	  wireGeoFile << channelid << "            " << xyzStart[1] << " " << xyzEnd[1] << "\n";
	  wireGeoFile << "                                " << xyzStart[2] << " " << xyzEnd[2] <<"\n";
	}
	wireGeoFile << "------------------------------------------\n\n";
      }
      wireGeoFile << "\n***********************************************\n\n";
    }
    wireGeoFile << "\n" << endl;
    wireGeoFile.close();
    */
}


//-----------------------------------------------------------------------
void CTree35t::saveChannelWireMap()
{
    std::ofstream mapfile;
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
    fOutFile->cd("/OpDet");
    if(fMakeAllPhotonsTree) fThePhotonTreeAll->Write();
    if(fMakeDetectedPhotonsTree) fThePhotonTreeDetected->Write();
    if(fMakeOpDetsTree) fTheOpDetTree->Write();
    if(fMakeOpDetEventsTree) fTheEventTree->Write();
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
    cout << "fNOpDet: " << fGeom->NOpDets() << endl;
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

    if (fProcessMCtruth)  processMC(event);
    processRaw(event);
    if (fProcessCalib) processCalib(event);
    if (fProcessHits)  processHits(event);
    if (fProcessReco) processRecoTracks(event);
    if (fProcessOpDet) processOpDet(event);
    processTiming(event);
    fEventTree->Fill();
    printEvent();

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

    fMC_trackPosition->Clear();
    fMC_trackMomentum->Clear();

    for (int i=0; i<MAX_HITS; i++) {
        hit_channel[i] = 0;
        hit_peakT[i] = 0;
        hit_charge[i] = 0;
	hit_wireID[i] = 0;
	hit_plane[i] = 0;
	hit_tpc[i] = 0;
    }

    fReco_trackPosition->Clear();

    averageWaveforms->Clear();
    waveformCount.clear();
    for (size_t i = 0; i < MAX_OPDET; ++i) {
      OpChannelToOpDet.at(i).clear();
      timestamp.at(i).clear();
    }

    fCount=0;
    fCountEventAll=0;
    fCountEventDetected=0;
    fOpChannel=0;
    for (int i=0; i<MAX_OPDET; i++) {
        fCountOpDetAll[i]=0;
        fCountOpDetDetected[i]=0;
    }
}


//-----------------------------------------------------------------------
void CTree35t::processRaw( const art::Event& event )
{

    //    unsigned int tpcNo, cryoNo;
    // get the objects holding all of the raw data information
    art::Handle< std::vector<raw::RawDigit> > rawdigit;
    event.getByLabel(fRawDigitLabel, rawdigit);
    std::cout << "raw digit label check: " << fRawDigitLabel << std::endl;
 
    // put it in a more easily usable form
    std::vector< art::Ptr<raw::RawDigit> >  rawhits;
    art::fill_ptr_vector(rawhits, rawdigit);

    // rawhits size should == Nchannels == 2048; (no hit channel has a flat 0-waveform)
    cout << "\n Raw Hits size: " << rawhits.size() << endl;

    //loop through all RawDigits (over entire channels)
    fRaw_Nhit = 0;
    for (auto const& hit : rawhits) {      
        int chanId = hit->Channel();
        int nSamples = hit->Samples();
        std::vector<short> uncompressed(nSamples);
	int pedestal = (int)hit->GetPedestal();
	//	std::cout << " channel " << chanId << " pedestal " << pedestal << std::endl;
	// uncompress the data
        if (fUncompressWithPed){
          raw::Uncompress(hit->ADCs(), uncompressed, pedestal, hit->Compression());
        }
        else{
          raw::Uncompress(hit->ADCs(), uncompressed, hit->Compression());
        }
	
	//        short thresh = pedestal + 1; // threshold set to 1 adc;
        short thresh = 1; // threshold set to 1 adc;
        bool isHit = false;
        for (auto const& adc : uncompressed) {
	  short ladc = adc-pedestal;
	  if (ladc > thresh) {
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
        // std::cout << " channel " << chanId << std::endl;
        bool hasTime = false;
        for (int i=0; i<nSamples; i++) {
	  short adc = uncompressed[i]-pedestal;
	  if (adc != 0) {
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
	std::cout<<"i = "<<i
		 <<"\nfMC_id[i]     = "<<fMC_id[i]
		 <<"\nfMC_pdg[i]    = "<<fMC_pdg[i]
		 <<"\nfMC_mother[i] = "<<fMC_mother[i]<<std::endl;
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

        TClonesArray *Lposition = new TClonesArray("TLorentzVector", numberTrajectoryPoints); 
        TClonesArray *Lmomentum = new TClonesArray("TLorentzVector", numberTrajectoryPoints);
        // Read the position and momentum along this particle track
        for(unsigned int j=0; j<numberTrajectoryPoints; j++) {
            new ((*Lposition)[j]) TLorentzVector(particle->Position(j));
            new ((*Lmomentum)[j]) TLorentzVector(particle->Momentum(j));
          // TLorentzVector *pos = new ((*Lposition)[j]) TLorentzVector(/*particle->Position(j)*/);
          // TLorentzVector *mom = new ((*Lmomentum)[j]) TLorentzVector(/*particle->Momentum(j)*/);
          // const TLorentzVector& tmp_pos = particle->Position(j);
          // const TLorentzVector& tmp_mom = particle->Momentum(j);
          // *pos = tmp_pos;
          // *mom = tmp_mom;
        }
        fMC_trackPosition->Add(Lposition);
        fMC_trackMomentum->Add(Lmomentum);


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
        hit_charge[i] = hit->Integral();
        hit_peakT[i] = hit->PeakTime();        
	hit_wireID[i] = hit->WireID().Wire;
	hit_tpc[i] = hit->WireID().TPC;
	hit_plane[i] = hit->WireID().Plane;
    }

    art::Handle< std::vector<recob::Hit> > tHitListHandle;
    if(! event.getByLabel("dcheat", tHitListHandle)) return;
    std::vector< art::Ptr<recob::Hit> > tHitlist;
    art::fill_ptr_vector(tHitlist, tHitListHandle);
    nthits = tHitlist.size();
    if (nthits>MAX_HITS) {
        mf::LogError("CTree35t") << "Event has " << nthits
        << " hits, MAX ALLOWED: " << MAX_HITS;
    }
    for (int i=0; i<nthits; i++) {
        art::Ptr<recob::Hit> hit = tHitlist[i];
        thit_channel[i] = hit->Channel();
        thit_charge[i] = hit->Integral();
        thit_peakT[i] = hit->PeakTime();        
	thit_wireID[i] = hit->WireID().Wire;
	thit_tpc[i] = hit->WireID().TPC;
	thit_plane[i] = hit->WireID().Plane;
    }


    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    if (! event.getByLabel(fClusterModuleLabel, clusterListHandle)) return;
    std::vector<art::Ptr<recob::Cluster> > clusterlist;
    art::fill_ptr_vector(clusterlist, clusterListHandle);
    art::FindManyP<recob::Hit> fm(clusterListHandle, event, fClusterModuleLabel);
    //if (fm.empty()) {cout<<"can't associate cluters with hits. aborting..."<<endl; return;}
    int nclusters = clusterlist.size();
    cout<<"# of clusters is "<<nclusters<<endl;
    nclhits=0;
    for(size_t i=0; i<clusterlist.size(); ++i){
      geo::PlaneID planeid = clusterlist[i]->Plane();
      std::vector< art::Ptr<recob::Hit> > hitslist = fm.at(i);
      size_t n_hits = hitslist.size();
      //for (auto theHit=hitslist.begin();theHit!=hitslist.end();theHit++){
      for(size_t j=0; j<n_hits; j++) {
	art::Ptr<recob::Hit> hit = hitslist[j];
	chit_channel[nclhits] = hit->Channel();
	chit_wire[nclhits] = hit->WireID().Wire;
	chit_peakT[nclhits] = hit->PeakTime();
	chit_plane[nclhits] = planeid.Plane;
	chit_tpc[nclhits] = planeid.TPC;
	chit_cryostat[nclhits] = planeid.Cryostat;
	chit_charge[nclhits] = clusterlist[i]->Charge(0);
	chit_cluster[nclhits] = i;
	nclhits++;
      }
      //cout<<"nclhits = "<<nclhits<<endl;
    }
}

void CTree35t::processRecoTracks( const art::Event& event )
{
    art::Handle< std::vector<recob::Track> > trackListHandle;
    if (! event.getByLabel(fTrackModuleLabel, trackListHandle)) return;
    std::vector<art::Ptr<recob::Track> > tracklist;
    art::fill_ptr_vector(tracklist, trackListHandle);

    reco_nTrack = tracklist.size();
    for (int i=0; i<reco_nTrack; i++) {
        art::Ptr<recob::Track> track = tracklist[i];
        if (track->NumberTrajectoryPoints() > 0) {
            // cout << "track momentum: " << track->VertexMomentum() << endl;
        }
        size_t numberTrajectoryPoints = track->NumberTrajectoryPoints();
        // cout << numberTrajectoryPoints << ", " << track->NumberTrajectoryPoints() << endl;
        
        TClonesArray *Lposition = new TClonesArray("TLorentzVector", numberTrajectoryPoints); 
        // Read the position and momentum along this particle track
        for(unsigned int j=0; j<numberTrajectoryPoints; j++) {
            new ((*Lposition)[j]) TLorentzVector(track->LocationAtPoint(j), 0);
        }
        fReco_trackPosition->Add(Lposition);
    }
}
  
void CTree35t::processOpDet(const art::Event& event)
{
    art::Handle< std::vector<raw::OpDetWaveform> > waveformHandle;
    event.getByLabel(fOpDetInputModule, fInstanceName, waveformHandle);
    art::Handle< std::vector< recob::OpHit > > OpHitHandle;
    event.getByLabel(fOpHitModule, OpHitHandle);
    std::cout << "OpHitHandle->size() = " << OpHitHandle->size() << std::endl;
    std::cout<<"waveformHandle->size() = "<<waveformHandle->size()<<std::endl;
    for (size_t i = 0; i < waveformHandle->size(); ++i) {
      art::Ptr<raw::OpDetWaveform> waveformPtr(waveformHandle, i);
      raw::OpDetWaveform pulse = *waveformPtr;
      double extractedTimestamp = pulse.TimeStamp();
      int channel = pulse.ChannelNumber();
      int OpDetNumber = fGeom->GeometryCore::OpDetFromOpChannel(channel);
      std::cout<<i<<": channel = "<<channel<<", timestamp = "<<std::setprecision(15)<<extractedTimestamp<<", OpDetNumber = "<<OpDetNumber<<", #Ticks = "<<pulse.size()<<", fSampleFreq = "<<fSampleFreq<<std::endl;
      //auto findchannel = find(OpChannelToOpDet.at(OpDetNumber).begin(), OpChannelToOpDet.at(OpDetNumber).end(), channel);
      //if (findchannel == OpChannelToOpDet.at(OpDetNumber).end()) {
      OpChannelToOpDet.at(OpDetNumber).push_back(channel);
      timestamp.at(OpDetNumber).push_back((int)extractedTimestamp);
	//}
      TH1D *hwaveform = new TH1D(Form("avgwaveform_channel_%i_%i", channel, waveformCount[channel]), Form("average waveform channel %i %i-waveform", channel, waveformCount[channel]), pulse.size(), extractedTimestamp, double(pulse.size())/fSampleFreq+extractedTimestamp);
      for (size_t tick = 0; tick < pulse.size(); ++tick) {
	hwaveform->Fill(extractedTimestamp+double(tick)/fSampleFreq, pulse[tick]);
      }
      waveformCount[channel]++;
      averageWaveforms->Add(hwaveform);

    }
    /*
      for (size_t i = 0; i != OpHitHandle->size(); ++i) {
      int channel = OpHitHandle->at(i).OpChannel();
      OpHitCount[channel]++;
      OpHitCountPerEvent[channel]++;
      }*/    
    
    /*
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    bool fUseLitePhotons = lgp->UseLitePhotons();

    // Service for determining opdet responses
    art::ServiceHandle<opdet::OpDetResponseInterface> odresponse;

    if (!fUseLitePhotons) {
      std::cout<<"Cannot get SimPhotonsLite. Aborting..."<<std::endl;
      return;
    }

    art::Handle< std::vector<sim::SimPhotonsLite> > photonHandle;
    event.getByLabel("largeant", photonHandle);
    std::vector<art::Ptr<sim::SimPhotonsLite> > photonlite;
    art::fill_ptr_vector(photonlite, photonHandle);

    // reset counters
    fCountEventAll=0;
    fCountEventDetected=0;
    //std::cout<<"Found OpDet hit collection of size "<< photonlite.size()<<std::endl;

    if (photonlite.size()>0) {
        int size = photonlite.size();
	//std::cout<<"photonlite.size() = "<<size<<std::endl;
        for (int k=0; k<size; k++) {
            //Get data from HitCollection entry
	    art::Ptr<sim::SimPhotonsLite> photon = photonlite[k];
	    fOpChannel=photon->OpChannel;
	    std::map<int, int> PhotonsMap = photon->DetectedPhotons;
	    
	    for (auto it = PhotonsMap.begin(); it!= PhotonsMap.end(); it++){
	        // Calculate wavelength in nm
	        fWavelength= 128;
		//Get arrival time from phot 
		fTime= it->first*2;
		//std::cout<<"Arrival time: " << fTime<<std::endl;
		for (int i = 0; i < it->second ; i++) { // second: number of photons
		    // Increment per OpDet counters and fill per phot trees      
		    fCountOpDetAll[fOpChannel]++;
		    if (fMakeAllPhotonsTree) fThePhotonTreeAll->Fill();
		    if (odresponse->detectedLite(fOpChannel)) {
		        if (fMakeDetectedPhotonsTree) fThePhotonTreeDetected->Fill();
			fCountOpDetDetected[fOpChannel]++;
			//std::cout<<"OpDetResponseInterface PerPhoton : Event "<<fEvent<<" OpChannel " <<fOpChannel << " Wavelength " << fWavelength << " Detected 1 "<<std::endl;
		    } else {
		        //std::cout<<"OpDetResponseInterface PerPhoton : Event "<<fEvent<<" OpChannel " <<fOpChannel << " Wavelength " << fWavelength << " Detected 0 "<<std::endl;
		    }
		}
	    }
	    // Incremenent per event and fill Per OpDet trees
	    fCountEventAll+=fCountOpDetAll[fOpChannel];
	    fCountEventDetected+=fCountOpDetDetected[fOpChannel];
	    //std::cout<<"OpDetResponseInterface PerOpDet : Event "<<fEvent<<" OpDet " << fOpChannel << " All " << fCountOpDetAll << " Det " <<fCountOpDetDetected<<std::endl;
	}
	if(fMakeOpDetsTree) fTheOpDetTree->Fill();
	if (fMakeOpDetEventsTree) fTheEventTree->Fill();
	//std::cout<<"OpDetResponseInterface PerEvent : Event "<<fEvent<<" All " << fCountEventAll << " Det " <<fCountEventDetected<<std::endl; 
    } else { 
        // if empty OpDet hit collection, 
        // add an empty record to the per event tree 
        if(fMakeOpDetsTree) fTheOpDetTree->Fill();
        if(fMakeOpDetEventsTree) fTheEventTree->Fill();
    }
    */
}

// -------------------------------------------------------------------
void CTree35t::processTiming(const art::Event& event)
{
  // RCE
  art::Handle<artdaq::Fragments> RCErawFragments;
  event.getByLabel(fRCERawDataLabel, fRCEFragType, RCErawFragments);
  bool RCEPresent = true;
  try { RCErawFragments->size(); }
  catch(std::exception e) {
    std::cout << "WARNING: Raw RCE data not found in event " << event.event() << std::endl;
    RCEPresent = false;
  }

  if (RCEPresent) {
    if(!RCErawFragments.isValid()){
      std::cerr << "Run: " << event.run() << ", SubRun: " << event.subRun()	<< ", Event: " << event.event() << " is NOT VALID" << std::endl;
      throw cet::exception("RCErawFragments NOT VALID");
      return;
    }
    
    artdaq::Fragments const& rawFragmentsRCE = *RCErawFragments;
    std::cout<<"rawFragmentsRCE.size() = "<<rawFragmentsRCE.size()<<std::endl;
    for ( size_t fragIndex=0; fragIndex < rawFragmentsRCE.size(); ++fragIndex ) {
      int ThisADCcount = 0;
      const artdaq::Fragment &singleFragment = rawFragmentsRCE[fragIndex];
      lbne::TpcMilliSliceFragment millisliceFragment(singleFragment);
      auto numMicroSlices = millisliceFragment.microSliceCount();
      for (unsigned int i_micro=0;i_micro<numMicroSlices;i_micro++) { // Loop through all MicroSlices
	std::unique_ptr <const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(i_micro);
	auto numNanoSlices = microSlice->nanoSliceCount();
	if (numNanoSlices) {
	  //ConsistRCE = 1;
	  ThisADCcount += numNanoSlices;
	  std::cout<<fragIndex<<": "<<numNanoSlices<<", "<<ThisADCcount<<std::endl;
	} // If numNanoSlices.                                                                    
	if ( fragIndex==0 && i_micro==0 ) {
	  std::cout << "Getting RCE time for event " << event.event();                        
	  if ( microSlice->nanoSliceCount() ) { // If this MicroSlice has any NanoSlices      
	    RCETimeBegin = microSlice->nanoSlice(0)->nova_timestamp();
	    std::cout << ", taking RCETime from first nanoslice " << RCETimeBegin << std::endl;      
	  } // NanoSlice                                                    
	  else {
	    RCETimeBegin = microSlice->software_message();
	    std::cout << ", taking RCETime from header " << RCETimeBegin << std::endl;  
	  }	  
	} // Looking at first MicroSlice of first Fragment                                 
	if ( fragIndex==rawFragmentsRCE.size()-1 && i_micro==numMicroSlices-1) { // last microslice of last fragment
	  RCETimeEnd = microSlice->nanoSlice(numNanoSlices-1)->nova_timestamp();
	} else {
	  std::unique_ptr <const lbne::TpcMicroSlice> prevMicroSlice = millisliceFragment.microSlice(i_micro-1);
	  RCETimeEnd = 2.*microSlice->software_message() - prevMicroSlice->software_message();
	}
      } // MicroSlice                                                                                                     
    }
  }

  // SSP Section ------------------------
  bool SSPPresent = true;
  art::Handle<artdaq::Fragments> SSPrawFragments;
  event.getByLabel(fSSPRawDataLabel, fSSPFragType, SSPrawFragments);
  try { SSPrawFragments->size(); }
  catch(std::exception e) {
    mf::LogWarning("SSPToOffline") << "WARNING: Raw SSP data not found in event " << event.event();
    SSPPresent = false;
  }

  if (SSPPresent) {
    if(!SSPrawFragments.isValid()){
      mf::LogError("SSPToOffline") << "Run: " << event.run() << ", SubRun: " << event.subRun() << ", Event: " << event.event() << " is NOT VALID";
      throw cet::exception("raw NOT VALID");
      return;
    }
    
    artdaq::Fragments const& rawFragmentsSSP = *SSPrawFragments;
    std::cout<<"rawFragmentsSSP.size() = "<<rawFragmentsSSP.size()<<std::endl;

    if ( rawFragmentsSSP.size() ) {
      const auto& frag(rawFragmentsSSP[0]);
      const SSPDAQ::MillisliceHeader* meta=0;
      if(frag.hasMetadata())
	{
	  meta = &(frag.metadata<lbne::SSPFragment::Metadata>()->sliceHeader);	  
	  std::cout << "=== SSP Metadata, Start time " << meta->startTime << ", End time " << meta->endTime << " Packet length " << meta->length << " Number of triggers " << meta->nTriggers << "===" << std::endl;
	  SSPTimeBegin = meta->startTime;
	}
    }
    std::cout << "Got SSP start time, it is " << SSPTimeBegin << std::endl;
    
    for (size_t fragIndex = 0; fragIndex < rawFragmentsSSP.size(); ++fragIndex) {
      const auto& frag(rawFragmentsSSP.at(fragIndex));
      const SSPDAQ::MillisliceHeader *meta = 0;
      if (frag.hasMetadata()) {
	meta = &(frag.metadata<lbne::SSPFragment::Metadata>()->sliceHeader);
	std::cout<<fragIndex<<", "<<meta->endTime<<std::endl;
	if (fragIndex == rawFragmentsSSP.size()-1) SSPTimeEnd = meta->endTime;
      }
    }
  } // SSP Present
}

// -------------------------------------------------------------------
void CTree35t::printEvent()
{
    cout << " Event: " << fEvent << endl;
    cout << "   Run: " << fRun << endl;
    cout << "SubRun: " << fSubRun << endl;

    cout << " Hit Channels: " << fRaw_Nhit << endl;
    cout << " MC  Ntracks:" << fMC_Ntrack;
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
    cout << "Number of reco tracks: " << reco_nTrack << endl;
    cout << "Number of all photons:      "<< fCountEventAll << endl;
    cout << "Number of photons detected: "<< fCountEventDetected << endl;

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
