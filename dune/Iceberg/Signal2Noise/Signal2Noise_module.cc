////////////////////////////////////////////////////////////////////////
// Class:       Signal2Noise
// Plugin Type: analyzer (art v3_05_00)
// File:        Signal2Noise_module.cc
//
// Generated at Mon May  4 10:06:28 2020 by Wanwei Wu using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"

using namespace std;


const int kMaxTracks = 1000;
const int kMaxHits = 2000;

class Signal2Noise;


class Signal2Noise : public art::EDAnalyzer {
public:
  explicit Signal2Noise(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Signal2Noise(Signal2Noise const&) = delete;
  Signal2Noise(Signal2Noise&&) = delete;
  Signal2Noise& operator=(Signal2Noise const&) = delete;
  Signal2Noise& operator=(Signal2Noise&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;

private:

  // Declare member data here.
  std::string fHitModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  bool fSaveRawWaveForm;
  std::vector<int> fSelectedWires;

  // reset
  void reset();

  // TTree
  TTree *fEventTree;

  // event information
  int event;
  int run;
  int subrun;
  double evttime;
  int year_month_date;
  int hour_min_sec;
  
  // track
  int ntrks;
  int trkid[kMaxTracks];
  float trkstart[kMaxTracks][3];
  float trkend[kMaxTracks][3];
  float trklen[kMaxTracks];
  float trkthetaxz[kMaxTracks];
  float trkthetayz[kMaxTracks];
  float trkstartcosxyz[kMaxTracks][3];
  float trkendcosxyz[kMaxTracks][3];

  // lifetime
  int ntrkhits[kMaxTracks][3];
  float trkdqdx[kMaxTracks][3][kMaxHits];
  //float trkdedx[kMaxTracks][3][kMaxHits];
  float trkx[kMaxTracks][3][kMaxHits];
  float trkt[kMaxTracks][3][kMaxHits];

  // hit
  double trkhitx[kMaxHits][3][kMaxHits];
  double trkhity[kMaxHits][3][kMaxHits];
  double trkhitz[kMaxHits][3][kMaxHits];

  int wireid[kMaxTracks][kMaxHits];
  int chid[kMaxTracks][kMaxHits];
  int tpcid[kMaxTracks][kMaxHits];
  float hit_plane[kMaxTracks][kMaxHits];

  double cosgma[kMaxTracks][kMaxHits];

  float amp[kMaxTracks][kMaxHits];
  int tamp[kMaxTracks][kMaxHits];
  float ped[kMaxTracks][kMaxHits];

  float noiserms[kMaxTracks][kMaxHits]; // calculate directly
  float noisermsfit[kMaxTracks][kMaxHits]; // rms from gaus fit
  float fsigma_ped[kMaxTracks][kMaxHits]; // fit for ped and sigma
  float fped[kMaxTracks][kMaxHits];
  float fchi2[kMaxTracks][kMaxHits];
  int fndf[kMaxTracks][kMaxHits];

  // waveform
  int fNticks; 
  int fNticksReadout;
  float fSampleRate;
  TH1I* fWaveForm[10];
};


Signal2Noise::Signal2Noise(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fHitModuleLabel         = p.get<std::string>("HitModuleLabel");
  fTrackModuleLabel       = p.get<std::string>("TrackModuleLabel");
  fCalorimetryModuleLabel = p.get<std::string>("CalorimetryModuleLabel");
  fSaveRawWaveForm        = p.get<bool>("SaveRawWaveForm");
  fSelectedWires          = p.get<std::vector<int>>("SelectedWires");

  // DetectorPropertiesService
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
  fNticks = detProp.NumberTimeSamples(); // number of clock ticks per event
  fNticksReadout = detProp.ReadOutWindowSize(); // number of clock ticks per readout window
  fSampleRate = sampling_rate(clockData); // period of the TPC readout electronics clock
  cout << "Numer of clock ticks per event: " << fNticks << endl;
  cout << "Numer of clock ticks per readout window: " << fNticksReadout << endl;
  cout << "Sampling rate: " << fSampleRate << endl;
}


void Signal2Noise::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  reset();

  // event info
  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  art::Timestamp ts = e.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();
  
  //cout << "tts " << tts << endl;
  //cout << "evttime " << evttime << endl;
  
  UInt_t year=0;
  UInt_t month=0;
  UInt_t day=0;

  year_month_date=tts.GetDate(kTRUE,0,&year,&month,&day);

  UInt_t hour=0;
  UInt_t min=0;
  UInt_t sec=0;

  hour_min_sec=tts.GetTime(kTRUE,0,&hour,&min,&sec);

  // channel status
  lariov::ChannelStatusProvider const& channelStatus
    = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

  // waveform
  art::Handle< std::vector<raw::RawDigit> > rawdigitListHandle;
  std::vector< art::Ptr<raw::RawDigit> > rawdigitlist;
  if (e.getByLabel("tpcrawdecoder", "daq", rawdigitListHandle)) {
    art::fill_ptr_vector(rawdigitlist, rawdigitListHandle);
  }

  //cout << "rawdigitlist.size():  " << rawdigitlist.size() << endl;

  // hit
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector< art::Ptr<recob::Hit> > hitlist;
  if (e.getByLabel(fHitModuleLabel, hitListHandle)) {
    art::fill_ptr_vector(hitlist, hitListHandle);
  }

  // track
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector< art::Ptr<recob::Track> > tracklist;
  if (e.getByLabel(fTrackModuleLabel,trackListHandle)) {
    art::fill_ptr_vector(tracklist, trackListHandle);
  }
 
  /*
  // wire
  art::Handle< std::vector<recob::Wire> > wireListHandle;
  std::vector< art::Ptr<recob::Track> > wirelist;
  if (e.getByLabel("digitwire",wireListHandle)) {
    art::fill_ptr_vector(wirelist, wireListHandle);
  }
  */

  art::ServiceHandle<geo::Geometry> geom;

  art::FindManyP<recob::Hit> fmtrkhit(trackListHandle, e, fTrackModuleLabel);
  art::FindManyP<anab::Calorimetry> fmtrkcalo(trackListHandle, e, fCalorimetryModuleLabel);
  //art::FindManyP<raw::RawDigit> fmhitrawdigit(hitListHandle, e, fHitModuleLabel);
  //art::FindManyP<recob::Wire> fmhitwire(hitListHandle, e, fHitModuleLabel);
  //art::FindManyP<raw::RawDigit> fmwirerawdigit(wireListHandle, e, "digitwire");
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmhittrkmeta(trackListHandle, e, fTrackModuleLabel);

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);

  // waveform: save several waveforms for check
  int nwaveform = 0; // only save 10 waveforms
  int nwaveform_plane_0 = 0; // only save 3 waveforms
  int nwaveform_plane_1 = 0; // only save 3 waveforms
  int nwaveform_plane_2 = 0; // only save 3 waveforms

  ntrks = 0;
  for (const auto& trk : tracklist) {

    trkid[ntrks] = trk->ID();
    trklen[ntrks] = trk->Length();

    // TVector3
    auto start = trk->Vertex();
    auto end = trk->End();
    auto start_dir = trk->VertexDirection();
    auto end_dir = trk->EndDirection();

    trkthetaxz[ntrks] = std::atan2(start_dir.X(), start_dir.Z());
    trkthetayz[ntrks]= std::atan2(start_dir.Y(), start_dir.Z());
    
    trkstart[ntrks][0] = start.X();
    trkstart[ntrks][1] = start.Y();
    trkstart[ntrks][2] = start.Z();

    trkend[ntrks][0] = end.X();
    trkend[ntrks][1] = end.Y();
    trkend[ntrks][2] = end.Z();

    trkstartcosxyz[ntrks][0] = start_dir.X();
    trkstartcosxyz[ntrks][1] = start_dir.Y();
    trkstartcosxyz[ntrks][2] = start_dir.Z();

    trkendcosxyz[ntrks][0] = end_dir.X();
    trkendcosxyz[ntrks][1] = end_dir.Y();
    trkendcosxyz[ntrks][2] = end_dir.Z();
    
    // calometry
    if (fmtrkcalo.isValid()) {
      std::vector<art::Ptr<anab::Calorimetry>> calos = fmtrkcalo.at(ntrks);
      for (size_t icalo=0; icalo<calos.size(); icalo++) {
        if (!calos[icalo]) continue;
        if (!calos[icalo]->PlaneID().isValid) continue;
        int planenum = calos[icalo]->PlaneID().Plane;
        if (planenum<0 || planenum>2) continue;

        const size_t NHits = calos[icalo] -> dEdx().size();
        ntrkhits[ntrks][planenum] = NHits;

        double minx = 1e9;
        for (size_t iHit=0; iHit<NHits; ++iHit) {
          cout << "plane: "<< planenum << "; pitch: " << (calos[icalo]->TrkPitchVec())[iHit] << endl;
          if ((calos[icalo]->TrkPitchVec())[iHit]>1) continue;
          const auto& TrkPos = (calos[icalo] -> XYZ())[iHit];
          if (TrkPos.X()<minx)
            minx = TrkPos.X();
        }// loop NHits

        for(size_t iHit = 0; iHit < NHits; ++iHit) {
          if ((calos[icalo]->TrkPitchVec())[iHit]>1) continue;
          const auto& TrkPos1 = (calos[icalo] -> XYZ())[iHit];
          double x = TrkPos1.X()-minx; //subtract the minx to get correct t0
          double XDriftVelocity = detProp.DriftVelocity()*1e-3; //cm/ns
          double t = x/(XDriftVelocity*1000); //change the velocity units to cm/ns to cm/us
          trkx[ntrks][planenum][iHit] = x;
          trkt[ntrks][planenum][iHit] = t;
          trkdqdx[ntrks][planenum][iHit] = (calos[icalo]->dQdx())[iHit];
        }

      } // end loop over icalo
    
    } //  end if fmtrkcalo

    // hits associated with trk
    std::vector<art::Ptr<recob::Hit>> allhits = fmtrkhit.at(ntrks); // use either trk.key() or ntrks 
    
    for (size_t ihit=0; ihit<allhits.size(); ihit++) {
      // wire plane information
      unsigned int wireplane = allhits[ihit]->WireID().Plane;
      if (wireplane <0 || wireplane>2) continue;
      unsigned int wire = allhits[ihit]->WireID().Wire;
      unsigned int tpc = allhits[ihit]->WireID().TPC;
      unsigned int channel = allhits[ihit]->Channel();
       
      if (channelStatus.IsBad(channel)) continue;

      // hit position: not all hits are associated with space points, using neighboring space points to interpolate
      double xyz[3] = {-9999.0, -9999.0, -9999.0};
      bool fBadhit = false;

      if (fmhittrkmeta.isValid()) {
        auto vhit = fmhittrkmeta.at(ntrks);
        auto vmeta = fmhittrkmeta.data(ntrks);

        for (size_t ii=0; ii<vhit.size(); ii++) {
          if (vhit[ii].key() == allhits[ihit].key()) {
            
            if (vmeta[ii]->Index() == std::numeric_limits<int>::max()) {
              fBadhit = true;
              continue;
            }

            if (vmeta[ii]->Index() >= trk->NumberTrajectoryPoints()) {
              throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<trk->NumberTrajectoryPoints()<<" for track index "<<ntrks<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
            }

            if (!trk->HasValidPoint(vmeta[ii]->Index())) {
              fBadhit = true;
              continue;
            }
            
            auto loc = trk->LocationAtPoint(vmeta[ii]->Index());
            xyz[0] = loc.X();
            xyz[1] = loc.Y();
            xyz[2] = loc.Z();
            
            break;
          } // if (vhit[ii].key() == allhits[ihit].key())
        
        } // end of for ii

      } // if fmhittrkmeta.isValid()

      if (fBadhit) continue;
      
      trkhitx[ntrks][wireplane][ihit] = xyz[0];
      trkhity[ntrks][wireplane][ihit] = xyz[1];
      trkhitz[ntrks][wireplane][ihit] = xyz[2];

      wireid[ntrks][ihit] = wire;
      chid[ntrks][ihit] = channel;
      tpcid[ntrks][ihit] = tpc;
      hit_plane[ntrks][ihit] = wireplane;

      // calculate track angle w.r.t. wire
      double angleToVert = geom->WireAngleToVertical(geom->View(allhits[ihit]->WireID()), allhits[ihit]->WireID().TPC, allhits[ihit]->WireID().Cryostat)-0.5*::util::pi<>();
      const auto& dir = trk->DirectionAtPoint(0);
      // angleToVert: return the angle w.r.t y+ axis, anti-closewise
      // dir: 3d track direction: u = (x,y,z);
      // vector that perpendicular to wires in yz plane v = (0, sin(angleToVert), cos(angleToVert))   
      // cos gamma = u.Dot(v)/(u.mag()*v.mag()) here, u.mag()=v.mag()=1
      double tmp_cosgamma = abs(sin(angleToVert)*dir.Y() + std::cos(angleToVert)*dir.Z());
      cosgma[ntrks][ihit] = tmp_cosgamma;

      //cout << "track direction: " << trkstartcosxyz[ntrks][0] << ", " << trkstartcosxyz[ntrks][1] << ", " << trkstartcosxyz[ntrks][2]  << endl;
      //cout << "angleToVert: " << angleToVert << endl;
      //cout << "dir: " << dir.X() << ", " << dir.Y() << ", " << dir.Z() << endl;

      /*
      // check wire direction on each plane
      cout << geom->Plane(wireplane).Wire(wire).ThetaZ(true) << endl;
      double wirestart[3];
      double wireend[3];
      geom->Plane(wireplane).Wire(wire).GetStart(wirestart);
      geom->Plane(wireplane).Wire(wire).GetEnd(wireend);
      cout << "wirestart: (" << wirestart[0] << ", " << wirestart[1] << ", "<<  wirestart[2] << ")" << endl;
      cout << "wireend: (" << wireend[0] << ", " << wireend[1] << ", "<<  wireend[2] << ")" << endl;
      */
      
      // raw waveform
      // loop over rawdigitlist to find the associated rawdighit with same channel of the hit
      int key_rawdigit = -1;
      for (size_t ird=0; ird<rawdigitlist.size(); ird++) {
        if (rawdigitlist[ird]->Channel() == channel) {
          key_rawdigit = ird;
          break;
        }
      }
      
      if (key_rawdigit == -1) continue; // www
      int datasize = rawdigitlist[key_rawdigit]->Samples();
      if (datasize==0) continue; // www: check on this
      // to use a compressed RawDigit, one has to create a new buffer, fill and use it
      //if (ihit==15) cout << "rawdigitlist.size(): " << rawdigitlist.size() << ";  rawdigitlist[key_rawdigit]->Channel(): " << rawdigitlist[key_rawdigit]->Channel() << endl;
      //if (ihit==15) cout << "channel: " << channel << "; tpc: " <<tpc << "; key_rawdigit: " << key_rawdigit<< "; datasize: " << datasize << endl;
      std::vector<short> rawadc(datasize); // create a buffer
      raw::Uncompress(rawdigitlist[key_rawdigit]->ADCs(), rawadc, rawdigitlist[key_rawdigit]->Compression());

      // pedestal
      ped[ntrks][ihit] = rawdigitlist[key_rawdigit]->GetPedestal(); // Pedestal level (ADC counts)

      if (fSaveRawWaveForm && nwaveform<10 && nwaveform_plane_0<4 && nwaveform_plane_1<4 && nwaveform_plane_2<5) {
        if (wireplane==0) nwaveform_plane_0++;
        if (wireplane==1) nwaveform_plane_1++;
        if (wireplane==2) nwaveform_plane_2++;

        fWaveForm[nwaveform]->SetNameTitle(Form("plane_%d_AdcChannel_%d", wireplane,  channel), Form("AdcChannel%d", channel));  
        
        for (int jj=0; jj<datasize; jj++) {
          fWaveForm[nwaveform]->SetBinContent(jj+1, rawadc[jj]-ped[ntrks][ihit]);
        }
        nwaveform++;
      }

      // ROI from the reconstructed hits
      int t0 = allhits[ihit]->PeakTime() - 3*(allhits[ihit]->RMS());
      if (t0<0) t0 = 0;
      int t1 = allhits[ihit]->PeakTime() + 3*(allhits[ihit]->RMS());
      if (t1>= datasize) t1 = datasize - 1;
      //cout << "t0: " << t0 << " ; t1: " << t1 << endl;
      
      // maximum pulse height of raw waveform
      float temp_max_pulseheight = -9999.;
      int temp_t_max_pulseheight = -1; // time in unit of ticks
      for (int itime=t0; itime <=t1; itime++) {
        if (rawadc[itime] > temp_max_pulseheight) {
          temp_max_pulseheight = rawadc[itime];
          temp_t_max_pulseheight = itime;
        }
      }

      amp[ntrks][ihit] = temp_max_pulseheight;
      tamp[ntrks][ihit] = temp_t_max_pulseheight;

      // noise rms calculation: ideally, this should be done for all wires, not only wires that have hits
        //std::cout << " ntrks: "<< ntrks << "; ihit: " << ihit << std::endl;
        //std::cout << "rawadc[0]: " << rawadc[0] << std::endl;
      // method 1: calculate rms directly
      int start_ped = 0; 
      int end_ped = 2000; 
      float temp_sum = 0.;
      int temp_number = 0;
      for (int iped=start_ped; iped<end_ped; iped++) {
        //std::cout << "iped: " << iped << "; ntrks: "<< ntrks << "; ihit: " << ihit << std::endl;
        //std::cout << "rawadc[iped]: " << rawadc[iped] << std::endl;
        //std::cout << "ped[ntrks][ihit]: "<< ped[ntrks][ihit]<<std::endl;
        //if (iped > tamp[ntrks][ihit]-10 && iped < tamp[ntrks][ihit]+10) continue; // ideally we should use this to skip ROI region
        float temp_res = rawadc[iped] - ped[ntrks][ihit];
        if (abs(temp_res) > 8.) continue; // skip ROI with a threshold
        temp_sum += temp_res*temp_res;
        temp_number++;
      }
      noiserms[ntrks][ihit] = sqrt(temp_sum/temp_number);
     
      // method 2: fit noise histogram with a gaus
      int nbin_noise = 40;
      double lbin_noise = -20.;
      double hbin_noise = 20.;
      
      TH1F *h1_noise = new TH1F(TString::Format("noise_trk%d_hit%d",ntrks, (int)ihit), TString::Format("noise_trk%d_hit%d",ntrks, (int)ihit), nbin_noise, lbin_noise, hbin_noise);
      // fill the readout datasize for each wire with hit: signal is included but would not affect the noise rms since signals are far way from the noise peak. One may also exclude signals by using ROI threshold cuts
      for (int jj=0; jj<datasize; jj++) {
        h1_noise->Fill(rawadc[jj]-ped[ntrks][ihit]);
      }
      TF1 *f1_noise = new TF1("f1_noise", "gaus" , -10., 10.);
      double par[3];
      h1_noise->Fit(f1_noise, "WWQ");
      f1_noise->GetParameters(&par[0]);
      noisermsfit[ntrks][ihit] = par[2]; // sigma from gaus fit

      delete h1_noise;
      delete f1_noise;

      // method 3: baseline projection and get sigma of pedestal distribution
      // To use this, one should be careful. pedestals are different for different plane and wires. plane 0, 1: ~2400  plane 2: ~950
      int nadc_ped = 3500; // include sig/noise and pedestal
      int rebin_ped = 10;
      int nbin_ped =nadc_ped/rebin_ped;

      TString name_ped = Form("ped1d_trk%d_hit%d_rebin", ntrks, (int)ihit);
      TH1 *phr = new TH1F(name_ped.Data(), name_ped.Data(), nbin_ped, 0, nadc_ped);
      phr->SetDirectory(0);
      for (int iped=start_ped; iped<datasize; iped++) {
        phr->Fill(rawadc[iped]);
      }
      
      int rbinmax1 = phr->GetMaximumBin();
      double radcmax1 = phr->GetBinCenter(rbinmax1);
      double adcmean = phr->GetMean();
      //cout << "rbinmax1: " << rbinmax1 << " ; radcmax1: " << radcmax1 << " ;  adcmean: " << adcmean << endl;
      
      // max may just be a sticky code. reduce it and find the next maximum
      double tmpval = 0.5*(phr->GetBinContent(rbinmax1-1)+phr->GetBinContent(rbinmax1+1));
      phr->SetBinContent(rbinmax1, tmpval);
      int rbinmax2 = phr->GetMaximumBin();

      // define max to be the first value if the two maxima are close or the average if they are far part
      // average peak position if there are twoo peaks
      int rbinmax = rbinmax1;
      if ( abs(rbinmax2-rbinmax1) > 1 ) {
        rbinmax = (rbinmax1 + rbinmax2)/2;
        adcmean = phr->GetMean();
      }

      double adcmax = phr->GetBinCenter(rbinmax);
      //cout << "adcmax: " << adcmax << " ; adcmean: " << adcmean << endl;

      delete phr;

      double wadc = 100.;
      // make sure the peak bin stays in range
      if ( abs(adcmax-radcmax1) > 0.5*wadc ) adcmax = radcmax1;
      double adc1 = adcmax - 0.5*wadc;
      adc1 = 10*int(adc1/10);
      if ( adcmean > adcmax + 10) adc1 += 10;
      double adc2 = adc1 + wadc;

      TH1* phf = new TH1F(name_ped.Data(), name_ped.Data(), wadc, adc1, adc2);
      phf->SetDirectory(nullptr);
      int countLo = 0;
      int countHi = 0;
      int count = 0;
      for (int isam=0; isam<datasize; ++isam ) {
        if ( isam < adc1 ) ++countLo;
        if ( isam >= adc2 ) ++countHi;
        ++count;
        phf->Fill(rawadc[isam]);
      }

      // Fetch the peak bin and suppress it for the fit if more than 20% (but not all) the data is in it
      int binmax = phf->GetMaximumBin();
      double valmax = phf->GetBinContent(binmax);
      double rangeIntegral = phf->Integral(1, phf->GetNbinsX());
      double peakBinFraction = valmax/rangeIntegral;
      bool allBin = peakBinFraction > 0.99;
      bool dropBin = valmax > 0.2*phf->Integral() && !allBin;
      int nbinsRemoved = 0;
      if ( dropBin ) {
        double tmpval = 0.5*(phf->GetBinContent(binmax-1)+phf->GetBinContent(binmax+1));
        phf->SetBinContent(binmax, tmpval);
        rangeIntegral = phf->Integral(1, phf->GetNbinsX());
        ++nbinsRemoved;
      }

      double amean = phf->GetMean() + 0.5;
      double alim1 = amean - 25.0;
      double alim2 = amean + 25.0;
      double FitRmsMin = 1.;
      double FitRmsMax = 20.;

      TF1 fitter("pedgaus", "gaus", adc1, adc2, TF1::EAddToList::kNo);
      fitter.SetParameters(0.1*rangeIntegral, amean, 5.0);
      fitter.SetParLimits(0, 0.01*rangeIntegral, rangeIntegral);
      fitter.SetParLimits(1, alim1, alim2);
      fitter.SetParLimits(2, FitRmsMin, FitRmsMax);
      if ( allBin ) fitter.FixParameter(1, amean);  // Fix posn.
      string fopt = "0";
      fopt = "WWBQ"; 
      // fit option:
      //   "WW" Set all weights to 1 including empty bins; ignore error bars
      //   "B"  fix one or more parameters and the fitting function is a predefined one, like polN, expo, landau, gaus
      //   "Q" Quiet mode (minimum printing)
      phf->Fit(&fitter, fopt.c_str());
      fsigma_ped[ntrks][ihit] = fitter.GetParameter(2);
      fped[ntrks][ihit] = fitter.GetParameter(1);
      fchi2[ntrks][ihit] = fitter.GetChisquare();
      fndf[ntrks][ihit] = fitter.GetNDF();
      
      /*
      cout <<"amp["<<ntrks<<"]["<<ihit<<"]: "<< amp[ntrks][ihit]<<"\n"
           <<"noiserms["<<ntrks<<"]["<<ihit<<": "<< noiserms[ntrks][ihit] <<"\n"           <<"s/n: " << (amp[ntrks][ihit] - ped[ntrks][ihit])/noiserms[ntrks][ihit] << "\n"
           << "chi2/ndf: " << fchi2[ntrks][ihit]<< "/"<<fndf[ntrks][ihit]<<"\n"
           << "fsigma_ped: " << fsigma_ped[ntrks][ihit]<< "\n"
           << "fped: " << fped[ntrks][ihit] << endl;
      */
      delete phf;

    } // end of for ihit
    
    ++ ntrks;
  } // end of for trk
  fEventTree->Fill();
}


void Signal2Noise::beginJob(){
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event");
  fEventTree->Branch("event", &event, "event/I");
  fEventTree->Branch("run", &run, "run/I");
  fEventTree->Branch("subrun", &subrun, "subrun/I");
  fEventTree->Branch("evttime", &evttime, "evttime/D");
  fEventTree->Branch("year_month_date", &year_month_date, "year_month_date/I");
  fEventTree->Branch("hour_min_sec", &hour_min_sec, "hour_min_sec/I");
  
  // track
  fEventTree->Branch("ntrks", &ntrks, "ntrks/I");
  fEventTree->Branch("trkid", trkid, "trkid[ntrks]/I");
  fEventTree->Branch("trkstart", trkstart, "trkstart[ntrks][3]/F");
  fEventTree->Branch("trkend", trkend, "trkend[ntrks][3]/F");
  fEventTree->Branch("trklen", trklen, "trklen[ntrks]/F");
  fEventTree->Branch("trkthetaxz", trkthetaxz, "trkthetaxz[ntrks]/F");
  fEventTree->Branch("trkthetayz", trkthetayz, "trkthetayz[ntrks]/F");
  fEventTree->Branch("trkstartcosxyz", trkstartcosxyz, "trkstartcosxyz[ntrks][3]/F");
  fEventTree->Branch("trkendcosxyz", trkendcosxyz, "trkendcosxyz[ntrks][3]/F");
 
  fEventTree->Branch("ntrkhits", ntrkhits, "ntrkhits[ntrks][3]/I");
  fEventTree->Branch("trkdqdx", trkdqdx, "trkdqdx[ntrks][3][1000]/F");
  //fEventTree->Branch("trkdedx", trkdedx, "trkdedx[ntrks][3][1000]/F");
  fEventTree->Branch("trkx", trkx, "trkx[ntrks][3][1000]/F");
  fEventTree->Branch("trkt", trkt, "trkt[ntrks][3][1000]/F");
  


  // hit
  fEventTree->Branch("trkhitx", trkhitx, "trkhitx[ntrks][3][1000]/D");
  fEventTree->Branch("trkhity", trkhity, "trkhity[ntrks][3][1000]/D");
  fEventTree->Branch("trkhitz", trkhitz, "trkhitz[ntrks][3][1000]/D");
 
  fEventTree->Branch("wireid", wireid, "wireid[ntrks][1000]/I");
  fEventTree->Branch("chid", wireid, "chid[ntrks][1000]/I");
  fEventTree->Branch("tpcid", wireid, "tpcid[ntrks][1000]/I");

  fEventTree->Branch("hit_plane", hit_plane, "hit_plane[ntrks][1000]/F");
  fEventTree->Branch("ped", ped, "ped[ntrks][1000]/F");
  fEventTree->Branch("amp", amp, "amp[ntrks][1000]/F");
  fEventTree->Branch("tamp", tamp, "tamp[ntrks][1000]/I");
  fEventTree->Branch("cosgma", cosgma, "cosgma[ntrks][1000]/D");

  fEventTree->Branch("noiserms", noiserms, "noiserms[ntrks][1000]/F");
  fEventTree->Branch("noisermsfit", noisermsfit, "noisermsfit[ntrks][1000]/F");
  fEventTree->Branch("fsigma_ped", fsigma_ped, "fsigma_ped[ntrks][1000]/F");
  fEventTree->Branch("fped", fped, "fped[ntrks][1000]/F");
  fEventTree->Branch("fchi2", fchi2, "fchi2[ntrks][1000]/F");
  fEventTree->Branch("fndf", fndf, "fndf[ntrks][1000]/I");



  // waveform
  if (fSaveRawWaveForm) {
    for (int i=0; i<10; i++) {
      fWaveForm[i] = tfs->make<TH1I>(Form("waveform_%d",i), "wire waveform", fNticksReadout, 0, fNticksReadout);
      fWaveForm[i]->SetStats(0);
      fWaveForm[i]->GetXaxis()->SetTitle("Time [ticks]");
      fWaveForm[i]->GetYaxis()->SetTitle("ADC");
    }
  }
}

void Signal2Noise::reset(){
  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  year_month_date = -99999;
  hour_min_sec = -99999;

  
  ntrks = -99999;
  for (size_t i=0; i<kMaxTracks; ++i) {
    trkid[i] = -1;
    trklen[i] = -1.0;
    trkthetaxz[i] = -9999.0;
    trkthetayz[i] = -9999.0;
    for (int j=0; j<3; j++) {
      trkstart[i][j] = -9999.0;
      trkend[i][j] = -9999.0;
      trkstartcosxyz[i][j] = -9999.0;
      trkendcosxyz[i][j] = -9999.0;
      
      ntrkhits[i][j] = -9999;

      for (int k=0; k<kMaxHits; k++) {
        trkhitx[i][j][k] = -9999.0;
        trkhity[i][j][k] = -9999.0;
        trkhitz[i][j][k] = -9999.0;

        trkdqdx[i][j][k] = -9999.0;
        //trkdedx[i][j][k] = -9999.0;
        trkx[i][j][k] = -9999.0;
        trkt[i][j][k] = -9999.0;
      }
    }

    for (int k=0; k<kMaxHits; k++) {
      wireid[i][k] = -99;
      chid[i][k] = -99;
      tpcid[i][k] = -99;

      hit_plane[i][k] = -1;
      ped[i][k] = -9999.0;
      amp[i][k] = -9999.0;
      tamp[i][k] = -1;
      
      noiserms[i][k] = -999.;
      noisermsfit[i][k] = -999.;
      fsigma_ped[i][k] = -99999.;
      fped[i][k] = -99999.;
      fchi2[i][k] = -99999.;
      fndf[i][k] = -99;

      cosgma[i][k] = -99.;
    
    }

  }


}

DEFINE_ART_MODULE(Signal2Noise)
