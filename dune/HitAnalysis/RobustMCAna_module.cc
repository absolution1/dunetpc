////////////////////////////////////////////////////////////////////////
// Class:       RobustMCAna
// Plugin Type: analyzer (art v2_05_00)
// File:        RobustMCAna_module.cc
//
// Generated at Sun Nov 20 14:28:42 2016 by Matthew Thiesse using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/AnalysisBase/T0.h"

#include <algorithm>
#include "TMath.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

//for the counter position map
#include "dune/daqinput35t/PennToOffline.h"

namespace hit {
  class RobustMCAna;
}

class hit::RobustMCAna : public art::EDAnalyzer {
public:
  explicit RobustMCAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RobustMCAna(RobustMCAna const &) = delete;
  RobustMCAna(RobustMCAna &&) = delete;
  RobustMCAna & operator = (RobustMCAna const &) = delete;
  RobustMCAna & operator = (RobustMCAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  raw::ChannelID_t OppositeTPCChannel(raw::ChannelID_t);
  bool ValidTrigger(std::vector<unsigned int> evtTriggers, unsigned int & c1arg, unsigned int & c2arg, unsigned int & trignumarg);
  void NumElectronsTofC(double & input);
  void SumADCTofC(double & input);
  void fCToSumADC(double & input);
  double Purity(double tp, double fp);
  double Efficiency(double tp, double fn);
  double MCC(double tp, double fp, double fn, double tn);
  double VMean(std::vector<double> v);

  // Struct to contain combined IDE information
  struct MCHit {
    raw::ChannelID_t channel;
    int pdg;
    std::vector<sim::IDE> idevec;
    std::vector<unsigned short> tdcvec;
    std::vector<double> chargevec;
    std::vector<const simb::MCParticle *> particlevec;
  };

  // FHICL parameters
  std::string fHitsModuleLabel;
  std::string fCounterT0ModuleLabel;
  bool fVerbose;
  double mcscale;
  int fNbins;
  double preampGain;
  double AmpToAreaFactor;
  double ADCtomV;

  // Histograms in the TFile
  TH1D* hEventPurity; // aka Precision
  TH1D* hEventEfficiency; // aka Recall
  TH1D* hEventChargePurity;
  TH1D* hEventChargeEfficiency;
  TH1D* hHitChargeRatio;
  TH1D* hdQ;
  TH2D* hChargeComp;
  TH1D* hdT;
  TH1D* hMCC;

  // Output TTrees
  TTree * fEventTree;
  int run;
  int subrun;
  int event;
  unsigned int c1;
  unsigned int c2;
  unsigned int trignum;
  double eventpurity;
  double eventefficiency;
  double chargepurity;
  double chargeefficiency;
  double eventchargeratio;
  double eventRecoQ;
  double eventMCQ;
  double eventdQ;
  double eventdT;
  double eventMCC;
  int tp;
  int fp;
  int fn;
  int tn;
  double tpc;
  double fpc;
  double fnc;

  TTree * fHitTree;
  double RecoQ;
  double MCQ;
  double RecoT;
  double MCT;
  bool foundBoth;

  // Services needed
  detinfo::DetectorProperties const * fDetProp;
  detinfo::DetectorClocks const * fClks;
  const lariov::ChannelStatusProvider* fCSP;

  // Counter position map, as taken from DAQToOffline
  std::map<unsigned int, std::pair<TVector3, std::vector<TVector3> > > fCounterPositionMap;
};

hit::RobustMCAna::RobustMCAna(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fClks(lar::providerFrom<detinfo::DetectorClocksService>())
{
  // Get counter position map for consistent ordering of c1 and c2 numbers
  DAQToOffline::MakeCounterPositionMap("","counterInformation.txt",fCounterPositionMap,0,0,0);

  // Read FHICL parameters
  preampGain = p.get<double>("PreampGainSetting",14.0); // mV/fC
  AmpToAreaFactor = p.get<double>("AmpToAreaFactor",7.615); // convert amplitude of shaping curve to area of pseudo-gaussian
  ADCtomV = p.get<double>("ADCtomV",2.808); // ADC/mV, digitization
  mcscale = p.get<double>("PreviousMCScale");
  fHitsModuleLabel = p.get<std::string>("HitsModuleLabel","robusthit");
  fVerbose = p.get<bool>("Verbose",false);
  fNbins = p.get<int>("NumberBins",22);
  fCounterT0ModuleLabel = p.get<std::string>("CounterT0ModuleLabel","t0counter");

  // Initialize channel status service
  fCSP = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

  // Initialize histograms
  art::ServiceHandle<art::TFileService> fTfs;
  hEventPurity = fTfs->make<TH1D>("hEventPurity","Purity of Reco Hit Selection; TP/(TP+FP); # Events",200,0,1.1);
  hEventEfficiency = fTfs->make<TH1D>("hEventEfficiency","Efficiency of Reco Hit Selection; TP/(TP+FN); # Events",200,0,1.1);
  hEventChargePurity = fTfs->make<TH1D>("hEventChargePurity","Purity of Reco Hit Charge Selection; ChgTP/(ChgTP+ChgFP); # Events",200,0,1.1);
  hEventChargeEfficiency = fTfs->make<TH1D>("hEventChargeEfficiency","Efficiency of Reco Hit Charge Selection; ChgTP/(ChgTP+ChgFN); # Events",200,0,1.1);
  hHitChargeRatio = fTfs->make<TH1D>("hHitChargeRatio","Ratio of reco hit charge to MC hit charge; Reco Charge / MC Charge; # Events",200,0,3);
  hdQ = fTfs->make<TH1D>("hdQ","; Charge Residual, RecoQ-MCQ; # Events",200,-3000,3000);
  hChargeComp = fTfs->make<TH2D>("hChargeComp","Comparison of reco and MC charges; MC Charge (fC); Reco Charge (fC)",400,0,30000,400,0,30000);
  hdT = fTfs->make<TH1D>("hdT","#DeltaT between MC hit and Reco Hit; Reco Hit Time - MC Hit Time (TDC ticks); # Events",200,-30,30);
  hMCC = fTfs->make<TH1D>("hMCC","Matthew's Correlation Coefficient; MCC; # Events",201,-1.1,1.1);

  // Make tree to contain bin-by-bin purity and efficiency calculations
  fHitTree = fTfs->make<TTree>("mcanahits","mcanahits");
  fHitTree->Branch("run",&run,"run/I");
  fHitTree->Branch("subrun",&subrun,"subrun/I");
  fHitTree->Branch("event",&event,"event/I");
  fHitTree->Branch("c1",&c1,"c1/i");
  fHitTree->Branch("c2",&c2,"c2/i");
  fHitTree->Branch("trignum",&trignum,"trignum/i");
  fHitTree->Branch("RecoQ",&RecoQ,"RecoQ/D");
  fHitTree->Branch("MCQ",&MCQ,"MCQ/D");
  fHitTree->Branch("RecoT",&RecoT,"RecoT/D");
  fHitTree->Branch("MCT",&MCT,"MCT/D");
  fHitTree->Branch("foundBoth",&foundBoth,"foundBoth/O");

  // Make tree to contain event-by-event purity and efficiency calculations
  fEventTree = fTfs->make<TTree>("mcanaevents","mcanaevents");
  fEventTree->Branch("run",&run,"run/I");
  fEventTree->Branch("subrun",&subrun,"subrun/I");
  fEventTree->Branch("event",&event,"event/I");
  fEventTree->Branch("c1",&c1,"c1/i");
  fEventTree->Branch("c2",&c2,"c2/i");
  fEventTree->Branch("trignum",&trignum,"trignum/i");
  fEventTree->Branch("purity",&eventpurity,"purity/D");
  fEventTree->Branch("efficiency",&eventefficiency,"efficiency/D");
  fEventTree->Branch("chargepurity",&chargepurity,"chargepurity/D");
  fEventTree->Branch("chargeefficiency",&chargeefficiency,"chargeefficiency/D");
  fEventTree->Branch("chargeratio",&eventchargeratio,"chargeratio/D");
  fEventTree->Branch("RecoQ",&eventRecoQ,"RecoQ/D");
  fEventTree->Branch("MCQ",&eventMCQ,"MCQ/D");
  fEventTree->Branch("dQ",&eventdQ,"dQ/D");
  fEventTree->Branch("dT",&eventdT,"dT/D");
  fEventTree->Branch("MCC",&eventMCC,"MCC/D");
  fEventTree->Branch("mcscale",&mcscale,"mcscale/D");
  fEventTree->Branch("tp",&tp,"tp/I");
  fEventTree->Branch("fp",&fp,"fp/I");
  fEventTree->Branch("fn",&fn,"fn/I");
  fEventTree->Branch("tn",&tn,"tn/I");
  fEventTree->Branch("tpc",&tpc,"tpc/D");
  fEventTree->Branch("fpc",&fpc,"fpc/D");
  fEventTree->Branch("fnc",&fnc,"fnc/D");
}

void hit::RobustMCAna::analyze(art::Event const & e)
{
  // Get Services
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<geo::Geometry> fGeom;

  // Get recob::Wire objects to check for dead channels
  art::Handle< std::vector<recob::Wire> > wireHandle;
  if (!e.getByLabel("caldata",wireHandle)) return;

  // Get reconstructed hits
  art::Handle< std::vector<recob::Hit> > hitHandle;
  if (!e.getByLabel(fHitsModuleLabel, hitHandle)) return;

  // Get all T0 products
  art::Handle< std::vector< anab::T0> > t0Handle;
  if (!e.getByLabel(fCounterT0ModuleLabel,t0Handle)) return;

  // get associated external triggers
  art::FindManyP<raw::ExternalTrigger> triggers(t0Handle,e,fCounterT0ModuleLabel);

  // Loop over t0's in this event
  for (size_t i_t0 = 0; i_t0 < t0Handle->size(); i_t0++)
    {
      eventRecoQ = 0;
      eventMCQ = 0;

      // Get this anab::T0 object
      art::Ptr<anab::T0> pt0(t0Handle,i_t0);

      // Find ExternalTriggers associated with this t0 object
      std::vector<art::Ptr<raw::ExternalTrigger> > trigvec = triggers.at(i_t0);

      // Collect all TrigIDs in a vector for processing
      std::vector<unsigned int> evtTriggers;
      for (auto const &trig : trigvec) evtTriggers.push_back(trig->GetTrigID());

      // Fill c1, c2, and trignum variables. Skip if invalid trigger, e.g. incomplete 112 trigger.
      if (!ValidTrigger(evtTriggers,c1,c2,trignum)) continue;

      // Make collection for building MCHit objects (note: not the art MCHit objects)
      std::vector<MCHit> MCHitVec;

      if (fVerbose) std::cout << "Found " << hitHandle->size() << " recob::Hits in total." << std::endl;

      // Start searching for true negatives by looking at channels without signal.
      // Here, start by collecting all kCollection channel numbers
      std::vector<raw::ChannelID_t> absentChannels;
      for (raw::ChannelID_t chan = 0; chan < fGeom->Nchannels(); ++chan)
        {
          if (fCSP->IsBad(chan)) continue;
          if (fGeom->SignalType(chan) != geo::kCollection) continue;
          absentChannels.push_back(chan);
        }

      // Remove channels from absentChannels list if it has a recob::Hit
      for (size_t i_hit = 0; i_hit < hitHandle->size(); ++i_hit)
        {
          art::Ptr<recob::Hit> hit(hitHandle,i_hit);
          std::vector<raw::ChannelID_t>::iterator cit = std::find(absentChannels.begin(),absentChannels.end(),hit->Channel());
          if (cit != absentChannels.end()) absentChannels.erase(cit);
        }

      for (auto sc : bt_serv->SimChannels())
        {
          bool channelexists = false;
          for (size_t i_wire = 0; i_wire < wireHandle->size(); ++i_wire)
            {
              art::Ptr<recob::Wire> pwire(wireHandle,i_wire);
              if (pwire->Channel() == sc->Channel())
                {
                  channelexists = true;
                  break;
                }
            }
          if (!channelexists) continue;

          // also remove channels from absentChannels list if it has a Sim Hit
          if (fCSP->IsBad(sc->Channel())) continue;
          if (fGeom->SignalType(sc->Channel()) != geo::kCollection) continue;
          std::vector<raw::ChannelID_t>::iterator cit = std::find(absentChannels.begin(),absentChannels.end(),sc->Channel());
          if (cit != absentChannels.end()) absentChannels.erase(cit);

          // also remove channels from absent channels list in the opposite drift volume.
          //I do this since I can assume there's one track and it can't be on both sides of the APA at the same Z coordinate.
          //Remember, this is for EW triggers only.
          std::vector<raw::ChannelID_t>::iterator citopp = std::find(absentChannels.begin(),absentChannels.end(),OppositeTPCChannel(sc->Channel()));
          if (citopp != absentChannels.end()) absentChannels.erase(citopp);

          // Initialize a MCHit object for this channel
          MCHit mchit;
          auto const& tdcidemap = sc->TDCIDEMap();
          // Loop over the track id and energies deposited on this channel
          for (auto const& tdcIt : tdcidemap)
            {
              // Get the vector of IDEs
              auto const& ideVec = tdcIt.second;

              // Get the TDC (time) of the deposit
              unsigned short tdc = tdcIt.first;

              // Loop over IDEs
              for (auto const& ideIt : ideVec)
                {
                  // We want only primary muons which cross the detector.
                  // Don't include deltas and other secondary particles because we want to
                  // see how well the reconstruction discounts these other tracks
                  if (abs(ideIt.trackID) != 1) continue;
                  const simb::MCParticle * part = pi_serv->TrackIdToParticle_P(abs(ideIt.trackID));
                  if (part->Mother() != 0) continue;

                  // If the track is a primary muon, then collect this IDE for this channel
                  mchit.idevec.push_back(ideIt);

                  // Sanity check: make sure channel numbers match
                  if (mchit.idevec.size() > 1 && mchit.channel != sc->Channel())
                    throw cet::exception("RobustMCAna") << "Channel IDs don't match!";

                  // Check : if the IDE vector contains different particle codes, the something is weird
                  mchit.channel = sc->Channel();
                  if (mchit.idevec.size() > 1 && mchit.pdg != part->PdgCode())
                    throw cet::exception("RobustMCAna") << "Pdg from previous IDE (" << mchit.pdg << ") doesn't equal this IDE Pdg (" << part->PdgCode() << ")";

                  // Fill the MCHit object, appending energy deposits in the case where more than one IDE exists on the channel
                  // Usually, it's either one or two IDEs per channel with signal. Very rarely more than 2
                  mchit.pdg = part->PdgCode();
                  mchit.particlevec.push_back(part);
                  mchit.tdcvec.push_back(tdc);
                  // Here, adjust the signal size for the MCScale value used in the DataOverlay step
                  // Otherwise, the charges wouldn't make any sense
                  mchit.chargevec.push_back(ideIt.numElectrons*mcscale);
                }
            }
          if (mchit.idevec.size()==0) continue;
          MCHitVec.push_back(mchit);
        }

      // If there are no MCHits and no recob::Hits, then assume an incomplete event,
      // and don't include it in the tracking efficiency denominator
      if (MCHitVec.size() == 0 && hitHandle->size() == 0) continue;

      run = e.run();
      subrun = e.subRun();
      event = e.event();

      // Start dealing with recob::Hits now
      // recoHitMap.first -> a unique ID for this hit
      std::map<size_t,art::Ptr<recob::Hit> > recoHitMap;
      // Collect all recob::Hits in my own format
      for (size_t i_hit = 0; i_hit < hitHandle->size(); ++i_hit)
        {
          art::Ptr<recob::Hit> hit(hitHandle,i_hit);
          recoHitMap[i_hit] = hit;
        }

      // Initialize collections to store various things
      std::vector<double> chargeratios;
      std::vector<double> dTs;
      std::vector<double> dQs;
      std::vector<art::Ptr<recob::Hit> > recoHits;
      std::vector<art::Ptr<recob::Hit> > otherwiseHits;


      // Loop over all MCHits
      for (auto mchit : MCHitVec)
        {
          // Collect information about MCHit
          raw::ChannelID_t channel = mchit.channel;
          int numIDEs = mchit.idevec.size();
          double meanSimTime = TMath::Mean(mchit.tdcvec.begin(),mchit.tdcvec.end(),mchit.chargevec.begin()); // charge weighted mean
          double sigmaSimTime = TMath::RMS(mchit.tdcvec.begin(),mchit.tdcvec.end(),mchit.chargevec.begin()); // charge weighted RMS
          double totalSimCharge = std::accumulate(mchit.chargevec.begin(),mchit.chargevec.end(),0);
          NumElectronsTofC(totalSimCharge); // convert # electrons to fC
          MCQ = totalSimCharge;
          fCToSumADC(MCQ);
          foundBoth = false;

          // Print information about sim hits
          if (fVerbose)
            {
              std::cout << "  Channel: " << channel << ", MeanIDETime: " << meanSimTime << ", RMSTime: " << sigmaSimTime << ", TotalCharge: " << totalSimCharge << ", NumIDEs: " << numIDEs << std::endl;
              std::cout << "::::::::::IDEs:::::::" << std::endl;
              for (size_t i = 0; i < mchit.idevec.size(); ++i)
                {
                  std::cout << "  TrackID: " << mchit.idevec[i].trackID << ", numElectrons: " << mchit.chargevec[i] << ", TDC: " << mchit.tdcvec[i] << ", Channel: " << mchit.channel << ", PDG: " << mchit.pdg << std::endl;
                  //std::cout << *mchit.particlevec[i] << std::endl;
                }
            }

          // Print information about recob::Hits
          if (fVerbose) std::cout << ":::::::::recob::Hits:::::::::" << std::endl;
          for (auto hit = recoHitMap.begin(); hit != recoHitMap.end(); )
            {
              // Bump iterator if recob::Hit is not on the same channel as MCHit
              if (hit->second->Channel() != channel)
                {
                  ++hit;
                  continue;
                }
              if (fVerbose) std::cout << "  Channel: " << hit->second->Channel() << ", PeakTime: " << hit->second->PeakTime() << ", RMS: " << hit->second->RMS() << ", Integral: " << hit->second->Integral() << std::endl;

              // If the sim hit falls within +/- 1 RMS of the recob::Hit, then consider it a good find
              if (hit->second->PeakTimePlusRMS() > meanSimTime && hit->second->PeakTimeMinusRMS() < meanSimTime)
                {
                  // Collect information about the recob::Hit, and the sim Hit
                  foundBoth = true;
                  recoHits.push_back(hit->second);
                  double totalRecoCharge = hit->second->Integral(); // sum ADC
                  RecoQ = totalRecoCharge; // ADC*tick
                  SumADCTofC(totalRecoCharge); // fC
                  double chargeratio = totalRecoCharge / totalSimCharge;
                  RecoT = hit->second->PeakTime();
                  MCT = meanSimTime;
                  double dT = RecoT-MCT;
                  chargeratios.push_back(chargeratio);
                  dTs.push_back(dT);
                  dQs.push_back(RecoQ-MCQ);
                  eventRecoQ += RecoQ;
                  eventMCQ += MCQ;

                  // Fill histograms to compare the reco and MC hits
                  hHitChargeRatio->Fill(chargeratio);
                  hChargeComp->Fill(totalSimCharge,totalRecoCharge); // units of fC
                  hdT->Fill(dT);
                  hdQ->Fill(RecoQ-MCQ);

                  // Remove this reco hit from the map so we can count the recob::Hits left over after associating all sim Hits
                  hit = recoHitMap.erase(hit);
                }
              else
                {
                  otherwiseHits.push_back(hit->second);
                  ++hit;
                }
            }
          fHitTree->Fill();

          if (fVerbose) std::cout << "::::::::::::::::::::::::::::" << std::endl;
        }

      // Calculate stuff and fill TTrees

      tp = recoHits.size();
      tn = absentChannels.size();
      fp = recoHitMap.size();
      fn = MCHitVec.size()-recoHits.size();

      eventpurity      = Purity(     (double)tp, (double)fp );
      eventefficiency  = Efficiency( (double)tp, (double)fn );
      eventMCC         = MCC(        (double)tp, (double)fp, (double)fn, (double)tn );
      eventdT          = VMean(dTs);
      eventchargeratio = VMean(chargeratios);
      eventdQ          = VMean(dQs);

      //if (fVerbose)
      {
        std::cout << "TP: " << tp << "  FP: " << fp << "  TN: " << tn << "  FN: " << fn << std::endl;
        std::cout << "recoHits.size()=" << tp << " MCHitVec.size()=" << MCHitVec.size() << " otherwiseHits.size()=" << otherwiseHits.size() << std::endl;
        std::cout << "Purity: " << eventpurity << "  Efficiency: " << eventefficiency << std::endl;
      }

      // Calculate charge purity and efficiency
      tpc = 0;
      for (auto h : recoHits) tpc += h->Integral();
      SumADCTofC(tpc);

      fpc = 0;
      for (auto h : recoHitMap) fpc += h.second->Integral();
      SumADCTofC(fpc);

      fnc = 0;
      for (auto h : MCHitVec) fnc += std::accumulate(h.chargevec.begin(),h.chargevec.end(),0);
      NumElectronsTofC(fnc); // fC
      fnc -= tpc;

      chargepurity               = Purity( (double)tpc, (double)fpc );
      chargeefficiency           = Efficiency( (double)tpc, (double)fnc );

      //if (fVerbose)
      {
        std::cout << "TPC: " << tpc << "  FPC: " << fpc << "  FNC: " << fnc << std::endl;
        std::cout << "Chg Purity: " << chargepurity << " Chg Efficiency: " << chargeefficiency << std::endl;
        std::cout << "Event RecoQ: " << eventRecoQ << " MCQ: " << eventMCQ << std::endl;
      }

      fEventTree->Fill();

      // Fill histograms
      if (eventpurity>=0) hEventPurity->Fill(eventpurity);
      if (eventefficiency>=0) hEventEfficiency->Fill(eventefficiency);
      if (chargepurity>=0) hEventChargePurity->Fill(chargepurity);
      if (chargeefficiency>=0) hEventChargeEfficiency->Fill(chargeefficiency);
      if (eventMCC>=-1) hMCC->Fill(eventMCC);

    }
}

// Calculate purity based on True Positives and False Positives
double hit::RobustMCAna::Purity(double tp, double fp)
{
  return (tp+fp<=0) ? -1.0 : tp/(tp+fp);
}

// Calculate efficiency based on True Positives and False Negatives
double hit::RobustMCAna::Efficiency(double tp, double fn)
{
  return (tp+fn<=0) ? -1.0 : tp/(tp+fn);
}

// Calculate Matthew's Correlation Coefficient
// Beware: need True Negatives, but not well defined for this case...
double hit::RobustMCAna::MCC(double tp, double fp, double fn, double tn)
{
  return ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)<=0) ? -2.0 : ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
}

// Convert Number of electrons to fC
void hit::RobustMCAna::NumElectronsTofC(double & input)
{
  input *= 1.602e-4;
}

// Convert Sum of ADCs (e.g. from recob::Hit integral) to fC
// Use information from electronics shaping and amplification factors
void hit::RobustMCAna::SumADCTofC(double & input)
{
  input *= (1/ADCtomV)*(1/(preampGain*AmpToAreaFactor));
}

// Convert fC to sum of ADCs (e.g. from simulated charge to units used in reco)
void hit::RobustMCAna::fCToSumADC(double & input)
{
  input *= ADCtomV*preampGain*AmpToAreaFactor;
}

// Robustly calculate the mean of a vector of doubles
// If empty, return -999
double hit::RobustMCAna::VMean(std::vector<double> v)
{
  return (v.size()==0) ? -999 : TMath::Mean(v.size(),v.data());
}

// Given a channel number, find the channel number of the wire which is at the
// same Z coordinate, but on the opposite side of the APA, in the other drift volume
raw::ChannelID_t hit::RobustMCAna::OppositeTPCChannel(raw::ChannelID_t chan)
{
  art::ServiceHandle<geo::Geometry> fGeom;
  std::vector<geo::WireID> attachedWires = fGeom->ChannelToWire(chan);
  for (auto const &w : attachedWires)
    {
      geo::TPCID::TPCID_t tpc = w.TPC;

      if (tpc % 2 == 0) tpc++;
      else tpc--;

      geo::WireID wid(w.Cryostat,tpc,w.Plane,w.Wire);

      raw::ChannelID_t newchan = fGeom->PlaneWireToChannel(wid);
      return newchan;
    }
  return chan;
}

// Test if this is a valid trigger. Remember: this method is a copy of the method in
// RobustHitFinder. If changed in either place, remember to change it in both places.
bool hit::RobustMCAna::ValidTrigger(std::vector<unsigned int> evtTriggers, unsigned int & c1arg, unsigned int & c2arg, unsigned int & trignumarg)
{
  c1arg=999; c2arg=999;
  int contains_111 = 0, contains_112 = 0, contains_113 = 0;
  int contains_Ntrigs = 0, contains_NU = 0, contains_NL = 0, contains_SU = 0, contains_SL = 0;
  int contains_EL = 0, contains_WU = 0, contains_TEL = 0;
  for (size_t i_c = 0; i_c < evtTriggers.size(); i_c++)
    {
      unsigned int trigID = evtTriggers[i_c];
      // for c2: trigID is an unsigned int and always >= 0
      //if (trigID >= 0 && trigID <= 5) contains_SL++;
      if (trigID <= 5) contains_SL++;
      if (trigID >= 6 && trigID <= 15) contains_EL++;
      if (trigID >= 16 && trigID <= 21) contains_NL++;
      if (trigID >= 22 && trigID <= 27) contains_NU++;
      if (trigID >= 28 && trigID <= 37) contains_WU++;
      if (trigID >= 38 && trigID <= 43) contains_SU++;
      if (trigID >= 44 && trigID <= 92) contains_TEL++;
      if (trigID == 111) contains_111++;
      if (trigID == 112) contains_112++;
      if (trigID == 113) contains_113++;
      contains_Ntrigs++;
    }
  if (contains_111 + contains_112 + contains_113 != 1) return false;        // too many/few coincidences!
  if (contains_TEL &&
      (contains_NU || contains_NL || contains_SU || contains_SL || contains_EL || contains_WU)) return false;        // track probably doesn't go through detector
  if (contains_Ntrigs != 3) return false;        // too much/little going on!
  if (contains_111 && (contains_NU || contains_NL || contains_SU || contains_SL)) return false;        // 111 should not have NU/NL/SU/SL
  if (contains_112 && (contains_EL || contains_WU || contains_SU || contains_NL)) return false;        // 112 should not have EL/WU/SU/NL
  if (contains_113 && (contains_EL || contains_WU || contains_NU || contains_SL)) return false;        // 113 should not have EL/WU/NU/SL
  if (contains_111 && (!contains_EL || !contains_WU)) return false;        // incomplete trigger
  if (contains_112 && (!contains_NU || !contains_SL)) return false;        // incomplete trigger
  if (contains_113 && (!contains_SU || !contains_NL)) return false;        // incomplete trigger

  std::vector<unsigned int> counterIDs;
  trignumarg = 0;
  for (size_t i_c = 0; i_c < evtTriggers.size(); i_c++)
    {
      unsigned int trigID = evtTriggers[i_c];
      if (trigID >= 44 && trigID <= 100) continue;
      if (trigID >= 111 && trigID <= 113)
        {
          trignumarg = trigID;
          continue;
        }
      counterIDs.push_back(trigID);
    }
  if (counterIDs.size() != 2) return false;
  if (trignumarg == 0) return false;

  if (trignumarg == 112 || trignumarg == 113)
    {
      if (fCounterPositionMap[counterIDs[0]].first.X() > fCounterPositionMap[counterIDs[1]].first.X())
        {
          c1arg = counterIDs[0];
          c2arg = counterIDs[1];
        }
      else
        {
          c1arg = counterIDs[1];
          c2arg = counterIDs[0];
        }
    }
  else if (trignumarg == 111)
    {
      if (fCounterPositionMap[counterIDs[0]].first.Z() > fCounterPositionMap[counterIDs[1]].first.Z())
        {
          c1arg = counterIDs[0];
          c2arg = counterIDs[1];
        }
      else
        {
          c1arg = counterIDs[1];
          c2arg = counterIDs[0];
        }
    }
  if (c1arg == c2arg) return false;
  if (c1arg == 999 || c2arg == 999) return false;

  return true;
}

DEFINE_ART_MODULE(hit::RobustMCAna)
