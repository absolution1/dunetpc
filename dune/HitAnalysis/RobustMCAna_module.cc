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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireGeo.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/sim.h"

#include <algorithm>
#include "TMath.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

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

  struct MCHit {
    raw::ChannelID_t channel;
    int pdg;
    std::vector<sim::IDE> idevec;
    std::vector<unsigned short> tdcvec;
    std::vector<double> chargevec;
    std::vector<const simb::MCParticle *> particlevec;
  };

  std::string fHitsModuleLabel;
  bool fVerbose;
  double fPreviousMCScale;

  const lariov::ChannelStatusProvider* fCSP;

  TH1D* hEventPurity;// aka Precision
  TH1D* hEventEfficiency;// aka Recall
  TH1D* hEventChargePurity;
  TH1D* hEventChargeEfficiency;
  TH1D* hHitChargeRatio;
  TH2D* hChargeComp;
  TH1D* hdT;
  TH1D* hMCC;

  double mcscale;

  TTree * fEventTree;
  double eventpurity;
  double eventefficiency;
  double chargepurity;
  double chargeefficiency;
  double eventchargeratio;
  double eventdT;
  double eventMCC;
  int tp;
  int fp;
  int fn;
  int tn;
  double tpc;
  double fpc;
  double fnc;

  TTree * fBinTree;
  int bin;
  double bincenter;
  double binpurity;
  double binefficiency;
  double binchargepurity;
  double binchargeefficiency;
  int bintp;
  int binfp;
  int binfn;
  double bintpc;
  double binfpc;
  double binfnc;
  double binchargeratio;

  int fNbins;

  detinfo::DetectorProperties const * fDetProp;
  detinfo::DetectorClocks const * fClks;

  double preampGain;
  double AmpToAreaFactor;
  double ADCtomV;
};


hit::RobustMCAna::RobustMCAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
  fClks(lar::providerFrom<detinfo::DetectorClocksService>())
{
  preampGain = p.get<double>("PreampGainSetting",14.0);// mV/fC
  AmpToAreaFactor = p.get<double>("AmpToAreaFactor",7.615);// convert amplitude of shaping curve to area of pseudo-gaussian
  ADCtomV = p.get<double>("ADCtomV",2.808);// ADC/mV, digitization

  fHitsModuleLabel = p.get<std::string>("HitsModuleLabel","robusthit");
  fVerbose = p.get<bool>("Verbose",false);
  fPreviousMCScale = p.get<double>("PreviousMCScale");
  fNbins = p.get<int>("NumberBins",22);

  fCSP = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

  art::ServiceHandle<art::TFileService> fTfs;
  hEventPurity = fTfs->make<TH1D>("hEventPurity","Purity of Reco Hit Selection; TP/(TP+FP); # Events",200,0,1.1);
  hEventEfficiency = fTfs->make<TH1D>("hEventEfficiency","Efficiency of Reco Hit Selection; TP/(TP+FN); # Events",200,0,1.1);
  hEventChargePurity = fTfs->make<TH1D>("hEventChargePurity","Purity of Reco Hit Charge Selection; ChgTP/(ChgTP+ChgFP); # Events",200,0,1.1);
  hEventChargeEfficiency = fTfs->make<TH1D>("hEventChargeEfficiency","Efficiency of Reco Hit Charge Selection; ChgTP/(ChgTP+ChgFN); # Events",200,0,1.1);
  hHitChargeRatio = fTfs->make<TH1D>("hHitChargeRatio","Ratio of reco hit charge to MC hit charge; Reco Charge / MC Charge; # Events",200,0,3);
  hChargeComp = fTfs->make<TH2D>("hChargeComp","Comparison of reco and MC charges; MC Charge (fC); Reco Charge (fC)",400,0,30000,400,0,30000);
  hdT = fTfs->make<TH1D>("hdT","#DeltaT between MC hit and Reco Hit; Reco Hit Time - MC Hit Time (TDC ticks); # Events",200,-30,30);
  hMCC = fTfs->make<TH1D>("hMCC","Matthew's Correlation Coefficient; MCC; # Events",201,-1.1,1.1);

  fBinTree = fTfs->make<TTree>("mcanabins","mcanabins");
  fBinTree->Branch("bin",&bin,"bin/I");
  fBinTree->Branch("bincenter",&bincenter,"bincenter/D");
  fBinTree->Branch("purity",&binpurity,"purity/D");
  fBinTree->Branch("efficiency",&binefficiency,"efficiency/D");
  fBinTree->Branch("chargepurity",&binchargepurity,"chargepurity/D");
  fBinTree->Branch("chargeefficiency",&binchargeefficiency,"chargeefficiency/D");
  fBinTree->Branch("tp",&bintp,"tp/I");
  fBinTree->Branch("fp",&binfp,"fp/I");
  fBinTree->Branch("fn",&binfn,"fn/I");
  fBinTree->Branch("tpc",&bintpc,"tpc/D");
  fBinTree->Branch("fpc",&binfpc,"fpc/D");
  fBinTree->Branch("fnc",&binfnc,"fnc/D");
  fBinTree->Branch("mcscale",&mcscale,"mcscale/D");
  fBinTree->Branch("chargeratio",&binchargeratio,"chargeratio/D");

  fEventTree = fTfs->make<TTree>("mcanaevents","mcanaevents");
  fEventTree->Branch("purity",&eventpurity,"purity/D");
  fEventTree->Branch("efficiency",&eventefficiency,"efficiency/D");
  fEventTree->Branch("chargepurity",&chargepurity,"chargepurity/D");
  fEventTree->Branch("chargeefficiency",&chargeefficiency,"chargeefficiency/D");
  fEventTree->Branch("chargeratio",&eventchargeratio,"chargeratio/D");
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
  art::ServiceHandle<cheat::BackTracker> bt;
  art::ServiceHandle<geo::Geometry> fGeom;

  art::Handle< std::vector<recob::Hit> > hitHandle;
  if (!e.getByLabel(fHitsModuleLabel, hitHandle))
    return;

  mcscale = fPreviousMCScale;

  std::vector<MCHit> MCHitVec;

  if (fVerbose) std::cout << "Found " << hitHandle->size() << " recob::Hits in total." << std::endl;

  std::vector<raw::ChannelID_t> absentChannels;
  for (raw::ChannelID_t chan = 0; chan < fGeom->Nchannels(); ++chan)
    {
      if (fCSP->IsBad(chan)) continue;
      if (fGeom->SignalType(chan) != geo::kCollection) continue;
      absentChannels.push_back(chan);
    }

  for (size_t i_hit = 0; i_hit < hitHandle->size(); ++i_hit)
    {
      art::Ptr<recob::Hit> hit(hitHandle,i_hit);
      std::vector<raw::ChannelID_t>::iterator cit = std::find(absentChannels.begin(),absentChannels.end(),hit->Channel());
      if (cit != absentChannels.end()) absentChannels.erase(cit);
    }

  for (auto sc : bt->SimChannels())
    {
      if (fCSP->IsBad(sc->Channel())) continue;
      if (fGeom->SignalType(sc->Channel()) != geo::kCollection) continue; 
      std::vector<raw::ChannelID_t>::iterator cit = std::find(absentChannels.begin(),absentChannels.end(),sc->Channel());
      if (cit != absentChannels.end()) absentChannels.erase(cit);

      std::vector<raw::ChannelID_t>::iterator citopp = std::find(absentChannels.begin(),absentChannels.end(),OppositeTPCChannel(sc->Channel()));
      if (citopp != absentChannels.end()) absentChannels.erase(citopp);

      MCHit mchit;
      auto const& tdcidemap = sc->TDCIDEMap();
      for (auto const& tdcIt : tdcidemap)
	{
	  auto const& ideVec = tdcIt.second;
	  unsigned short tdc = tdcIt.first;
	  for (auto const& ideIt : ideVec)
	    {
	      if (abs(ideIt.trackID) != 1) continue;
	      const simb::MCParticle * part = bt->TrackIDToParticle(abs(ideIt.trackID)); 
	      if (part->Mother() != 0) continue;
	      mchit.idevec.push_back(ideIt);
	      if (mchit.idevec.size() > 1 && mchit.channel != sc->Channel()) 
		throw cet::exception("RobustMCAna") << "Channel IDs don't match!";
	      mchit.channel = sc->Channel();
	      if (mchit.idevec.size() > 1 && mchit.pdg != part->PdgCode())
		throw cet::exception("RobustMCAna") << "Pdg from previous IDE (" << mchit.pdg << ") doesn't equal this IDE Pdg (" << part->PdgCode() << ")";
	      mchit.pdg = part->PdgCode();
	      mchit.particlevec.push_back(part);
	      mchit.tdcvec.push_back(tdc);
	      mchit.chargevec.push_back(ideIt.numElectrons*fPreviousMCScale);
	    }
	}
      MCHitVec.push_back(mchit);
      
    }

  if (MCHitVec.size() == 0 && hitHandle->size() == 0) return;

  std::map<size_t,art::Ptr<recob::Hit> > recoHitMap;
  for (size_t i_hit = 0; i_hit < hitHandle->size(); ++i_hit)
    {
      art::Ptr<recob::Hit> hit(hitHandle,i_hit);
      recoHitMap[i_hit] = hit;
    }

  std::vector<double> chargeratios;
  std::map<int,std::vector<double> > chargeratiosbins;
  std::vector<double> dTs;
  std::vector<art::Ptr<recob::Hit> > recoHits;
  std::vector<MCHit> MCHits;
  int numMCHits = 0, numRecoHits = 0;
  for (auto mchit : MCHitVec)
    {
      if (mchit.idevec.size() == 0) continue;
      raw::ChannelID_t channel = mchit.channel;
      int numIDEs = mchit.idevec.size();
      double meanSimTime = TMath::Mean(mchit.tdcvec.begin(),mchit.tdcvec.end(),mchit.chargevec.begin());
      double sigmaSimTime = TMath::RMS(mchit.tdcvec.begin(),mchit.tdcvec.end(),mchit.chargevec.begin());
      double totalSimCharge = std::accumulate(mchit.chargevec.begin(),mchit.chargevec.end(),0);
      //totalSimCharge *= fDetProp->ElectronsToADC();
      totalSimCharge *= 1.602e-4; // convert # electrons to fC
      
      if (fVerbose) 
	{
	  std::cout << "  Channel: " << channel << ", MeanIDETime: " << meanSimTime << ", RMSTime: " << sigmaSimTime << ", TotalCharge: " << totalSimCharge << ", NumIDEs: " << numIDEs << std::endl;
	  std::cout << "::::::::::IDEs:::::::" << std::endl;
	  for (size_t i = 0; i < mchit.idevec.size(); ++i)
	    {
	      std::cout << "  TrackID: " << mchit.idevec[i].trackID << ", numElectrons: " << mchit.chargevec[i] << ", TDC: " << mchit.tdcvec[i] << ", Channel: " << mchit.channel << ", PDG: " << mchit.pdg << std::endl;
	      std::cout << *mchit.particlevec[i] << std::endl;
	    }
	}
      MCHits.push_back(mchit);
      ++numMCHits;
      if (fVerbose) std::cout << ":::::::::recob::Hits:::::::::" << std::endl;
      for (auto hit = recoHitMap.begin(); hit != recoHitMap.end(); )
	{
	  if (hit->second->Channel() != channel) 
	    {
	      ++hit;
	      continue;
	    }
	  
	  if (fVerbose) std::cout << "  Channel: " << hit->second->Channel() << ", PeakTime: " << hit->second->PeakTime() << ", RMS: " << hit->second->RMS() << ", Integral: " << hit->second->Integral() << std::endl;
	  if (hit->second->PeakTimePlusRMS() > meanSimTime && hit->second->PeakTimeMinusRMS() < meanSimTime)
	    {
	      recoHits.push_back(hit->second);
	      numRecoHits++;
	      double totalRecoCharge = hit->second->Integral(); // sum ADC
	      totalRecoCharge *= (1/ADCtomV); // sum ADC -> sum mV
	      totalRecoCharge *= (1/(preampGain*AmpToAreaFactor)); // sum mV -> fC
	      hHitChargeRatio->Fill(totalRecoCharge / totalSimCharge);
	      chargeratios.push_back(totalRecoCharge / totalSimCharge);
	      hChargeComp->Fill(totalSimCharge,totalRecoCharge);// units of fC
	      hdT->Fill(hit->second->PeakTime()-meanSimTime);
	      dTs.push_back(hit->second->PeakTime()-meanSimTime);
	      hit = recoHitMap.erase(hit);
	      for (bin = 0; bin < fNbins; ++bin)
		{
		  double min = bin*(2012.0/fNbins);
		  double max = (bin+1)*(2012.0/fNbins);
		  double t = fClks->TPCTick2TrigTime(hit->second->PeakTime());
		  if (t >= min && t < max)
		    {
		      chargeratiosbins[bin].push_back(totalRecoCharge / totalSimCharge);
		    }
		}
	    }
	  else
	    {
	      ++hit;
	    }
	}
      if (fVerbose) std::cout << "::::::::::::::::::::::::::::" << std::endl;
    }

  tp = recoHits.size();
  tn = absentChannels.size();
  fp = recoHitMap.size();
  fn = MCHits.size()-recoHits.size();

  double purity = (double)tp/((double)(tp+fp));
  double efficiency = (double)tp/((double)(tp+fn));
  double mcc = ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
  std::cout << "TP: " << tp << "  FP: " << fp << "  TN: " << tn << "  FN: " << fn << "  MCC: " << mcc << std::endl;

  eventpurity = (tp+fp==0) ? -1 : purity;
  eventefficiency = (tp+fn==0) ? -1 : efficiency;
  eventMCC = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)==0) ? -1 : mcc;
  eventdT = (dTs.size()==0) ? -999 : TMath::Mean(dTs.size(),dTs.data());
  eventchargeratio = (chargeratios.size()==0) ? -999 : TMath::Mean(chargeratios.size(),chargeratios.data());

  tpc = 0;
  for (auto h : recoHits)
    {
      tpc += h->Integral();
    }
  tpc *= (1/ADCtomV);
  tpc *= (1/(preampGain*AmpToAreaFactor));

  fpc = 0;
  for (auto h : recoHitMap)
    {
      fpc += h.second->Integral();
    }
  fpc *= (1/ADCtomV);
  fpc *= (1/(preampGain*AmpToAreaFactor));

  fnc = 0;
  for (auto h : MCHits)
    {
      fnc += std::accumulate(h.chargevec.begin(),h.chargevec.end(),0);
    }
  fnc *= 1.602e-4; // fC
  fnc -= tpc;

  chargepurity = (tpc+fpc==0) ? -1 : tpc/(tpc+fpc);
  chargeefficiency = (tpc+fnc==0) ? -1 : tpc/(tpc+fnc);
  std::cout << "TPC: " << tpc << "  FPC: " << fpc << "  FNC: " << fnc << "  CPur: " << chargepurity << "  CEff: " << chargeefficiency << std::endl;
  
  fEventTree->Fill();

  if (tp+fp != 0) hEventPurity->Fill(purity);
  if (tp+fn != 0) hEventEfficiency->Fill(efficiency);
  if (tpc+fpc != 0) hEventChargePurity->Fill(chargepurity);
  if (tpc+fnc != 0) hEventChargeEfficiency->Fill(chargeefficiency);
  if ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)!=0) hMCC->Fill(mcc);

  for (bin = 0; bin < fNbins; ++bin)
    {
      double min = bin*(2012.0/fNbins);
      double max = (bin+1)*(2012.0/fNbins);
      bincenter = (min+max)/2.0;
      bintp = 0;
      bintpc = 0;
      for (auto hit : recoHits)
	{
	  double t = fClks->TPCTick2TrigTime(hit->PeakTime());
	  if (t >= min && t < max) 
	    {
	      bintp++;
	      bintpc += hit->Integral();
	    }
	}
      bintpc *= (1/ADCtomV);
      bintpc *= (1/(preampGain*AmpToAreaFactor));

      binfp = 0;
      binfpc = 0;
      for (auto hit : recoHitMap)
	{
	  double t = fClks->TPCTick2TrigTime(hit.second->PeakTime());
	  if (t >= min && t < max) 
	    {
	      binfp++;
	      binfpc += hit.second->Integral();
	    }
	}
      binfpc *= (1/ADCtomV);
      binfpc *= (1/(preampGain*AmpToAreaFactor));

      binfn = 0;
      binfnc = 0;
      
      for (auto hit : MCHits)
	{
	  double meanSimTime = TMath::Mean(hit.tdcvec.begin(),hit.tdcvec.end(),hit.chargevec.begin());
	  double t = fClks->TPCTick2TrigTime(meanSimTime);
	  if (t >= min && t < max) 
	    {
	      binfn++;
	      binfnc += std::accumulate(hit.chargevec.begin(),hit.chargevec.end(),0);
	    }
	}
      binfn -= bintp;
      binfnc *= 1.602e-4;
      binfnc -= bintpc;

      binpurity = (bintp+binfp<=0) ? -1.0 : ((double)bintp)/((double)(bintp+binfp));
      binefficiency = (bintp+binfn<=0) ? -1.0 : ((double)bintp)/((double)(bintp+binfn));
      binchargepurity = (bintpc+binfpc<=0) ? -1.0 : bintpc/(bintpc+binfpc);
      binchargeefficiency = (bintpc+binfnc<=0) ? -1.0 : bintpc/(bintpc+binfnc);
      binchargeratio = (chargeratiosbins[bin].size()==0) ? -999 : TMath::Mean(chargeratiosbins[bin].size(),chargeratiosbins[bin].data());

      
      std::cout << "Event " << e.event() << "  bin " << bin << "  fp=" << binfp << " fn=" << binfn << " tp=" << bintp << " pur=" << binpurity << " eff=" << binefficiency << " cpur=" << binchargepurity << " ceff=" << binchargeefficiency << "  cratio=" << binchargeratio << std::endl;

      fBinTree->Fill();
    }
}

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

DEFINE_ART_MODULE(hit::RobustMCAna)
