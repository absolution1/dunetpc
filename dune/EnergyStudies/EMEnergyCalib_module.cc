////////////////////////////////////////////////////////////////////////
// Class:       EMEnergyCalib
// Module Type: analyzer
// File:        EMEnergyCalib_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), August 2015
//
// Analyser module to produce useful information for characterising
// em showers.
////////////////////////////////////////////////////////////////////////

// Framework includes:
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDAnalyzer.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "MCCheater/BackTracker.h"
#include "Utilities/AssociationUtil.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/MCParticle.h"
#include "Simulation/ParticleList.h"
#include "Simulation/sim.h"

// ROOT & C++ includes
#include <string>
#include <vector>
#include <map>
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TVector3.h"

const int kMaxHits = 10000;

namespace emshower {
  class EMEnergyCalib;
}

class emshower::EMEnergyCalib : public art::EDAnalyzer {
public:

  EMEnergyCalib(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& evt);
  void reconfigure(fhicl::ParameterSet const& p);
  void reset();
  int FindTrackID(art::Ptr<recob::Hit> const& hit);

private:

  // Variables for the tree
  TTree* fTree;
  double trueEnergy;
  double depositU;
  double depositV;
  double depositZ;
  double vertexDetectorDist;
  int    nhits;
  int    hit_tpc        [kMaxHits];
  int    hit_plane      [kMaxHits];
  int    hit_wire       [kMaxHits];
  int    hit_channel    [kMaxHits];
  double hit_peakT      [kMaxHits];
  double hit_charge     [kMaxHits];
  int    hit_truetrackid[kMaxHits];
  int    hit_clusterid  [kMaxHits];

  std::string fHitsModuleLabel;
  std::string fClusterModuleLabel;
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::BackTracker> backtracker;
  art::ServiceHandle<geo::Geometry> geom;

};

emshower::EMEnergyCalib::EMEnergyCalib(fhicl::ParameterSet const& pset) : EDAnalyzer(pset) {
  this->reconfigure(pset);
  fTree = tfs->make<TTree>("EMEnergyCalib","EMEnergyCalib");
  fTree->Branch("TrueEnergy",        &trueEnergy);
  fTree->Branch("DepositU",          &depositU);
  fTree->Branch("DepositV",          &depositV);
  fTree->Branch("DepositZ",          &depositZ);
  fTree->Branch("VertexDetectorDist",&vertexDetectorDist);
  fTree->Branch("NHits",             &nhits);
  fTree->Branch("Hit_TPC",           hit_tpc,        "hit_tpc[NHits]/I");
  fTree->Branch("Hit_Plane",         hit_plane,      "hit_plane[NHits]/I");
  fTree->Branch("Hit_Wire",          hit_wire,       "hit_wire[NHits]/I");
  fTree->Branch("Hit_Channel",       hit_channel,    "hit_channel[NHits]/I");
  fTree->Branch("Hit_PeakT",         hit_peakT,      "hit_peakT[NHits]/D");
  fTree->Branch("Hit_Charge",        hit_charge,     "hit_charge[NHits]/D");
  fTree->Branch("Hit_TrueTrackID",   hit_truetrackid,"hit_truetrackid[NHits]/I");
  fTree->Branch("Hit_ClusterID",     hit_clusterid,  "hit_clusterid[NHits]/I");
}

void emshower::EMEnergyCalib::reconfigure(fhicl::ParameterSet const& pset) {
  fHitsModuleLabel = pset.get<std::string>("HitsModuleLabel");
  fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
}

void emshower::EMEnergyCalib::analyze(art::Event const& evt) {

  /// Analyse function to save information for calibrating shower energies
  /// This is written assuming single particle per event

  this->reset();

  // Get the hits out of the event record
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel(fHitsModuleLabel,hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // Get the clusters out of the event record
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  std::vector<art::Ptr<recob::Cluster> > clusters;
  if (evt.getByLabel(fClusterModuleLabel,clusterHandle))
    art::fill_ptr_vector(clusters, clusterHandle);

  art::FindManyP<recob::Cluster> fmc(hitHandle, evt, fClusterModuleLabel);

  // Look at the hits
  for (unsigned int hitIt = 0; hitIt < hits.size(); ++hitIt) {

    if (hitIt >= kMaxHits) continue;

    // Get the hit
    art::Ptr<recob::Hit> hit = hits.at(hitIt);

    // Fill hit level info
    hit_tpc    [hitIt] = hit->WireID().TPC;
    hit_plane  [hitIt] = hit->WireID().Plane;
    hit_wire   [hitIt] = hit->WireID().Wire;
    hit_peakT  [hitIt] = hit->PeakTime();
    hit_charge [hitIt] = hit->Integral();
    hit_channel[hitIt] = hit->Channel();

    // Find the true track this hit is associated with
    hit_truetrackid[hitIt] = this->FindTrackID(hit);

    // Find the cluster index this hit it associated with (-1 if unclustered)
    if (fmc.isValid()) {
      std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(hitIt);
      if (clusters.size() != 0) {
	hit_clusterid[hitIt] = clusters.at(0)->ID();
      }
      else hit_clusterid[hitIt] = -1;
    }

  }

  // Event level information

  nhits = hits.size();

  const sim::ParticleList& trueParticles = backtracker->ParticleList();
  const simb::MCParticle* trueParticle = trueParticles.Primary(0);

  trueEnergy = trueParticle->Momentum().E();

  // Find the energy deposited on each plane in the TPC
  const std::vector<const sim::SimChannel*>& simChannels = backtracker->SimChannels();
  for (std::vector<const sim::SimChannel*>::const_iterator channelIt = simChannels.begin(); channelIt != simChannels.end(); ++channelIt) {
    int plane = geom->View((*channelIt)->Channel());
    const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = (*channelIt)->TDCIDEMap();
    for (std::map<unsigned short, std::vector<sim::IDE> >::const_iterator tdcIt = tdcidemap.begin(); tdcIt != tdcidemap.end(); ++tdcIt) {
      const std::vector<sim::IDE>& idevec = tdcIt->second;
      for (std::vector<sim::IDE>::const_iterator ideIt = idevec.begin(); ideIt != idevec.end(); ++ideIt) {
	switch (plane) {
	case 0:
	  depositU += ideIt->energy;
	  break;
	case 1:
	  depositV += ideIt->energy;
	  break;
	case 2:
	  depositZ += ideIt->energy;
	  break;
	}
      }
    }
  }

  // Find the distance between the particle vertex and the edge of the detector
  TVector3 vertex = TVector3(trueParticle->Vx(),trueParticle->Vy(),trueParticle->Vz());
  TVector3 end = TVector3(trueParticle->EndX(),trueParticle->EndY(),trueParticle->EndZ());
  TVector3 direction = TVector3(trueParticle->Px(),trueParticle->Py(),trueParticle->Pz()).Unit();

  int distanceStep = 1, steps = 0;
  TVector3 pos;
  bool inTPC = true;
  while (inTPC) {
    pos = end + ( (steps*distanceStep) * direction );
    double currentPos[3]; currentPos[0] = pos.X(); currentPos[1] = pos.Y(); currentPos[2] = pos.Z();
    if (!geom->FindTPCAtPosition(currentPos).isValid)
      inTPC = false;
    ++steps;
  }
  vertexDetectorDist = (end - pos).Mag();

  // Put energies in GeV units
  depositU /= 1000;
  depositV /= 1000;
  depositZ /= 1000;

  fTree->Fill();

  return;

}

int emshower::EMEnergyCalib::FindTrackID(art::Ptr<recob::Hit> const& hit) {
  double particleEnergy = 0;
  int likelyTrackID = 0;
  std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackID(hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }
  return likelyTrackID;
}

void emshower::EMEnergyCalib::reset() {
  trueEnergy = 0;
  depositU = 0;
  depositV = 0;
  depositZ = 0;
  vertexDetectorDist = 0;
  nhits = 0;
  for (int hit = 0; hit < kMaxHits; ++hit) {
    hit_tpc[hit] = 0;
    hit_plane[hit] = 0;
    hit_wire[hit] = 0;
    hit_channel[hit] = 0;
    hit_peakT[hit] = 0;
    hit_charge[hit] = 0;
    hit_clusterid[hit] = 0;
  }
}

DEFINE_ART_MODULE(emshower::EMEnergyCalib)
