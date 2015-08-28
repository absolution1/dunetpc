////////////////////////////////////////////////////////////////////////
// Class:       MultiEMShowers
// Module Type: analyzer
// File:        MultiEMShowers_module.cc
// Author: dorota.stefan@cern.ch robert.sulej@cern.ch
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"

#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Shower.h"
#include "MCCheater/BackTracker.h"

#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include <memory>

#include "DirOfGamma/DirOfGamma.h"

// ROOT includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMathBase.h"

namespace ems {
	class MCinfo;
  class MultiEMShowers;
}

class ems::MCinfo
{
	public:
	MCinfo(const art::Event& evt);
	void Info(const art::Event& evt);

	int GetNgammas(void) const { return fNgammas; }

	double GetMompi0(void) const { return fMompi0; }
	double GetMomGamma1(void) const { return fGammamom1; }
	double GetMomGamma2(void) const { return fGammamom2; }

	double GetCosine(void) { return fCosine; }

	TVector3 GetPospi0(void) const & { return fPi0pos; }
	TVector3 GetPosgamma1(void) const & { return fConvgamma1; }
	TVector3 GetPosgamma2(void) const & { return fConvgamma2; }

	TVector3 GetDirgamma1(void) const & { return fDirgamma1; }
	TVector3 GetDirgamma2(void) const & { return fDirgamma2; }

	private:
	int fNgammas;

	double fMompi0;
	double fGammamom1;
	double fGammamom2;
	
	double fCosine;

	TVector3 fPi0pos;
	TVector3 fConvgamma1;
	TVector3 fConvgamma2;
	TVector3 fDirgamma1;
	TVector3 fDirgamma2;
};

ems::MCinfo::MCinfo(const art::Event& evt)
{
	Info(evt);
}

void ems::MCinfo::Info(const art::Event& evt)
{
	fMompi0 = 0.0; fPi0pos.SetXYZ(0,0,0); 
	fNgammas = 0;
	fCosine = 0.0;

	fGammamom1 = 0.0; fGammamom2 = 0.0;
	fConvgamma1.SetXYZ(0,0,0); fConvgamma2.SetXYZ(0,0,0); 
	fDirgamma1.SetXYZ(0,0,0); fDirgamma2.SetXYZ(0,0,0);

	art::ServiceHandle< cheat::BackTracker > bt;
	const sim::ParticleList& plist = bt->ParticleList();
	for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar)
	{
		simb::MCParticle* particle = ipar->second;

		if ((particle->Process() == "primary") && (particle->PdgCode() == 111))
		{
			fMompi0 = particle->P();
			
			TLorentzVector posvec3 = particle->Position();
			TVector3 pospi0(posvec3.X(), posvec3.Y(), posvec3.Z());
			fPi0pos =  pospi0;

			if (particle->NumberDaughters() == 2)
			{
				fNgammas = particle->NumberDaughters();
				TLorentzVector mom1 = bt->TrackIDToParticle(particle->Daughter(0))->Momentum();
				TLorentzVector mom2 = bt->TrackIDToParticle(particle->Daughter(1))->Momentum();

				TVector3 mom1vec3(mom1.Px(), mom1.Py(), mom1.Pz());
				fGammamom1 = bt->TrackIDToParticle(particle->Daughter(0))->P();
				TVector3 mom2vec3(mom2.Px(), mom2.Py(), mom2.Pz());
				fGammamom2 = bt->TrackIDToParticle(particle->Daughter(1))->P();

				TLorentzVector pos1 = bt->TrackIDToParticle(particle->Daughter(0))->EndPosition();
				TLorentzVector pos2 = bt->TrackIDToParticle(particle->Daughter(1))->EndPosition();
				
				TVector3 pos1vec3(pos1.X(), pos1.Y(), pos1.Z());
				fConvgamma1 = pos1vec3;
				TVector3 pos2vec3(pos2.X(), pos2.Y(), pos2.Z());
				fConvgamma2 = pos2vec3;

				TVector3 vecnorm1 = mom1vec3 * (1.0 / mom1vec3.Mag());
				fDirgamma1 = vecnorm1;
				TVector3 vecnorm2 = mom2vec3 * (1.0 / mom2vec3.Mag());
				fDirgamma2 = vecnorm2;
		
				fCosine = fDirgamma1 * fDirgamma2;
				break;
			}
			else
			{
				fNgammas = particle->NumberDaughters();
			}
		}	
	}
}

class ems::MultiEMShowers : public art::EDAnalyzer {
public:
  explicit MultiEMShowers(fhicl::ParameterSet const & p);

  MultiEMShowers(MultiEMShowers const &) = delete;
  MultiEMShowers(MultiEMShowers &&) = delete;
  MultiEMShowers & operator = (MultiEMShowers const &) = delete;
  MultiEMShowers & operator = (MultiEMShowers &&) = delete;

	void beginJob() override;
  
  void analyze(art::Event const & e) override;

	void reconfigure(fhicl::ParameterSet const& p);

private:
	int fConvGood;
	int fConvWrong;
	int fConvBothGood;

	// ROOT
	TTree* fEvTree; 
	int fEvNumber;
	int fNGroups;
	int fNShs;

	// mc 
	double fPi0mom;
	int fNgammas;

	TTree* fShTree; 	
	double fStartX; double fStartY; double fStartZ;
	double fMCrecovtx; double fMCrecoTh;
	double fDistConvrecomc1; double fDistConvrecomc2;

  std::string fHitsModuleLabel;
	std::string fCluModuleLabel;
	std::string fTrk3DModuleLabel;
	std::string fVtxModuleLabel;
	std::string fShsModuleLabel;
};


ems::MultiEMShowers::MultiEMShowers(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  
{
	fConvGood = 0;
	fConvWrong = 0;
	fConvBothGood = 0;
	reconfigure(p);
}

void ems::MultiEMShowers::reconfigure(fhicl::ParameterSet const& p)
{
  fHitsModuleLabel = p.get< std::string >("HitsModuleLabel");
  fCluModuleLabel = p.get< std::string >("ClustersModuleLabel");
  fTrk3DModuleLabel = p.get< std::string >("Trk3DModuleLabel");
  fVtxModuleLabel = p.get< std::string >("VtxModuleLabel");
	fShsModuleLabel = p.get< std::string >("ShsModuleLabel");
  return;
}

void ems::MultiEMShowers::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;

	fEvTree = tfs->make<TTree>("MultiShowers", "showers3d");
	fEvTree->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fEvTree->Branch("fNGroups", &fNGroups, "fNGroups/I");
	fEvTree->Branch("fPi0mom", &fPi0mom, "fPi0mom/D");
	fEvTree->Branch("fNgammas", &fNgammas, "fNgammas/I");
	fEvTree->Branch("fMCrecovtx", &fMCrecovtx, "fMCrecovtx/D");
	fEvTree->Branch("fMCrecoTh", &fMCrecoTh, "fMCrecoTh/D");
	fEvTree->Branch("fConvGood", &fConvGood, "fConvGood/I");
	fEvTree->Branch("fConvWrong", &fConvWrong, "fConvWrong/I");
	fEvTree->Branch("fConvBothGood", &fConvBothGood, "fConvBothGood/I");

	fShTree = tfs->make<TTree>("Shower", "conv point");
	fShTree->Branch("fNShs", &fNShs, "fNShs/I");
	fShTree->Branch("fStartX", &fStartX, "fStartX/D");
	fShTree->Branch("fStartY", &fStartY, "fStartY/D");
	fShTree->Branch("fStartZ", &fStartZ, "fStartZ/D");
}

void ems::MultiEMShowers::analyze(art::Event const & e)
{
	fEvNumber = e.id().event();
	fNGroups = 0; fNShs = 0;
	fStartX = 0.0; fStartY = 0.0; fStartZ = 0.0;
	fPi0mom = 0.0; fNgammas = 0;
	fDistConvrecomc1 = 0.0; fDistConvrecomc2 = 0.0;
	fMCrecovtx = -10.0;
	fMCrecoTh = -400.0;

	ems::MCinfo mc(e);
	fPi0mom = mc.GetMompi0();
	fNgammas = mc.GetNgammas();
	TVector3 pospi0 = mc.GetPospi0();

	double cosinemc = mc.GetCosine();
	TVector3 convp[2]; 
	convp[0] = mc.GetPosgamma1();
	convp[1] = mc.GetPosgamma2();
	const double maxdist = 2.0; //cm

	art::Handle< std::vector<recob::Shower> > shsListHandle;
  art::Handle< std::vector<recob::Track> > trkListHandle;
	art::Handle< std::vector<recob::Vertex> > vtxListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	art::Handle< std::vector<recob::Hit> > hitListHandle;
	if (e.getByLabel(fShsModuleLabel, shsListHandle) &&
			e.getByLabel(fTrk3DModuleLabel, trkListHandle) &&
	    e.getByLabel(fVtxModuleLabel, vtxListHandle) &&
	    e.getByLabel(fCluModuleLabel, cluListHandle) &&
	    e.getByLabel(fHitsModuleLabel, hitListHandle))
	{
		art::FindManyP< recob::Cluster > cluFromShs(shsListHandle, e, fShsModuleLabel);
		art::FindManyP< recob::Cluster > cluFromTrk(trkListHandle, e, fTrk3DModuleLabel);
		art::FindManyP< recob::Vertex > vtxFromTrk(trkListHandle, e, fVtxModuleLabel);
		art::FindManyP< recob::Hit > hitFromClu(cluListHandle, e, fCluModuleLabel);

		fNGroups = shsListHandle->size();

		int countph = 0;
		if (fNgammas == 2)
		{
			int idph = -1; 
			for (size_t s = 0; s < shsListHandle->size(); ++s)
			{
				const recob::Shower& sh = (*shsListHandle)[s];
				double mindist = maxdist; bool found = false; 
			
				for (int i = 0; i < fNgammas; ++i)
				{
					double dist = sqrt(pma::Dist2(sh.ShowerStart(), convp[i]));
					if ((dist < mindist) && (idph != i)) 
					{ mindist =  dist; idph = i; found = true; }
				}
				if (found) { fConvGood++; countph++; }
				else { fConvWrong++; }
			}
			if (countph == 2) fConvBothGood++;
		
			if (countph == 2)
			{
				for (size_t s = 0; s < shsListHandle->size(); ++s)	
				{
					const recob::Shower& sh = (*shsListHandle)[s];
					TVector3 pos = sh.ShowerStart(); 
					fStartX = pos.X(); fStartY = pos.Y(); fStartZ = pos.Z();
					fNShs++;
		
					fShTree->Fill();	
				}
			}
		}
		// compute the crossing point
		if (countph == 2)
		{
			std::vector< std::pair<TVector3, TVector3> > lines;
			const recob::Shower& sh1 = (*shsListHandle)[0];
			const recob::Shower& sh2 = (*shsListHandle)[1];

			std::pair<TVector3, TVector3> frontback1(sh1.ShowerStart(), sh1.ShowerStart() + sh1.Direction());
			std::pair<TVector3, TVector3> frontback2(sh2.ShowerStart(), sh2.ShowerStart() + sh2.Direction());
			lines.push_back(frontback1); lines.push_back(frontback2);

			TVector3 result;
			pma::SolveLeastSquares3D(lines, result); // mse.

			double dist1_0 = pma::Dist2(result, sh1.ShowerStart());
			double dist2_0 = pma::Dist2(result, sh1.ShowerStart() + sh1.Direction());
			double dist1_1 = pma::Dist2(result, sh2.ShowerStart());
			double dist2_1 = pma::Dist2(result, sh2.ShowerStart() + sh2.Direction());
			if ((dist1_0 > dist2_0) || (dist1_1 > dist2_1)) {}
			else
			{
				fMCrecovtx = std::sqrt(pma::Dist2(pospi0, result));

				double cosine_reco = sh1.Direction() * sh2.Direction();
				double threco = 180.0F * (std::acos(cosine_reco)) / TMath::Pi();
				double thmc = 180.0F * (std::acos(cosinemc)) / TMath::Pi();

				fMCrecoTh = threco - thmc;
			}
		}

	}

	fEvTree->Fill();		
}

DEFINE_ART_MODULE(ems::MultiEMShowers)

