#ifndef ChargeDist_H
#define ChargeDist_H

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EventSelector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"

// C++ includes
#include <vector>
#include <string>
#include <iostream>

namespace {
    class ChargeDist : public art::EDAnalyzer {
        public:
        
        //Standard constructor for an ART module
        explicit ChargeDist(fhicl::ParameterSet const&);
        virtual ~ChargeDist();
        
        //The analyzer routine, called once per event
        void analyze(const art::Event& event) override;
        
        //Functions
        float rawDigitChargeDist(const art::Event& event, std::string labelToGetRawDigits);
        float rawDigitChargeDistThreshold(const art::Event& event, std::string labelToGetRawDigits, float threshold);
        
        double getPhotonFlashTime(const art::Event& event, std::string OpFlashLabel);
        std::vector<double> getListOfUniqueMCParticles(const art::Event& event, std::string MCParticleLabel, int PDGcode);
        std::vector<double> getPositionDifferences(const art::Event& event, std::string labelToGetTracks, double PhotonFlashTime, double truthXPos);
        
        private:
        std::string prim;
        
        //Parameters read from the fhicl file
        std::string fDetSimProducerLabel;
        std::string fCalDataModuleLabel;
        std::string MCParticleLabel;
        std::string HitSpacePointLabel;
        std::string OpFlashLabel;
        std::string HitClusterLabel;
        std::string HitTrackLabel;
        std::string CaloLabel;
        
        //Tree + branches to be filled in the module
        TTree* tr;
        
        //event information
        double fEvent;
        
        //rawdigit information
        double numRawDigits;
        double numPositiveSignals;
        double numNegativeSignals;
        
        //charge distribution stuff
        double charge;
        double chargeCorrect;
        double chargeRecobCorrect;
        double chargeRecob25PerCentCorrect;
        double chargeRecob50PerCentCorrect;
        double chargeRecob75PerCentCorrect;
        double charge_1MeVThresh;
        double charge_15MeVThresh;
        double charge_2MeVThresh;
                
        //amount of correction
        double DCFactor_Truth;
        double FractionalDC_Truth;
        double DCFactor_Reco;
        double FractionalDC_Reco;
        
        //efficiency counter
        double efficiencyCounter;
        
        //primary truth information
        double nuEnergy;
        double summedTruthEnergyTotal;
        double numNeutrons;
        double summedTruthEnergyNeutron;
        double numProtons;
        double summedTruthEnergyProton;
        double numElectrons;
        double summedTruthEnergyElectron;
        double numGammas;
        double summedTruthEnergyGamma;
        double meanTruthEnergyGamma;
        double numTotalGammas;
        
        //conversion factors to MeV
        float slope = 591.27;
        float Eslope = 1.32336;
        float intercept = -295.984;
        float Eintercept = 24.899;
        //define thresholds calculated from conversion factors above
        float threshold_1MeV = 295.286;
        float threshold_15MeV = 590.921;
        float threshold_2MeV = 886.556;
        
        //position resolution study
        int indexOfSmallestPosDiff;
        double truthXPos;
        double truthYPos;
        double truthZPos;
        //position of trajectory point
        double recoTP_XPos;
        double recoTP_YPos;
        double recoTP_ZPos;
        //X position of first hit
        double recoXPos;
        
        double photonFlashTime;
        double posDiffFirstHit; //first hit
        double posDiffFirstTP; //first trajectory point
        
        double posDiffUltimate; //track direction not flipped, first track is prim electron
        
        //distance from closest photon detector
        double posDiffDistanceOpDet;
        double posDiffDistanceXFirstHit;
        double posDiffDistanceXFirstTP;
        //3d position resolution: first trajectory point
        double posDiffDist3D;
        double posDiffYFirstTP;
        double posDiffZFirstTP;
        
        //hit charge
        double totalHitCharge;
        
        //drift times
        double driftTime;
        double driftRecobTime;
        
        
        
      };
    }

    #endif

    namespace {
        DEFINE_ART_MODULE(ChargeDist)
    }

    namespace {
        //---------------------------------------------------------------------------
        //Constructor
        //---------------------------------------------------------------------------
        ChargeDist::ChargeDist(fhicl::ParameterSet const& parameterSet): EDAnalyzer(parameterSet){
            
            prim = "primary";
            //Read the fcl-file
            fDetSimProducerLabel = parameterSet.get< std::string >("DetSimLabel");
            fCalDataModuleLabel = parameterSet.get< std::string >("CalDataModuleLabel");
            MCParticleLabel = parameterSet.get<std::string>("MCParticleLabel");
            HitSpacePointLabel = parameterSet.get< std::string >("HitSpacePointLabel");
            OpFlashLabel = parameterSet.get< std::string >("OpFlashLabel");
            HitClusterLabel = parameterSet.get<std::string>("HitClusterLabel");
            HitTrackLabel = parameterSet.get<std::string>("HitTrackLabel");
            CaloLabel = parameterSet.get<std::string>("CaloLabel");

            //Get a service that creates ROOT histograms and trees
            art::ServiceHandle< art::TFileService > tfs;
                        
            //initialize the tree, branches
            tr = tfs->make<TTree>("DataStoring", "DataStoring");
            
            //event information
            tr->Branch("fEvent", &fEvent, "fEvent/D");
            
            //raw digit information
            tr->Branch("numRawDigits", &numRawDigits, "numRawDigits/D");
            tr->Branch("numPositiveSignals", &numPositiveSignals, "numPositiveSignals/D");
            tr->Branch("numNegativeSignals", &numNegativeSignals, "numNegativeSignals/D");
            
            //charge 
            tr->Branch("Charge", &charge, "Charge/D");
            tr->Branch("ChargeTruthDriftCorrect", &chargeCorrect, "ChargeTruthDriftCorrect/D");
            tr->Branch("ChargeRecobDriftCorrect", &chargeRecobCorrect, "ChargeRecobDriftCorrect/D");         
            tr->Branch("chargeRecob25PerCentCorrect", &chargeRecob25PerCentCorrect, "chargeRecob25PerCentCorrect/D");
            tr->Branch("chargeRecob50PerCentCorrect", &chargeRecob50PerCentCorrect, "chargeRecob50PerCentCorrect/D");
            tr->Branch("chargeRecob75PerCentCorrect", &chargeRecob75PerCentCorrect, "chargeRecob75PerCentCorrect/D");
			tr->Branch("charge_1MeVThresh", &charge_1MeVThresh, "charge_1MeVThresh/D");
            tr->Branch("charge_15MeVThresh", &charge_15MeVThresh, "charge_15MeVThresh/D");
            tr->Branch("charge_2MeVThresh", &charge_2MeVThresh, "charge_2MeVThresh/D");           
            
            //correction factors
            tr->Branch("DCFactor_Truth", &DCFactor_Truth, "DCFactor_Truth/D");
            tr->Branch("FractionalDC_Truth", &FractionalDC_Truth, "FractionalDC_Truth/D");
            tr->Branch("DCFactor_Reco", &DCFactor_Reco, "DCFactor_Reco/D");
            tr->Branch("FractionalDC_Reco", &FractionalDC_Reco, "FractionalDC_Reco/D");
            
            tr->Branch("efficiencyCounter", &efficiencyCounter, "efficiencyCounter/D");
            
            //drift times
            tr->Branch("PhotonTime", &photonFlashTime, "PhotonTime/D");
            tr->Branch("driftTime", &driftTime, "driftTime/D");
            tr->Branch("driftRecobTime", &driftRecobTime, "driftRecobTime/D");
            
            //position resolution study
            tr->Branch("TruthXPos", &truthXPos, "TruthXPos/D");
            tr->Branch("TruthYPos", &truthYPos, "TruthYPos/D");
            tr->Branch("TruthZPos", &truthZPos, "TruthZPos/D");
            //position of trajectory point
            tr->Branch("recoTP_XPos", &recoTP_XPos, "recoTP_XPos/D");
            tr->Branch("recoTP_YPos", &recoTP_YPos, "recoTP_YPos/D");
            tr->Branch("recoTP_ZPos", &recoTP_ZPos, "recoTP_ZPos/D");
            //x position of first hit
            tr->Branch("recoXPos", &recoXPos, "recoXPos/D");
            
            //position resolution things
            tr->Branch("indexOfSmallestPosDiff", &indexOfSmallestPosDiff, "indexOfSmallestPosDiff/I");
            tr->Branch("posDiffFirstHit", &posDiffFirstHit, "posDiffFirstHit/D");
            tr->Branch("posDiffFirstTP", &posDiffFirstTP, "posDiffFirstTP/D");
            tr->Branch("posDiffDistanceOpDet", &posDiffDistanceOpDet, "posDiffDistanceOpDet/D");
            tr->Branch("posDiffDistanceXFirstHit", &posDiffDistanceXFirstHit, "posDiffDistanceXFirstHit/D");
            tr->Branch("posDiffDistanceXFirstTP", &posDiffDistanceXFirstTP, "posDiffDistanceXFirstTP/D");
            tr->Branch("posDiffYFirstTP", &posDiffYFirstTP, "posDiffYFirstTP/D");
            tr->Branch("posDiffZFirstTP", &posDiffZFirstTP, "posDiffZFirstTP/D");
            tr->Branch("posDiffDist3D", &posDiffDist3D, "posDiffDist3D/D");
            tr->Branch("posDiffUltimate", &posDiffUltimate, "posDiffUltimate/D");
            
            //primary truth information
            tr->Branch("nuEnergy", &nuEnergy, "nuEnergy/D");
            tr->Branch("summedTruthEnergyTotal", &summedTruthEnergyTotal, "summedTruthEnergyTotal/D");
            tr->Branch("numNeutrons", &numNeutrons, "numNeutrons/D");
            tr->Branch("summedTruthEnergyNeutron", &summedTruthEnergyNeutron, "summedTruthEnergyNeutron/D");
            tr->Branch("numProtons", &numProtons, "numProtons/D");
            tr->Branch("summedTruthEnergyProton", &summedTruthEnergyProton, "summedTruthEnergyProton/D");
            tr->Branch("numElectrons", &numElectrons, "numElectrons/D");
            tr->Branch("summedTruthEnergyElectron", &summedTruthEnergyElectron, "summedTruthEnergyElectron/D");
            tr->Branch("numGammas", &numGammas, "numGammas/D");
            tr->Branch("summedTruthEnergyGamma", &summedTruthEnergyGamma, "summedTruthEnergyGamma/D");
            tr->Branch("meanTruthEnergyGamma", &meanTruthEnergyGamma, "meanTruthEnergyGamma/D");
            tr->Branch("numTotalGammas", &numTotalGammas, "numTotalGammas/D");
            
        }

    ChargeDist::~ChargeDist(){}
    void ChargeDist::analyze(const art::Event& event){
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////// SETUP ///////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
        //Get the property services
        auto const* detectorProperties =
            lar::providerFrom< detinfo::DetectorPropertiesService >();
        auto const* detectorClocks = lar::providerFrom< detinfo::DetectorClocksService >();
        //grab electron lifetime, drift velocity
        //August 7 2018: tau = 3000.0 us
        float tau = detectorProperties->ElectronLifetime(); //microseconds    
        //August 20 2018: driftVelocity = 0.160563 cm/us    
        float driftVelocity = detectorProperties->DriftVelocity();
        
        //get neutrino truth energy
		art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
		std::vector<art::Ptr<simb::MCTruth>> mclist;
		//if(event.getByLabel("generator", mctruthListHandle)) 
		//	art::fill_ptr_vector(mclist, mctruthListHandle);
		if(event.getByLabel("marley", mctruthListHandle)) 
			art::fill_ptr_vector(mclist, mctruthListHandle);
			
		for(unsigned int iList = 0; iList < mclist.size(); ++iList)
			nuEnergy = mclist[iList]->GetNeutrino().Nu().E()*1000; //MeV      
			
		//std::cout << "Energy: " << nuEnergy << " " << mclist.size() << std::endl;
		
		//fill the event information here
		fEvent = event.id().event();
                    
    ////////////////////////////////////////////////////////////////////////
    /////////////////////// TRUTH DRIFT CORRECTION /////////////////////////
    ////////////////////////////////////////////////////////////////////////
        auto particleHandle
            = event.getValidHandle<std::vector<simb::MCParticle>>(MCParticleLabel);
                    
        for(simb::MCParticle const& part: *particleHandle) { 
			//primary electron
            if(prim.compare(part.Process()) == 0 && part.PdgCode() == 11){
                const TLorentzVector& PositionStart = part.simb::MCParticle::Position(0);
                driftTime = std::abs(PositionStart.X()/driftVelocity);
                truthXPos = PositionStart.X();
                truthYPos = PositionStart.Y();
                truthZPos = PositionStart.Z();
                //std::cout << "Primary electron starts at (" << truthXPos << ", " << truthYPos << ", " << truthZPos << ")" << std::endl;
                
            }
            
        }
        //define vector to hold vertex
		TVector3 truth;
		truth.SetXYZ(truthXPos, truthYPos, truthZPos);
              
    ////////////////////////////////////////////////////////////////////////
    //////////////////////// CHARGE DISTRIBUTION ///////////////////////////
    ////////////////////////////////////////////////////////////////////////
        //float totalCharge = rawDigitChargeDist(event, fDetSimProducerLabel);
        charge = rawDigitChargeDist(event, fDetSimProducerLabel);
        
        std::cout << "Event " << event.id().event() << ": " << charge << std::endl;
        
        //get geometry service
        auto const& geometry(*lar::providerFrom< geo::Geometry >());
		// Get a raw waveform handle
		// (it behaves like a pointer to a vector of waveforms) from the event
		art::ValidHandle< std::vector< raw::RawDigit > > rawDigitHandle
		 = event.getValidHandle< std::vector< raw::RawDigit > >(fDetSimProducerLabel);

        // Only use collection wires
        unsigned short wirePlane = geo::kZ;
        
        std::vector<double> pedestals;
        std::vector<double> rawChargesPos;
        std::vector<double> rawChargesNeg;
        
        // Loop over all wires
        for (auto const &rawDigit : (*rawDigitHandle))
        {
            
          // Have to uncompress RawDigit first in case it is compressed
          raw::RawDigit::ADCvector_t adcvec(rawDigit.Samples());
          raw::Uncompress(rawDigit.ADCs(), adcvec, rawDigit.GetPedestal(), rawDigit.Compression());
          
          //define the channel stuff here to use in signal loop
          unsigned int channelNumber = rawDigit.Channel();
          std::vector< geo::WireID > wireID = geometry.ChannelToWire(channelNumber);
                                      
          // Subtract the pedestal
          for (short const& signal : adcvec){
            // ADC counts are non-negative
              if (signal > 0){ 
				
				float sigToConsider = static_cast< float >(signal - rawDigit.GetPedestal());
                
				for (auto const& wireIter : wireID)
					if (wireIter.Plane == wirePlane){
						pedestals.emplace_back(rawDigit.GetPedestal());
						
						if(sigToConsider > 0) rawChargesPos.emplace_back(sigToConsider);
						else if(sigToConsider < 0) rawChargesNeg.emplace_back(sigToConsider);
					}
                
              }
        }//end loop over signals

        }//end loop over raw digits
        
        numRawDigits = rawDigitHandle->size();
        numPositiveSignals = rawChargesPos.size();
        numNegativeSignals = rawChargesNeg.size();
        
        
        //get events for mcc11 charge exploratory study
        //Double_t chargeCut = 2000; //4.25
        //Double_t chargeCut = 1000; //5.25
        //Double_t chargeCut = 0; //6.25
        //Double_t chargeCut = 7500; //7.25
        //Double_t chargeCut = 2000; //8.25
        //Double_t chargeCut = 5000; //10.75
        //Double_t chargeCut = 7500; //12.75
        //Double_t chargeCut = 6000; //14.75
        //Double_t chargeCut = 5000; //16.75
        
        //if(charge < chargeCut) std::cout << "Event " << event.id().event() << " is below the charge cut: " << charge << " < " << chargeCut << std::endl;
        //else std::cout << "Event " << event.id().event() << " is above the charge cut: " << charge << " > " << chargeCut << std::endl;
        
        //apply thresholds
        charge_1MeVThresh = rawDigitChargeDistThreshold(event, fDetSimProducerLabel, threshold_1MeV);
        charge_15MeVThresh = rawDigitChargeDistThreshold(event, fDetSimProducerLabel, threshold_15MeV);
        charge_2MeVThresh = rawDigitChargeDistThreshold(event, fDetSimProducerLabel, threshold_2MeV);
        
    ////////////////////////////////////////////////////////////////////////
    /////////////////////// RECOB DRIFT CORRECTION /////////////////////////
    ////////////////////////////////////////////////////////////////////////
        photonFlashTime = getPhotonFlashTime(event, OpFlashLabel);
          
        //declare variales
        std::vector<double> hitPeakTimes;
        std::vector<geo::WireID> hitWireIDs;
        std::vector<double> allHitPeakTimes;
        std::vector<geo::WireID> allHitWireIDs;
        std::vector<double> hitStartTicks;
            
        //get track information
        auto trackListHandle = event.getValidHandle<std::vector<recob::Track>>(HitTrackLabel);
        art::FindMany<recob::Hit> fmh(trackListHandle, event, HitTrackLabel);
          
        std::vector<double> trajTruthDists;
         
        //loop through track objects
        for(size_t iTrack = 0; iTrack < trackListHandle->size(); ++iTrack){
            art::Ptr<recob::Track> trPtr(trackListHandle, iTrack);
            
            //find smallest dist between 1st traj point of each track
            auto trajAll = trPtr->TrajectoryPoint(0);
            TVector3 trajAllPos;
            trajAllPos.SetXYZ(trajAll.position.X(), trajAll.position.Y(), trajAll.position.Z());
            
            TVector3 diff = trajAllPos - truth;
            double dist = diff.Mag();
            trajTruthDists.emplace_back(dist);
            //get the hits associated with every track in event
            std::vector<const recob::Hit*> allHits = fmh.at(iTrack);
            //get ALL the PeakTimes associated with the hits
            for(std::size_t iHit = 0; iHit < allHits.size(); ++iHit){
                //std::cout << "Track #" << iTrack << ", Hit #" << iHit << ": " << allHits[iHit]->StartTick() << ", " << allHits[iHit]->PeakTime() << std::endl;
                allHitPeakTimes.emplace_back(allHits[iHit]->PeakTime());
                allHitWireIDs.emplace_back(allHits[iHit]->WireID());
            }
            
            //std::cout << "Track #" << iTrack << ": " << *std::min_element(testPT.begin(), testPT.end()) << std::endl;
            
            //looking at the first track in the event for primary electron
            if(iTrack == (unsigned int)0){
                
                //POS DIFF: FIRST TRAJECTORY POINT IN FIRST TRACK
                auto traj = trPtr->TrajectoryPoint(0);
                double x = traj.position.X();
                
                recoTP_XPos = traj.position.X();
                recoTP_YPos = traj.position.Y();
                recoTP_ZPos = traj.position.Z();
                
                //find the distance to the nearest photon detector
                auto const& geometry(*lar::providerFrom< geo::Geometry >());
                unsigned int closestChannel = geometry.GetClosestOpDet(traj.position);
                auto const& closestOpDet = geometry.OpDetGeoFromOpChannel(closestChannel);
                posDiffDistanceOpDet = closestOpDet.DistanceToPoint(traj.position);
                
                posDiffDistanceXFirstTP = std::abs(x);
                
                posDiffYFirstTP = std::abs(traj.position.Y()) - std::abs(truthYPos);
                posDiffZFirstTP = std::abs(traj.position.Z()) - std::abs(truthZPos);
                
                TVector3 recob;
                recob.SetXYZ(traj.position.X(), traj.position.Y(), traj.position.Z());
                TVector3 diff = recob-truth;
                posDiffDist3D = diff.Mag();
                
                //std::cout << "Event " << event.id().event() << ": " << posDiffDist3D << " for (" << truthXPos << ", " << truthYPos << ", " << truthZPos << ") vs. (" << traj.position.X() << ", " << traj.position.Y() << ", " << traj.position.Z() << ")" << std::endl;
                
                //std::cout << closestOpDet.GetCenter().X() << ", " << closestOpDet.GetCenter().Y() << ", " << closestOpDet.GetCenter().Z()<< std::endl;
                //std::cout << traj.position.X() << ", " << traj.position.Y() << ", " << traj.position.Z() << std::endl;
                
                posDiffFirstTP = std::abs(x) - std::abs(truthXPos);
                
                //std::cout << posDiffFirstTP << " vs. " << posDiff3D << std::endl;
                
                //get hits associated with first track
                std::vector<const recob::Hit*> hits = fmh.at(iTrack);
                for(std::size_t iHit = 0; iHit < hits.size(); ++iHit){
                    hitPeakTimes.emplace_back(hits[iHit]->PeakTime()); //units ticks
                    hitWireIDs.emplace_back(hits[iHit]->WireID());
                    
                    auto tmp = hits[iHit]->StartTick();
                    
                    hitStartTicks.emplace_back( (double)tmp );
                }
            }
        }//end loop through tracks in event
          
        //check that there IS a track/hits associated with track
        if(hitPeakTimes.size() > 0 && allHitPeakTimes.size() > 0){
            //std::cout << "COMPARE TO " << *std::min_element(allHitPeakTimes.begin(), allHitPeakTimes.end()) << std::endl;
            
            //FIRST HIT IN FIRST TRACK
            double wireDepoTicks3 = hitPeakTimes[0];
            geo::PlaneID pid3(hitWireIDs[0].Cryostat, hitWireIDs[0].TPC, hitWireIDs[0].Plane);
            double driftRecobTime3 = detectorClocks->TPCTick2Time(wireDepoTicks3) - photonFlashTime; //units us
            
            //also define driftRecobTime here
            driftRecobTime = detectorClocks->TPCTick2Time(wireDepoTicks3) - photonFlashTime;
            
            double driftRecobTicks3 = detectorClocks->TPCG4Time2Tick(driftRecobTime3*1000); //units ticks
            double xPosRecob3 = detectorProperties->ConvertTicksToX(driftRecobTicks3, pid3);
            posDiffDistanceXFirstHit = std::abs(xPosRecob3);
            
            posDiffFirstHit = std::abs(xPosRecob3)-std::abs(truthXPos);
            
            recoXPos = xPosRecob3;
            
            //SMALLEST DIFF IN FIRST HIT OUT OF ALL TRACKS
            std::vector<double> posDiffToCheck = getPositionDifferences(event, HitTrackLabel, photonFlashTime, truthXPos);
            //okay so now I have all the pos diff for all the tracks in the event
            //now I think I'll output in a branch the track with the smallest abs(pos diff)
            std::vector<double> absPosDiff;
            for(std::size_t i = 0; i < posDiffToCheck.size(); ++i){
                absPosDiff.emplace_back(abs(posDiffToCheck[i]));
            }
            indexOfSmallestPosDiff = std::distance(absPosDiff.begin(), std::min_element( absPosDiff.begin(), absPosDiff.end() ) );
            
            //if index is 0, then from first track
            if(indexOfSmallestPosDiff == 0){
                //now need to see whether the track is flipped or not
                //how would we determine that? if the pos diff is small enough?
                //could also find position diff from end, see if better
                
                double wireDepoTicks4 = hitPeakTimes[hitPeakTimes.size()-1];
                geo::PlaneID pid4(hitWireIDs[hitPeakTimes.size()-1].Cryostat, hitWireIDs[hitPeakTimes.size()-1].TPC, hitWireIDs[hitPeakTimes.size()-1].Plane);
                double driftRecobTime4 = detectorClocks->TPCTick2Time(wireDepoTicks4) - photonFlashTime; //units us
                double driftRecobTicks4 = detectorClocks->TPCG4Time2Tick(driftRecobTime4*1000); //units ticks
                double xPosRecob4 = detectorProperties->ConvertTicksToX(driftRecobTicks4, pid4);
                double posDiffLastHit = std::abs(xPosRecob4)-std::abs(truthXPos);
                
                //if first hit is smaller, then closer to front
                if( abs(posDiffFirstHit) <= abs(posDiffLastHit) ){
                    posDiffUltimate = posDiffFirstHit;
                }
                
            }

            
        }//check to make sure there are hits/tracks assn
        
    ////////////////////////////////////////////////////////////////////////
    /////////////////////// MC TRUTH INFORMATION ///////////////////////////
    ////////////////////////////////////////////////////////////////////////

		//now I need to loop through all the MC particles in the event
        std::vector<int> PDGcodes;  
        for(simb::MCParticle const& part: *particleHandle) 
			PDGcodes.emplace_back(part.PdgCode());  
        //sort and delete duplicates
        sort( PDGcodes.begin(), PDGcodes.end() );
        PDGcodes.erase( unique( PDGcodes.begin(), PDGcodes.end() ), PDGcodes.end() );
        
        double totalEnergy = 0.0;
        for(unsigned int i = 0; i < PDGcodes.size(); i++){
            std::vector<double> energies = getListOfUniqueMCParticles(event, MCParticleLabel, PDGcodes[i]);
            
            double tmp = std::accumulate(energies.begin(), energies.end(), 0.0);
            
            //std::cout << PDGcodes[i] << " " << tmp << std::endl;
            
            totalEnergy += tmp;
        }
        summedTruthEnergyTotal = totalEnergy;
            
        std::vector<double> uniqueNeutrons = getListOfUniqueMCParticles(event, MCParticleLabel, 2112);
        numNeutrons = uniqueNeutrons.size();
        summedTruthEnergyNeutron = std::accumulate(uniqueNeutrons.begin(), uniqueNeutrons.end(), 0.0);
        
        std::vector<double> uniqueProtons = getListOfUniqueMCParticles(event, MCParticleLabel, 2212);
        numProtons = uniqueProtons.size();
        summedTruthEnergyProton = std::accumulate(uniqueProtons.begin(), uniqueProtons.end(), 0.0);
        
        std::vector<double> uniqueElectrons = getListOfUniqueMCParticles(event, MCParticleLabel, 11);
        numElectrons = uniqueElectrons.size();
        summedTruthEnergyElectron = std::accumulate(uniqueElectrons.begin(), uniqueElectrons.end(), 0.0);
        
        std::vector<double> uniqueGammas = getListOfUniqueMCParticles(event, MCParticleLabel, 22);
        numGammas = uniqueGammas.size();
        summedTruthEnergyGamma = std::accumulate(uniqueGammas.begin(), uniqueGammas.end(), 0.0);
          
        if(numGammas > 0) meanTruthEnergyGamma = summedTruthEnergyGamma/numGammas;
        else meanTruthEnergyGamma = 0;

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////// FILL HISTS, TREE ////////////////////////////
    ////////////////////////////////////////////////////////////////////////
        
        //add values to tree
        //charge = totalCharge;
        chargeCorrect = charge*exp(driftTime/tau);
        chargeRecobCorrect = charge*exp(driftRecobTime/tau);
        
        //std::cout << "Drift Correct: " << charge << " * exp( " << driftTime << " / " << tau << " )" << std::endl;
        //std::cout << "Drift Time: " << driftTime << " = " << truthXPos << " / " << driftVelocity << std::endl;
        
        double totalCorrectFactor = driftRecobTime/tau;
        double factor25 = totalCorrectFactor*0.25;
        double factor50 = totalCorrectFactor*0.50;
        double factor75 = totalCorrectFactor*0.75;
        chargeRecob25PerCentCorrect = charge*exp(factor25);
        chargeRecob50PerCentCorrect = charge*exp(factor50);
        chargeRecob75PerCentCorrect = charge*exp(factor75);
        
        DCFactor_Reco = exp(driftRecobTime/tau);
        FractionalDC_Reco = charge/chargeRecobCorrect;
        DCFactor_Truth = exp(driftTime/tau);
        FractionalDC_Truth = charge/chargeCorrect;
        
        //efficiency factor
        if(trackListHandle->size() > 0) efficiencyCounter = 1;
        else efficiencyCounter = 0;
        
        //fill no matter what?
        //worried that this will be filled with duplicate information
        tr->Fill();
        
        //set charge to some null value to check
        charge = -99999.0;
        
        //only fill tree if there is appropriate reconstruction information
        //if(trackListHandle->size() > 0) tr->Fill();
        //else std::cout << "Event " << event.id().event() << ": no track objects. Tree not filled." << std::endl;
        
        
      }//end main code
      
    ////////////////////////////////////////////////////////////////////////
    ///////////////////////////// FUNCTIONS ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    
	float ChargeDist::rawDigitChargeDist(const art::Event& event, std::string labelToGetRawDigits){
		//get geometry service
        auto const& geometry(*lar::providerFrom< geo::Geometry >());
		// Get a raw waveform handle
		// (it behaves like a pointer to a vector of waveforms) from the event
		art::ValidHandle< std::vector< raw::RawDigit > > rawDigitHandle
		 = event.getValidHandle< std::vector< raw::RawDigit > >(labelToGetRawDigits);

        // Only use collection wires
        unsigned short wirePlane = geo::kZ;

        // Total charge collected in the event
        float totalCharge = 0.0;
        
        // Loop over all wires
        for (auto const &rawDigit : (*rawDigitHandle))
        {
          // Total charge on this wire
          float rawCharge = 0.0;
            
          // Have to uncompress RawDigit first in case it is compressed
          raw::RawDigit::ADCvector_t adcvec(rawDigit.Samples());
          raw::Uncompress(rawDigit.ADCs(), adcvec, rawDigit.GetPedestal(), rawDigit.Compression());
          
          
          // Subtract the pedestal
          for (short const& signal : adcvec){
            // ADC counts are non-negative
              if (signal > 0){ 
			//if (signal > rawDigit.GetPedestal()){

                rawCharge += static_cast< float >(signal - rawDigit.GetPedestal());
              }
        }//end loop over signals

          // Check whether this is a collection wire
          // If it is, add its integrated charge to the totalCharge variable
          unsigned int channelNumber = rawDigit.Channel();
          std::vector< geo::WireID > wireID =
                                      geometry.ChannelToWire(channelNumber);
          for (auto const& wireIter : wireID)
            if (wireIter.Plane == wirePlane) totalCharge += rawCharge;
            
        }//end loop over raw digits
        
        return totalCharge;

    }
        
        float ChargeDist::rawDigitChargeDistThreshold(const art::Event& event, std::string labelToGetRawDigits, float threshold){
            
            //get geometry service
            auto const& geometry(*lar::providerFrom< geo::Geometry >());
            // Get a raw waveform handle
            // (it behaves like a pointer to a vector of waveforms) from the event
            art::ValidHandle< std::vector< raw::RawDigit > > rawDigitHandle
            = event.getValidHandle< std::vector< raw::RawDigit > >
            (labelToGetRawDigits);
            
            // Only use collection wires
            unsigned short wirePlane = geo::kZ;
            
            // Total charge collected in the event
            float totalCharge = 0.0;
            
            // Loop over all wires
            for (auto const &rawDigit : (*rawDigitHandle))
            {
                // Total charge on this wire
                float rawCharge = 0.0;
                
                // Have to uncompress RawDigit first in case it is compressed
                raw::RawDigit::ADCvector_t adcvec(rawDigit.Samples());
                //raw::Uncompress(rawDigit.ADCs(), adcvec, rawDigit.Compression());
                
                raw::Uncompress(rawDigit.ADCs(), adcvec, rawDigit.GetPedestal(), rawDigit.Compression());
                
                // Subtract the pedestal
                for (short const& signal : adcvec){
                    // ADC counts are non-negative
                    if (signal > 0){
                        // if (signal > rawDigit.GetPedestal()){
                        
                        rawCharge += static_cast< float >(signal - rawDigit.GetPedestal());
                    }
                }//end loop over signals
                
                
                //if(rawCharge > threshold)std::cout << rawCharge << " " << threshold << std::endl;
                
                // Check whether this is a collection wire
                // If it is, add its integrated charge to the totalCharge variable
                unsigned int channelNumber = rawDigit.Channel();
                std::vector< geo::WireID > wireID =
                geometry.ChannelToWire(channelNumber);
                for (auto const& wireIter : wireID)
                    //if (wireIter.Plane == wirePlane) totalCharge += rawCharge;
                    if (wireIter.Plane == wirePlane && rawCharge > threshold) totalCharge += rawCharge;
                
            }//end loop over raw digits
            
            return totalCharge;
        }

    double ChargeDist::getPhotonFlashTime(const art::Event& event, std::string OpFlashLabel){
        double photonFlashTime = 0.0;
        
        //handle for the opFlash
        auto opflashHandle = event.getValidHandle<std::vector<recob::OpFlash>>(OpFlashLabel);
        std::vector<double> totalPEs; std::vector<float> flashtimes;
        for (std::size_t iFlash = 0; iFlash < opflashHandle->size(); ++iFlash) {
            art::Ptr<recob::OpFlash> opFlashPtr(opflashHandle, iFlash);
            totalPEs.emplace_back(opFlashPtr->TotalPE());
            flashtimes.emplace_back(opFlashPtr->Time()); //units us
            
            //std::cout << iFlash << ": " << opFlashPtr->Time() << " us " << std::endl;
        }
        //find max opflash
        if(opflashHandle->size() > 0){
            int index = std::distance(totalPEs.begin(), std::max_element(totalPEs.begin(), totalPEs.end()));
            photonFlashTime = flashtimes[index];
        }
        
        return photonFlashTime;
    }

    std::vector<double> ChargeDist::getListOfUniqueMCParticles(const art::Event& event, std::string MCParticleLabel, int PDGcode){
        auto particleHandle
        = event.getValidHandle<std::vector<simb::MCParticle>>(MCParticleLabel);
        
        std::vector<double> truthEnergy;
        std::vector<std::string> beginProcess;
        std::vector<std::string> endProcess;
        std::vector<int> trackIDs;
        std::vector<int> motherIDs;
        
        for(simb::MCParticle const& part: *particleHandle) {
            //if(part.PdgCode() == PDGcode){
            if(part.PdgCode() == PDGcode && prim.compare(part.Process()) == 0){
                
                double KE = 1000*(part.E() - part.Mass());
                
                truthEnergy.emplace_back(KE);
                beginProcess.emplace_back(part.Process());
                endProcess.emplace_back(part.EndProcess());
                trackIDs.emplace_back(part.TrackId());
                motherIDs.emplace_back(part.Mother());
                
            }
        }
        
        std::vector<double> uniqueParticles;
        for(unsigned int i = 0; i < motherIDs.size(); i++){
            int motherToCheck = motherIDs[i];
            //see if this mother ID is in the list of track IDs
            if( std::find( trackIDs.begin(), trackIDs.end(), motherToCheck ) != trackIDs.end() ){
                //found it in the list of track IDs
                //check that the beginning/end process matches
                
                int index = std::distance(trackIDs.begin(), std::find( trackIDs.begin(), trackIDs.end(), motherToCheck ) );
                //this index is the track ID of the mother particle
                //i tells you the track ID of the daughter particle
                
                std::string beginToCheck = beginProcess[i];
                std::string endToCheck = endProcess[index];
                
                if(beginToCheck.compare(endToCheck) == 0){
                    //nothing
                }
                else{
                    //if the processes don't match, then probably new neutron
                    //need to keep the index so we can sum up truth energy
                    uniqueParticles.emplace_back(i);
                }
            }
            else{
                //did not find the mother ID in list of track IDs
                //so keep the index so we can sum up the truth energy
                uniqueParticles.emplace_back(i);
            }
        }
        
        std::vector<double> uniqueEnergies;
        for(unsigned int i = 0; i < uniqueParticles.size(); i++){
            int ind = uniqueParticles[i];
            double energy = truthEnergy[ind];
            uniqueEnergies.emplace_back(energy);
            
        }
        
        return uniqueEnergies;
        
    }

    std::vector<double> ChargeDist::getPositionDifferences(const art::Event& event, std::string labelToGetTracks, double PhotonFlashTime, double truthXPos){
        auto const* detectorProperties =
          lar::providerFrom< detinfo::DetectorPropertiesService >();
        auto const* detectorClocks = lar::providerFrom< detinfo::DetectorClocksService >();
     
        std::vector<double> positionDifferences;
        auto trackListHandle = event.getValidHandle<std::vector<recob::Track>>(labelToGetTracks);
        
        art::FindMany<recob::Hit> fmh(trackListHandle, event, labelToGetTracks);
        
        std::vector<std::vector<double>> hitPeakTimes;
        std::vector<std::vector<geo::WireID>> hitWireIDs;
        
        //loop through track objects
        for(size_t iTrack = 0; iTrack < trackListHandle->size(); ++iTrack){
            art::Ptr<recob::Track> tkPtr(trackListHandle, iTrack);
            
            //get hits associated with track
            std::vector<const recob::Hit*> allHits = fmh.at(iTrack);
            
            //define vector to hold peak times for the hits in this track
            std::vector<double> tmpHitPeakTimes;
            std::vector<geo::WireID> tmpHitWireIDs;
            
            //loop through the hits, grabbing the peaktimes
            for(std::size_t iHit = 0; iHit < allHits.size(); ++iHit){
                tmpHitPeakTimes.emplace_back(allHits[iHit]->PeakTime());
                tmpHitWireIDs.emplace_back(allHits[iHit]->WireID());
            }
            
            //now, add the tmp vector to the 2D vector
            hitPeakTimes.emplace_back(tmpHitPeakTimes);
            hitWireIDs.emplace_back(tmpHitWireIDs);
            
        }
        
        //check to make sure that there are tracks
        if(hitPeakTimes.size() > 0){
            for(std::size_t i = 0; i < hitPeakTimes.size(); ++i){
                //okay, now we grab the peak time
                double firstHitPeakTime = hitPeakTimes[i][0];
                geo::PlaneID pid(hitWireIDs[i][0].Cryostat, hitWireIDs[i][0].TPC, hitWireIDs[i][0].Plane);
                
                //now calculate the pos diff
                double driftRecobTime = detectorClocks->TPCTick2Time(firstHitPeakTime) - photonFlashTime; //units us
                double driftRecobTicks = detectorClocks->TPCG4Time2Tick(driftRecobTime*1000);
                double xPosRecob = detectorProperties->ConvertTicksToX(driftRecobTicks, pid);
                
                double posDiff = std::abs(xPosRecob) - std::abs(truthXPos);
                
                positionDifferences.emplace_back(posDiff);
                
            }
            
        }
        
        return positionDifferences;
    }


    }

