/*
 * \file: CalibrationTree_module.cc
 * \author: JStock (jason.stock@mines.sdsmt.edu)
 * \brief: This is a small analysis Tree made for use in the Calibration group to investigate radiological backgrounds and calibration sources.
 *
 */

#include "CalibrationTree.h"

namespace{}

namespace CalibrationTree {

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  CalibrationTree::CalibrationTree(fhicl::ParameterSet const& pSet)
    :EDAnalyzer(pSet),
    private_HitLabel(pSet.get<art::InputTag>("HitLabel","gausshit")),
    private_OpHitLabel(pSet.get<art::InputTag>("OpHitLabel","ophit"))
  {  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  /*
  CalibrationTree::CalibrationTree(fhiclConfig const& config)
    :EDAnalyzer(config),
    private_HitLabel(config.HitLabel()),
    private_OpHitLabel(config.OpHitLabel())
    {  }
    */ //Commented out until I have a subtable for EDAnalyzer

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void CalibrationTree::FillJHRV(int trackId, art::Ptr<recob::Hit> hit, int hitIndex, double split){ //Factor out everything that doesn't require hit (all the trackID stuff) into it's own function (as it is duplicated in CalibrationTree::FillJHRV for ophits) or turn this into a template.
    const simb::MCParticle* tIdToPart = PIS->TrackIdToParticle_P(trackId);
    const simb::MCParticle* eve       = PIS->TrackIdToParticle_P(PIS->TrackIdToEveTrackId(trackId));
    //bool findPart = [&toFind](auto& a){ return a.particle==tIdToPart;}
    //bool findEve  = [](auto& a){ return a.particle==eve;}
    std::vector<jointHitRecord>::iterator partRec = std::find(jhrV.begin(), jhrV.end(), jointHitRecord(tIdToPart));
    std::vector<jointHitRecord>::iterator eveRec = std::find(jhrV.begin(), jhrV.end(), jointHitRecord(eve));
    //std::vector<jointHitRecord>::iterator partRec = std::find_if(jhrV.begin(), jhrV.end(), jointHitRecord(tIdToPart));
    //std::vector<jointHitRecord>::iterator eveRec = std::find_if(jhrV.begin(), jhrV.end(), jointHitRecord(eve));
    auto epos = std::distance(jhrV.begin(),eveRec); 
    if(partRec==jhrV.end()){ //make new Entry
      if(eveRec==jhrV.end()){ //Eve must also be added.
        //auto itr = jhrV.emplace_back(tIdToPart);
        jhrV.emplace_back(tIdToPart);
        jhrV.at(jhrV.size()-1).eveIndex=jhrV.size();
        jhrV.emplace_back(eve);
      }else{//Eve exists, part does not. Make part, assign correct eve.
        jhrV.emplace_back(tIdToPart);
        jhrV.at(jhrV.size()-1).eveIndex=epos;
      }
    }else{//Part exists. Does part know it has an eve? If no, tell it where it's eve lives.
      if( !(partRec->eveIndex) && (eveRec!=jhrV.end() ) ){
        partRec->eveIndex = epos;
      }else if( !(partRec->eveIndex) && (eveRec!=jhrV.end() ) ){ 
        partRec->eveIndex=jhrV.size();
        jhrV.emplace_back(eve);
      }else{//particle exists and knows it's eve
        //Only care if it doesn't agree with the eve we just found.
        if(partRec->eveIndex!=epos){throw;}
      }
    }
    //Do the actual filling.
    partRec = std::find(jhrV.begin(), jhrV.end(), tIdToPart);//update partRec (in case we added the particle ourselves.
    if(partRec==jhrV.end()){throw;} //This should never happen by this point, as I explicitly added it if it didn't exist.
    partRec->AddHit(hit, split);
  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void CalibrationTree::FillJHRV(int trackId, art::Ptr<recob::OpHit> hit, int ophitIndex, double split){
    const simb::MCParticle* tIdToPart = PIS->TrackIdToParticle_P(trackId);
    const simb::MCParticle* eve       = PIS->TrackIdToParticle_P(PIS->TrackIdToEveTrackId(trackId));
//    bool findPart = [](auto& a){ return a.particle==tIdToPart;}
//    bool findEve  = [](auto& a){ return a.particle==eve;}
    std::vector<jointHitRecord>::iterator partRec = std::find(jhrV.begin(), jhrV.end(), tIdToPart);
    std::vector<jointHitRecord>::iterator eveRec = std::find(jhrV.begin(), jhrV.end(), eve);
    auto epos = std::distance(jhrV.begin(),eveRec); 
    if(partRec==jhrV.end()){ //make new Entry
      if(eveRec==jhrV.end()){ //Eve must also be added.
        jhrV.emplace_back(tIdToPart);
        jhrV.at(jhrV.size()-1).eveIndex=jhrV.size();
        jhrV.emplace_back(eve);
      }else{//Eve exists, part does not. Make part, assign correct eve.
        jhrV.emplace_back(tIdToPart);
        jhrV.at(jhrV.size()-1).eveIndex=epos;
      }
    }else{//Part exists. Does part know it has an eve? If no, tell it where it's eve lives.
      if( !(partRec->eveIndex) && (eveRec!=jhrV.end() ) ){
        partRec->eveIndex = epos;
      }else if( !(partRec->eveIndex) && (eveRec!=jhrV.end() ) ){ 
        partRec->eveIndex=jhrV.size();
        jhrV.emplace_back(eve);
      }else{//particle exists and knows it's eve
        //Only care if it doesn't agree with the eve we just found.
        if(partRec->eveIndex!=epos){throw;}
      }
    }
    //Do the actual filling.
    partRec = std::find(jhrV.begin(), jhrV.end(), tIdToPart);//update partRec (in case we added the particle ourselves.
    if(partRec==jhrV.end()){throw;} //This should never happen by this point, as I explicitly added it if it didn't exist.
    partRec->AddHit(hit, split);
  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void CalibrationTree::beginJob(){
    private_CalibrationTree = private_service_tfs->make<TTree>("CalibrationTree", "A TTree for short records of the hits and ophits produced by given particles in the simulation");
    private_CalibrationParticle = private_CalibrationTree->Branch("particles", &private_particleBuffer, "pid/i:eveIndex/i:event/i:subrun/i:run/i");
    //private_CalibrationCharge = private_CalibrationTree->Branch(private_CalibrationParticle, "hits",&private_hitBuffer,"charge/d:time/d:width/d:wire/i:hitIndex/i");
    //private_CalibrationLight = private_CalibrationTree->Branch(private_CalibrationParticle, "ophits", &private_lightBuffer, "pes/d:time/d:width/d:opchan/i:opdet/i:opHitIndex/i");
    private_CalibrationCharge = private_service_tfs->make<TBranch>(private_CalibrationParticle, "hits",&private_hitBuffer,"charge/d:time/d:width/d:wire/i:hitIndex/i");
    private_CalibrationLight = private_service_tfs->make<TBranch>(private_CalibrationParticle, "ophits", &private_lightBuffer, "pes/d:time/d:width/d:opchan/i:opdet/i:opHitIndex/i");
  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void CalibrationTree::analyze(const art::Event& evt) {
    //Set up BackTracking //Now in header
/*    art::ServiceHandle<cheat::ParticleInventoryService> PIS;
    art::ServiceHandle<cheat::BackTrackerService> BTS;
    art::ServiceHandle<cheat::PhotonBackTrackerService> PBS;*/

    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitList;
    if(evt.getByLabel(private_HitLabel,hitHandle) )
      art::fill_ptr_vector(hitList, hitHandle);
    art::Handle<std::vector<recob::OpHit>> opHitHandle;
    std::vector<art::Ptr<recob::OpHit>> opHitList;
    if(evt.getByLabel(private_OpHitLabel,opHitHandle) )
      art::fill_ptr_vector(opHitList, opHitHandle);

    unsigned int hitCounter=0;
    for ( auto const& hit : hitList){
      hitCounter++;
      std::vector<const sim::IDE*> IDEs_Ps 
        = BTS->HitToSimIDEs_Ps(hit);
      if(IDEs_Ps.empty()){
        //I don't know what this would be. Maybe make a dummy particle? A dummy jhr?
        //I am leaning toward a dummy JHR that just holds all of the hits that aren't from tracks. This singleton doesn't need to be in the map, but can be handeled as a special case. If I do want it in the map, I have to make a dummy MCParticle.
        //That way nothing is lost and noise studies are possible.
      }else{//This needs to work more like TrackIDEs, but really it needs to be particleIDEs, with appropriate particle energy fractions. Perhaps I should add this to the BackTracker (can't add to SimChannels. Too much work, breaking change. But doing it in the backtracker should be fine.
        for( auto const& ide_P : IDEs_Ps){
          //This can still be a fake hit. Include a check for an unbelievably low amount of energy deposited (if the hit is noise but happened to corespond to a stray electron that did in fact hit the wires).
          //I chose to scale by nelectrons, this can be changed later by flipping which line is commented below.
          int trackId = ide_P->trackID;
          double split = ide_P->numElectrons;
          this->CalibrationTree::FillJHRV(trackId, hit, hitCounter, split);
        }//end for ide_P in IDEs_Ps
      }//end IDEs_Ps not empty
    }//end for hit in hitlist


    hitCounter=0;
    for ( auto const& ophit : opHitList ){
      ++hitCounter;
      std::vector<const sim::SDP*> SDPs_Ps 
        = PBS->OpHitToSimSDPs_Ps(ophit);
      if(SDPs_Ps.empty()){
      }else{
        for( auto const& sdp_P : SDPs_Ps){
          int trackId = sdp_P->trackID;
          double split = sdp_P->numPhotons;
          this->CalibrationTree::FillJHRV(trackId, ophit, hitCounter, split);
        }//end for sdp_P
      }//end if SDPs_Ps not empty
    } //Add ophits to hitmap. 
    //The normalization could be broken out into it's own function easily.
    int partTracker=0;
    for ( auto jhr : jhrV){
      ++partTracker;

      double total=0.0;
      for(auto& htuple : jhr.hitsRec){
        total+=std::get<2>(htuple);
      }
      for(auto& htuple : jhr.hitsRec){
        //htuple->second/=total; //why bother fixing this in the entry. I can just fix it as it goes into the Tree. (No assignment call needed for the record.
        auto const& hit = std::get<0>(htuple);
        private_hitBuffer.update(hit->Integral(), hit->PeakTime(), hit->SigmaPeakTime(), (std::get<2>(htuple)/total), hit->Channel(), std::get<1>(htuple), partTracker);
        private_CalibrationCharge->Fill();
      }
      total=0.0;
      for(auto& htuple : jhr.ohitsRec){
        total+=std::get<2>(htuple);
      }
      for(auto& htuple : jhr.ohitsRec){
        auto const& hit = std::get<0>(htuple);
        private_lightBuffer.update(hit->PE(), hit->PeakTime(), hit->Width(), (std::get<2>(htuple)/total), hit->OpChannel(),  std::get<1>(htuple), partTracker);
        private_CalibrationLight->Fill();
      }
      //Add this entry to the TTree
    }

  } //end analyse 

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void CalibrationTree::endJob(){
  }

}

DEFINE_ART_MODULE(CalibrationTree::CalibrationTree)
