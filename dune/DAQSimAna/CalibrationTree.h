/*
 * \file: CalibrationTree.h
 * \author: Jason Stock (jason.stock@mines.sdsmt.edu
 * \brief: This is a small Tree made for use in calibration analysis.
 *
 */

#ifndef CALIBRATIONTREE_H
#define CALIBRATIONTREE_H

//Includes

//LArSoft Includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"

//FrameworkIncludes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"

//CPP includes
#include <vector>
#include <map>


namespace{}//

namespace CalibrationTree {


  class CalibrationTree : public art::EDAnalyzer
  {
    public:

      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/*      struct fhiclConfig{
        fhicl::Atom<art::InputTag> HitLabel{
          fhicl::Name("HitLabel"),
            fhicl::Comment("The module label to be used for reading Hits."),
            "gaushit"
        };// end HitLabel
        fhicl::Atom<art::InputTag> OpHitLabel{
          fhicl::Name("OpHitLabel"),
            fhicl::Comment("The module label to be used for reading OpHits."),
            "ophit"
        };//end OpHitLabel
        fhicl::Atom<art::InputTag> TrackLabel{
          fhicl::Name("TrackLabel"),
            fhicl::Comment("The module label to be used for reading Tracks."),
            "pmtrack"
        };//end TrackLabel
      };//end fhiclConfig
*/
      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      struct jointHitRecord{ //To be used in a map with Particles
        const simb::MCParticle* particle;
        //const simb::MCParticle* eve;
        int eveIndex;
        std::vector<std::tuple<art::Ptr<recob::Hit>, unsigned int, double>> hitsRec;
        std::vector<std::tuple<art::Ptr<recob::OpHit>, unsigned int, double>> ohitsRec;
        jointHitRecord(const simb::MCParticle* part )
          :particle(part)
        {}
        void AddHit(art::Ptr<recob::Hit> hit, unsigned int hitIndex, double split=1.0 ){ 
          hitsRec.push_back(std::make_tuple (hit, hitIndex, split)); 
        }
        void AddHit(art::Ptr<recob::OpHit> ophit, unsigned int hitIndex, double split=1.0){
          ohitsRec.push_back(std::make_tuple (ophit, hitIndex, split));
        }
        std::vector<std::tuple<art::Ptr<recob::Hit>, unsigned int, double>> const& GetHits(){ return hitsRec; }
        std::vector<std::tuple<art::Ptr<recob::OpHit>, unsigned int, double>> const& GetOpHits(){ return ohitsRec;}
        void clear(){
          hitsRec.clear();
          ohitsRec.clear();
        }
        friend bool operator<( jointHitRecord const& lhs, jointHitRecord const& rhs) { return lhs.particle<rhs.particle; }
        friend bool operator>(  jointHitRecord const& lhs,jointHitRecord const& rhs) { return rhs<lhs; }
        friend bool operator<=( jointHitRecord const& lhs, jointHitRecord const& rhs) { return !(lhs>rhs);}
        friend bool operator>=( jointHitRecord const& lhs, jointHitRecord const& rhs) { return !(lhs<rhs);}
        friend bool operator==( jointHitRecord const& lhs, jointHitRecord const& rhs) { return lhs.particle==rhs.particle;}
        friend bool operator!=( jointHitRecord const& lhs, jointHitRecord const& rhs) { return ! (lhs==rhs); }
      };

      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      struct particleBuffer{ 
        int pnum;
        int pid;
        int eveIndex;
        int event;
        int subrun;
        int run;
        //buffer constructor
        particleBuffer()
          :pnum(0), pid(0), eveIndex(0), event(0), subrun(0), run(0)
        { }
        particleBuffer(int mpnum, int mpid, int meveIndex, int mevent, int msubrun, int mrun)
          :pnum(mpnum), pid(mpid), eveIndex(meveIndex), event(mevent), subrun(msubrun), run(mrun)
        { }
        //buffer update
        void update(int mpnum, int mpid, int meveIndex, int mevent, int msubrun, int mrun){
          pnum=mpnum; pid=mpid; eveIndex=meveIndex; event=mevent; subrun=msubrun; run=mrun;
        }
      };

      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      struct hitBuffer{
        double charge;
        double time;
        double width;
        double split;
        int wire;
        int hitIndex;
        int pnum;
        hitBuffer()
          :charge(0),time(0), width(0), split(0), wire(0), hitIndex(0), pnum(0)
        { }
        hitBuffer(double mcharge, double mtime, double mwidth, double msplit, int mwire, int mhitIndex, int mpnum)
          :charge(mcharge), time(mtime), width(mwidth), split(msplit), wire(mwire), hitIndex(mhitIndex), pnum(mpnum)
        {}
        void update(double mcharge, double mtime, double mwidth, double msplit, int mwire, int mhitIndex, int mpnum){
          charge=mcharge; time=mtime; width=mwidth; split=msplit; wire=mwire; hitIndex=mhitIndex; pnum=mpnum;
        }
      };

      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      struct lightBuffer{
        double PEs;
        double time;
        double width;
        double split;
        int opChannel;
        int opHitIndex;
        int pnum;
        lightBuffer()
          :PEs(0),time(0),width(0), split(0),opChannel(0),  opHitIndex(0), pnum(0)
        {}
        lightBuffer(double mpes, double mtime, double mwidth, double msplit, int mopchan,  int mopindex, int mpnum)
          :PEs(mpes), time(mtime), width(mwidth), split(msplit), opChannel(mopchan),  opHitIndex(mopindex), pnum(mpnum)
        {}
        void update(double mpes, double mtime, double mwidth, double msplit, int mopchan,  int mopindex, int mpnum) {
          PEs=mpes; time=mtime; width=mwidth; split=msplit; opChannel=mopchan; opHitIndex=mopindex; pnum=mpnum;
        }
      };

      CalibrationTree(fhicl::ParameterSet const& pSet);
      // CalibrationTree(fhiclConfig const& config); //commented out because of EDAnalyzer (too much work to make another table for now
      //CalibrationTree(fhicl::ParameterSet const& parameterSet);
      //CalibrationTree(fhiclConfig const& config);

      virtual void beginJob();
      virtual void analyze(const art::Event& evt) override;
      virtual void endJob();

        void FillJHRV(int trackId, art::Ptr<recob::Hit> hit, int hitIndex, double split);
        void FillJHRV(int trackId, art::Ptr<recob::OpHit> hit, int hitIndex, double split);


    private:
  art::ServiceHandle<cheat::ParticleInventoryService> PIS;
  art::ServiceHandle<cheat::BackTrackerService> BTS;
  art::ServiceHandle<cheat::PhotonBackTrackerService> PBS;

      std::vector<jointHitRecord> jhrV;

      art::InputTag private_HitLabel;
      art::InputTag private_OpHitLabel;

      void FillJHRV(int trackId, art::Ptr<recob::Hit> hit, double split=1.0);
      void FillJHRV(int trackId, art::Ptr<recob::OpHit> hit, double split=1.0);

      art::ServiceHandle<art::TFileService> private_service_tfs;
      TTree* private_CalibrationTree;
      TBranch* private_CalibrationParticle;
      TBranch* private_CalibrationCharge;
      TBranch* private_CalibrationLight;

      particleBuffer private_particleBuffer;
      hitBuffer      private_hitBuffer;
      lightBuffer  private_lightBuffer;



      //Also need buffers

  };// end class CalibrationTree

}

#endif//endif CALIBRATIONTREE_H
