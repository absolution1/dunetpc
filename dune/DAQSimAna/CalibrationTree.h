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
#include "art/Framework/Core/EDAnalyser.h"
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
      struct fhiclConfig{
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

      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      struct jointHitRecord{ //To be used in a map with Particles
        const simb::MCParticle* particle;
        //const simb::MCParticle* eve;
        int eveIndex;
        std::vector<std::tuple<art::Ptr<recob::Hit>, unsigned int, double>> hitsRec;
        std::vector<std::tuple<art::Ptr<recob::OpHit>, unsigned int, double>> ohitsRec;
        jointHitRecord::jointHitRecord(simb::MCParticle* part )
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
        bool operator<(const jointHitRecord& rhs) const{ return particle<rhs.particle; }
        bool operator>(const jointHitRecord& rhs) const{ return rhs<this; }
        bool operator<=(const jointHitRecord& rhs) const{ return !(this>rhs);}
        bool operator>=(const jointHitRecord& rhs) const{ return !(this<rhs);}
        bool operator==(const jointHitRecord& rhs) const{ return particle==rhs.particle;}
        bool operator!=(const jointHitRecords& rhs) const{ return ! this==rhs; }
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
        update(int mpnum, int mpid, int meveIndex, int mevent, int msubrun, int mrun)
          :pnum(mpnum), pid(mpid), eveIndex(meveIndex), event(mevent), subrun(msubrun), run(mrun)
        { }
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
        hitBuffer(double mcharge, double mtime, double mwidth, double msplit, int mwire, int mhitIndex, mpnum)
          :charge(mcharge), time(mtime), width(mwidth), split(msplit), wire(mwire), hitIndex(mhitIndex), pnum(mpnum)
        {}
        update(double mcharge, double mtime, double mwidth, double msplit, int mwire, int mhitIndex, int mpnum)
          :charge(mcharge), time(mtime), width(mwidth), split(msplit), wire(mwire), hitIndex(mhitIndex), pnum(mpnum)
        { }
      };

      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      struct lightBuffer{
        double PEs;
        double time;
        double width;
        int opChannel;
        int opHitIndex;
        lightBuffer()
          :PEs(0),time(0),width(0), split(0),opChannel(0), opDet(0), opHitIndex(0), pnum(0)
        {}
        lightBuffer(double mpes, double mtime, double mwidth, double msplit, int mopchan,  int mopindex, int mpnum)
          :PEs(mpes), time(mtime), width(mwidth), split(msplit), opChannel(mopchan),  opHitIndex(mopindex), pnum(mpnum)
        {}
        update(double mpes, double mtime, double mwidth, double msplit, int mopchan,  int mopindex, int mpnum)
          :PEs(mpes), time(mtime), width(mwidth), split(msplit), opChannel(mopchan),  opHitIndex(mopindex), pnum(mpnum)
        {}
      };

      explicit CalibrationTree(fhicl::ParameterSet const& parameterSet);
      explicit CalibrationTree(fhiclConfig const& config);

      virtual void beginJob();
      virtual void analyse(const art::Event& evt) override;
      virtual void endJob();

    private:

      std::vector<jointHitRecord> jhrV;

      void FillJHRV(int trackId, art::Ptr<recob::Hit> hit, double split=1.0);
      void FillJHRV(int trackId, art::Ptr<recob::OpHit> hit, double split=1.0);

      art::ServiceHandle<art::TFileService> private_service_tfs;
      TTree* private_CalibrationTree;
      TBranch* private_CalibrationParticle;
      TBranch* private_CalibrationCharge;
      TBranch* private_CalibrationLight;

      particleBuffer private_particleBuffer;
      hitBuffer      private_hitBuffer;
      opticalBuffer  private_lightBuffer;

      //Also need buffers

  };// end class CalibrationTree

}

#endif//endif CALIBRATIONTREE_H
