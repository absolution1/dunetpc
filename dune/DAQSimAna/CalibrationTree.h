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

      explicit CalibrationTree(fhicl::ParameterSet const& parameterSet);
      explicit CalibrationTree(fhiclConfig& config);

      virtual void beginJob();
      virtual void analyse(const art::Event& evt) override;
      virtual void endJob();

    private:
      
      TTree* private_CalibrationTree;
      TBranch* private_CalibrationParticle;
      TBranch* private_CalibrationCharge;
      TBranch* private_CalibrationLight;

      //Also need buffers

  }// end class CalibrationTree

}

#endif//endif CALIBRATIONTREE_H
