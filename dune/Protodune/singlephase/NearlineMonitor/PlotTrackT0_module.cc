////////////////////////////////////////////////////////////////////////
//
//  PlotTrackT0_module.cc
//
//  Author: Leigh Whitehead (leigh.howard.whitehead@cern.ch)
//
//  Make a histogram of the T0 of those tracks that have it measured
//
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

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <fstream>

#include "TH1F.h"

namespace nlana {
  class PlotTrackT0;
}

class nlana::PlotTrackT0 : public art::EDAnalyzer {
public:
  explicit PlotTrackT0(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PlotTrackT0(PlotTrackT0 const &) = delete;
  PlotTrackT0(PlotTrackT0 &&) = delete;
  PlotTrackT0 & operator = (PlotTrackT0 const &) = delete;
  PlotTrackT0 & operator = (PlotTrackT0 &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;

private:

  std::string fTrackProducerLabel;
  
  TH1F *fT0Hist;
  TH1F *fT0HistHiRes;
  TH1F *fCosmicT0Hist;
  TH1F *fCosmicT0HistHiRes;

};


nlana::PlotTrackT0::PlotTrackT0(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)  // ,
 // More initializers here.
{
  reconfigure(pset);
  
}

//--------------------------------------------------------------------
void nlana::PlotTrackT0::beginJob()
{
  // Implementation of optional member function here.
  
  art::ServiceHandle<art::TFileService> tfs;
  fT0Hist = tfs->make<TH1F>("TrackT0",";T0 (us)",100,-4000,4000);
  fT0HistHiRes = tfs->make<TH1F>("TrackT0HiRes",";T0 (us)",1000,-4000,4000);
  fCosmicT0Hist = tfs->make<TH1F>("CosmicTrackT0",";T0 (us)",100,-4000,4000);
  fCosmicT0HistHiRes = tfs->make<TH1F>("CosmicTrackT0HiRes",";T0 (us)",1000,-4000,4000);

} // beginJob

//--------------------------------------------------------------------
void nlana::PlotTrackT0::reconfigure(fhicl::ParameterSet const & pset)
{
  // The name of the module that produced the tracks
  fTrackProducerLabel  = pset.get<std::string>("TrackProducerLabel");
} // reconfigure


//--------------------------------------------------------------------
void nlana::PlotTrackT0::endJob()
{

} // endJob

//--------------------------------------------------------------------
void nlana::PlotTrackT0::analyze(art::Event const & evt)
{
  
  // Try finding some tracks
  art::ValidHandle< std::vector<recob::Track> > trackHandle
          = evt.getValidHandle<std::vector<recob::Track> >(fTrackProducerLabel);

  // Find the associations between tracks and T0
  const art::FindManyP<anab::T0> findTrackT0(trackHandle,evt,fTrackProducerLabel);

  // Also look for cosmic tags so we can make a T0 plot for cosmic tagged events only
  const art::FindManyP<anab::CosmicTag> findCosmicTag(trackHandle,evt,fTrackProducerLabel);

  for ( size_t track_index = 0; track_index != trackHandle->size(); ++track_index )
  {
    const auto thisTrack = (*trackHandle)[track_index];

    // Did this track have an associated T0?
    auto const& t0s = findTrackT0.at(track_index);
    if(t0s.size() != 0){
      fT0Hist->Fill(t0s[0]->Time());
      fT0HistHiRes->Fill(t0s[0]->Time());
      // Did this track also have a cosmic tag?
      auto const& tag = findCosmicTag.at(track_index);
      if(tag.size() != 0){
        fCosmicT0Hist->Fill(t0s[0]->Time());
        fCosmicT0HistHiRes->Fill(t0s[0]->Time());
      }
    }
  }

} // analyze


DEFINE_ART_MODULE(nlana::PlotTrackT0)
