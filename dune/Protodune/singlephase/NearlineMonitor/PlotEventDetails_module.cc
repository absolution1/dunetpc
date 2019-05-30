////////////////////////////////////////////////////////////////////////
//
//  PlotEventDetails_module.cc
//
//  Author: Leigh Whitehead (leigh.howard.whitehead@cern.ch)
//
//  Module to provide basic event information for protoDUNE
//  Nearline Monitoring
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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TH1F.h"

namespace nlana {
  class PlotEventDetails;
}

class nlana::PlotEventDetails : public art::EDAnalyzer {
public:
  explicit PlotEventDetails(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PlotEventDetails(PlotEventDetails const &) = delete;
  PlotEventDetails(PlotEventDetails &&) = delete;
  PlotEventDetails & operator = (PlotEventDetails const &) = delete;
  PlotEventDetails & operator = (PlotEventDetails &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;

private:

  std::string fParticleProducerLabel;
  std::string fTrackProducerLabel;
  std::string fShowerProducerLabel;

  float fLongTrackCut;

  // Number of tracks
  TH1F *fNParticlesHist;
  TH1F *fNTracksHist;
  TH1F *fNLongTracksHist;
  TH1F *fNShowersHist;

};


nlana::PlotEventDetails::PlotEventDetails(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)  // ,
 // More initializers here.
{
  reconfigure(pset);
  
}

//--------------------------------------------------------------------
void nlana::PlotEventDetails::beginJob()
{
  // Implementation of optional member function here.
  
  art::ServiceHandle<art::TFileService> tfs;

  fNParticlesHist = tfs->make<TH1F>("NParticles",";Number of particles",50,0,200);
  fNTracksHist = tfs->make<TH1F>("NTracks",";Number of tracks",50,0,200);
  fNLongTracksHist = tfs->make<TH1F>("NLongTracks",";Number of long tracks",50,0,200);
  fNShowersHist = tfs->make<TH1F>("NShowers",";Number of showers",50,0,200);

} // beginJob

//--------------------------------------------------------------------
void nlana::PlotEventDetails::reconfigure(fhicl::ParameterSet const & pset)
{
  // The names of the module that produced the reconstructed objects
  fParticleProducerLabel  = pset.get<std::string>("ParticleProducerLabel");
  fTrackProducerLabel  = pset.get<std::string>("TrackProducerLabel");
  fShowerProducerLabel  = pset.get<std::string>("ShowerProducerLabel");

  // Threshold cut for considering a track as long
  fLongTrackCut = pset.get<float>("LongTrackThreshold");

} // reconfigure


//--------------------------------------------------------------------
void nlana::PlotEventDetails::endJob()
{

} // endJob

//--------------------------------------------------------------------
void nlana::PlotEventDetails::analyze(art::Event const & evt)
{
  
  // Try finding some tracks
  art::ValidHandle< std::vector<recob::PFParticle> > particleHandle = evt.getValidHandle<std::vector<recob::PFParticle> >(fParticleProducerLabel);
  art::ValidHandle< std::vector<recob::Track> >  trackHandle = evt.getValidHandle<std::vector<recob::Track> >(fTrackProducerLabel);
  art::ValidHandle< std::vector<recob::Shower> > showerHandle = evt.getValidHandle<std::vector<recob::Shower> >(fShowerProducerLabel);

  fNParticlesHist->Fill(particleHandle->size());

  fNTracksHist->Fill(trackHandle->size());
  unsigned int nLongTracks = 0;
  for(auto const track : (*trackHandle)){
    if(track.Length() > fLongTrackCut) ++nLongTracks; 
  }
  fNLongTracksHist->Fill(nLongTracks);

  fNShowersHist->Fill(showerHandle->size());

} // analyze


DEFINE_ART_MODULE(nlana::PlotEventDetails)
