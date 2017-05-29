////////////////////////////////////////////////////////////////////////
// Class:       Lifetime
// Plugin Type: analyzer (art v2_06_03)
// File:        Lifetime_module.cc
//
// Generated at Mon May 29 09:34:51 2017 by Bruce Baller using cetskelgen
// from cetlib version v2_03_00.
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

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"


namespace nlana {
  class Lifetime;
}


class nlana::Lifetime : public art::EDAnalyzer {
public:
  explicit Lifetime(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Lifetime(Lifetime const &) = delete;
  Lifetime(Lifetime &&) = delete;
  Lifetime & operator = (Lifetime const &) = delete;
  Lifetime & operator = (Lifetime &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  std::string fClusterModuleLabel;

};


nlana::Lifetime::Lifetime(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)  // ,
 // More initializers here.
{
  reconfigure(pset);
}

//--------------------------------------------------------------------
void nlana::Lifetime::beginJob()
{
  // Implementation of optional member function here.
  std::cout<<"beginJob: Make ntuples here\n";
} // beginJob

//--------------------------------------------------------------------
void nlana::Lifetime::reconfigure(fhicl::ParameterSet const & pset)
{
  // Implementation of optional member function here.
  fClusterModuleLabel         = pset.get<std::string>("ClusterModuleLabel");
} // reconfigure

//--------------------------------------------------------------------
void nlana::Lifetime::endJob()
{
  // Implementation of optional member function here.
} // endJob

//--------------------------------------------------------------------
void nlana::Lifetime::analyze(art::Event const & evt)
{
  
  int event  = evt.id().event(); 
  int run    = evt.run();
  int subrun = evt.subRun();
  std::cout<<"Inside analyze "<<run<<" "<<" subrun "<<subrun<<" event "<<event<<"\n";
  
  art::ValidHandle<std::vector<recob::Cluster>> clsVecHandle = evt.getValidHandle<std::vector<recob::Cluster>>(fClusterModuleLabel);
  
  for(unsigned int icl = 0; icl < clsVecHandle->size(); ++icl) {
    art::Ptr<recob::Cluster> cls = art::Ptr<recob::Cluster>(clsVecHandle, icl);
    if(cls->NHits() < 200) continue;
    std::cout<<"cls "<<icl<<" "<<cls->Plane()<<" "<<cls->NHits()<<"\n";
  } // icl

} // analyze


DEFINE_ART_MODULE(nlana::Lifetime)
