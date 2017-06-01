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

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
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
  float fFiducialCut;
  float fTickLo[12] {0};
  float fTickHi[12] {0};

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
  fFiducialCut                = pset.get<float>("FiducialCut");
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

  

  static bool first = true;
  if(first) {
    first = false;
    // Get the low and high tick range for plane 2 in each TPC
    const geo::GeometryCore* geom = lar::providerFrom<geo::Geometry>();
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    
    for (const geo::TPCID& tpcid: geom->IterateTPCIDs()) {
      geo::TPCGeo const& tpc = geom->TPC(tpcid);
      unsigned short not03 = (tpcid.TPC % 4);
      if(not03 == 0 || not03 == 3) continue;
      tpc.LocalToWorld(local,world);
      double xx = world[0]-geom->DetHalfWidth(tpcid.TPC, tpcid.Cryostat) + fFiducialCut;
      fTickLo[tpcid.TPC] = detprop->ConvertXToTicks(xx, 2, tpcid.TPC, tpcid.Cryostat);
      xx = world[0]+geom->DetHalfWidth(tpcid.TPC, tpcid.Cryostat) - fFiducialCut;
      fTickHi[tpcid.TPC] = detprop->ConvertXToTicks(xx, 2, tpcid.TPC, tpcid.Cryostat);
      if(fTickLo[tpcid.TPC] > fTickHi[tpcid.TPC]) std::swap(fTickLo[tpcid.TPC], fTickHi[tpcid.TPC]);
      std::cout<<"TPC "<<tpcid<<" Lo "<<fTickLo[tpcid.TPC]<<" Hi "<<fTickHi[tpcid.TPC]<<"\n";
    } // tpcid
  }
  
  const geo::GeometryCore* geom = lar::providerFrom<geo::Geometry>();
  for (const geo::TPCID& tpcid: geom->IterateTPCIDs()) {
    unsigned short not03 = (tpcid.TPC % 4);
    if(not03 == 0 || not03 == 3) continue;
    std::cout<<"tpc "<<tpcid<<" Lo "<<fTickLo[tpcid.TPC]<<" Hi "<<fTickHi[tpcid.TPC]<<"\n";
  }

  art::ValidHandle<std::vector<recob::Cluster>> clsVecHandle = evt.getValidHandle<std::vector<recob::Cluster>>(fClusterModuleLabel);
  
  for(unsigned int icl = 0; icl < clsVecHandle->size(); ++icl) {
    art::Ptr<recob::Cluster> cls = art::Ptr<recob::Cluster>(clsVecHandle, icl);
    // only consider the collection plane
    if(cls->Plane().Plane != 2) continue;
    float sTick = cls->StartTick();
    float eTick = cls->EndTick();
    if(sTick > eTick) std::swap(sTick, eTick);
    if(cls->StartTick() > fTickLo[cls->Plane().TPC] || cls->EndTick() < fTickHi[cls->Plane().TPC]) continue;
    std::cout<<"cls "<<icl<<" "<<cls->Plane().TPC<<" "<<(int)cls->StartTick()<<" "<<(int)cls->EndTick()<<"\n";
  } // icl

} // analyze


DEFINE_ART_MODULE(nlana::Lifetime)
