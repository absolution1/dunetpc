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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TH1F.h"


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
  
  TH1F *fdTick;

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
  
  art::ServiceHandle<art::TFileService> tfs;
  fdTick = tfs->make<TH1F>("dTick","dTick", 430, 0, 4300);

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
  
  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double msPerTick = 1E-6 * detprop->SamplingRate();
  std::cout<<"Inside analyze "<<run<<" "<<" subrun "<<subrun<<" event "<<event<<" msPerTick "<<msPerTick<<"\n";

  art::ValidHandle<std::vector<recob::Cluster>> clsVecHandle = evt.getValidHandle<std::vector<recob::Cluster>>(fClusterModuleLabel);
  art::FindManyP<recob::Hit> clsHitsFind(clsVecHandle, evt, fClusterModuleLabel);
  
  // This should be about 10 cm
  constexpr float ticksPerHist = 200;
  
  double bigLife = 0;
  double bigLifeErr = 0;
  double bigLifeCnt = 0;
  for(unsigned int icl = 0; icl < clsVecHandle->size(); ++icl) {
    art::Ptr<recob::Cluster> cls = art::Ptr<recob::Cluster>(clsVecHandle, icl);
    // only consider the collection plane
    if(cls->Plane().Plane != 2) continue;
    float sTick = cls->StartTick();
    float eTick = cls->EndTick();
    if(sTick > eTick) std::swap(sTick, eTick);
    float dTick = eTick - sTick;
    fdTick->Fill(dTick);
    unsigned short nhist = dTick / ticksPerHist;
    if(nhist < 5) continue;
    // Get the hits
    std::vector<art::Ptr<recob::Hit> > clsHits;
    clsHitsFind.get(icl, clsHits);
    if(clsHits.size() < 100) continue;
    // Find the average and maximum charge of these histograms
    std::vector<double> ave(nhist), cnt(nhist), err(nhist);
    std::vector<float> maxChg(nhist, 10000);
    for(unsigned short nit = 0; nit < 2; ++nit) {
      for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
        ave[ihist] = 0;
        cnt[ihist] = 0;
        err[ihist] = 0;
      } // ihist
      // Sum to get the (truncated) average
      for(auto& pht : clsHits) {
        unsigned short ihist = (pht->PeakTime() - sTick) / ticksPerHist;
        // require single hits
        if(pht->Multiplicity() != 1) continue;
        if(ihist > nhist - 1) continue;
        float chg = pht->Integral();
        if(chg > maxChg[ihist]) continue;
        ave[ihist] += chg;
        err[ihist] += chg * chg;
        ++cnt[ihist];
      } // ii
      for(unsigned short ihist = 0; ihist < nhist; ++ihist) if(cnt[ihist] > 0) ave[ihist] /= cnt[ihist];
      if(nit == 0) {
        // Calculate a max charge cut on the first iteration
        for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
          maxChg[ihist] = 0;
          if(cnt[ihist] < 5) continue;
          maxChg[ihist] = 1.2 * ave[ihist];
        } // ihist
      } else {
        // Calculate the error on the average on the second iteration
        for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
          if(cnt[ihist] < 3) continue;
          double arg = err[ihist] - cnt[ihist] * ave[ihist] * ave[ihist];
          if(arg < 0) {
            cnt[ihist] = 0;
            continue;
          }
          // calculate rms
          err[ihist] = sqrt(arg / (cnt[ihist] - 1));
          // convert to error on the mean
          err[ihist] /= sqrt(cnt[ihist]);
        } // ihist
      } // second iteration
    } // nit
    
    // require a minimum count in each histogram
    unsigned short nok = 0;
    for(auto& icnt : cnt) if(icnt > 4) ++nok;
    if(nok != cnt.size()) continue;
    
    auto& sht = clsHits[0];
    std::cout<<"Cls "<<sht->WireID().TPC<<":"<<sht->WireID().Plane<<":"<<sht->WireID().Wire<<":"<<(int)sht->PeakTime();
    auto& eht = clsHits[clsHits.size()-1];
    std::cout<<"  "<<eht->WireID().TPC<<":"<<eht->WireID().Plane<<":"<<eht->WireID().Wire<<":"<<(int)eht->PeakTime();
    std::cout<<" sTick "<<(int)sTick<<" "<<(int)eTick<<"\n";
    for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
      float xx = (ihist + 0.5) * ticksPerHist + sTick;
      std::cout<<"ihist "<<ihist<<" tick "<<(int)xx<<" ave "<<(int)ave[ihist]<<" +/- "<<(int)err[ihist]<<" cnt "<<(int)cnt[ihist]<<"\n";
    } // ihist

    // fit to find the lifetime
    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;
    double xx, yy, wght, arg;
    
    for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
      if(err[ihist] == 0) continue;
      xx = (ihist + 0.5) * ticksPerHist;
      yy = log(ave[ihist]);
      // error on log(x) = dx / x
      arg = ave[ihist] / err[ihist];
      wght = arg * arg;
      sum += wght;
      sumx += wght * xx;
      sumy += wght * yy;
      sumx2 += wght * xx * xx;
      sumxy += wght * xx * yy;
      sumy2 += wght * yy * yy;
    } // ihist
    // calculate coefficients and std dev
    double delta = sum * sumx2 - sumx * sumx;
    if(delta == 0.) continue;
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    // slope = 1 / lifetime
    double B = (sumxy * sum  - sumx * sumy) / delta;
    if(B > 0) {
      std::cout<<"Positive lifetime "<<B<<"\n";
      continue;
    }
    B = - 1 / B;
//    std::cout<<"B "<<B;
    double life = B * msPerTick;
    
    // calculate chisq
    double chi = 0;
    for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
      xx = (ihist + 0.5) * ticksPerHist;
      yy = exp(A - xx / B);
      arg = (yy - ave[ihist]) / err[ihist];
      chi += arg * arg;
//      std::cout<<"chk "<<ihist<<" xx "<<xx<<" yy "<<yy<<" ave "<<ave[ihist]<<" arg "<<arg<<"\n";
    }
    chi /= (double)(nhist - 2);
    std::cout<<"life "<<life<<" chi "<<chi<<"\n";
    if(chi > 5) continue;
    ++bigLifeCnt;
    bigLife += life;
    bigLifeErr += life * life;

  } // icl
  
  if(bigLifeCnt < 2) return;
  bigLife /= bigLifeCnt;
  bigLifeErr = bigLifeErr - bigLifeCnt * bigLife * bigLife;
  if(bigLifeErr < 0) return;
  bigLifeErr = sqrt(bigLifeErr / (bigLifeCnt - 1));
  // convert to error on the mean
  bigLifeErr /= sqrt(bigLifeCnt);
  std::cout<<"bigLife "<<bigLife<<" +/- "<<bigLifeErr<<" bigLifeCnt "<<bigLifeCnt<<"\n";

} // analyze


DEFINE_ART_MODULE(nlana::Lifetime)
