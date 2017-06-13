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
  class SPLifetime;
}


class nlana::SPLifetime : public art::EDAnalyzer {
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
  float fChiCut;
  float fLandauCut;
  int lastRun;
  double bigLife[12];
  double bigLifeErr[12];
  double bigLifeCnt[12];
  
  TH1F *fLife;
  TH1F *fFracSelHits;
  TH1F *fChiDOF;
  TH1F *fChgHist[11];

};


nlana::SPLifetime::SPLifetime(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)  // ,
 // More initializers here.
{
  reconfigure(pset);
  
}

//--------------------------------------------------------------------
void nlana::SPLifetime::beginJob()
{
  // Implementation of optional member function here.
  
  art::ServiceHandle<art::TFileService> tfs;
  fLife = tfs->make<TH1F>("Lifetime","Lifetime", 100, 0, 10);
  fFracSelHits = tfs->make<TH1F>("FracSelHits","FracSelHits", 50, 0, 1);
  fChiDOF = tfs->make<TH1F>("ChiDOF","ChiDOF", 80, 0, 800);
  
  unsigned short nbins = 40;
  float maxbin = 800;
  fChgHist[0] = tfs->make<TH1F>("Chg0","Chg0", nbins, 0, maxbin);
  fChgHist[1] = tfs->make<TH1F>("Chg1","Chg0", nbins, 0, maxbin);
  fChgHist[2] = tfs->make<TH1F>("Chg2","Chg0", nbins, 0, maxbin);
  fChgHist[3] = tfs->make<TH1F>("Chg3","Chg0", nbins, 0, maxbin);
  fChgHist[4] = tfs->make<TH1F>("Chg4","Chg0", nbins, 0, maxbin);
  fChgHist[5] = tfs->make<TH1F>("Chg5","Chg0", nbins, 0, maxbin);
  fChgHist[6] = tfs->make<TH1F>("Chg6","Chg0", nbins, 0, maxbin);
  fChgHist[7] = tfs->make<TH1F>("Chg7","Chg0", nbins, 0, maxbin);
  fChgHist[8] = tfs->make<TH1F>("Chg8","Chg0", nbins, 0, maxbin);
  fChgHist[9] = tfs->make<TH1F>("Chg9","Chg0", nbins, 0, maxbin);
  fChgHist[10] = tfs->make<TH1F>("Chg10","Chg0", nbins, 0, maxbin);
  
  // initialize an invalid run number
  lastRun = -1;
  for(unsigned short tpc = 0; tpc < 12; ++tpc) {
    bigLife[tpc] = 0;
    bigLifeErr[tpc] = 0;
    bigLifeCnt[tpc] = 0;
  }

} // beginJob

//--------------------------------------------------------------------
void nlana::SPLifetime::reconfigure(fhicl::ParameterSet const & pset)
{
  // Implementation of optional member function here.
  fClusterModuleLabel  = pset.get<std::string>("ClusterModuleLabel");
  fLandauCut           = pset.get<float>("LandauCut");
  fChiCut              = pset.get<float>("ChiCut");
} // reconfigure

//--------------------------------------------------------------------
void nlana::SPLifetime::endJob()
{
  
  std::cout<<"Run  tpc   lifetime  error count\n";
  for(unsigned short tpc = 0; tpc < 12; ++tpc) {
    if(bigLifeCnt[tpc] < 2) continue;
    bigLife[tpc] /= bigLifeCnt[tpc];
    bigLifeErr[tpc] = bigLifeErr[tpc] - bigLifeCnt[tpc] * bigLife[tpc] * bigLife[tpc];
    if(bigLifeErr[tpc] < 0) continue;
    bigLifeErr[tpc] = sqrt(bigLifeErr[tpc] / (bigLifeCnt[tpc] - 1));
    // convert to error on the mean
    bigLifeErr[tpc] /= sqrt(bigLifeCnt[tpc]);
    std::cout<<lastRun<<" "<<tpc<<" "<<std::fixed<<std::setprecision(3)<<bigLife[tpc]<<" "<<bigLifeErr[tpc]<<" "<<(int)bigLifeCnt[tpc]<<"\n";
  }
} // endJob

//--------------------------------------------------------------------
void nlana::SPLifetime::analyze(art::Event const & evt)
{
  
//  int event  = evt.id().event(); 
  int run    = evt.run();
//  int subrun = evt.subRun();
  
  if(lastRun < 0) lastRun = run;
  
  if(run != lastRun) mf::LogVerbatim("LIFE")<<"The run number has changed from "<<lastRun<<" to "<<run<<" This might not be good...";
  
  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double msPerTick = 1E-6 * detprop->SamplingRate();
//  std::cout<<"Inside analyze "<<run<<" "<<" subrun "<<subrun<<" event "<<event<<" msPerTick "<<msPerTick<<"\n";

  art::ValidHandle<std::vector<recob::Cluster>> clsVecHandle = evt.getValidHandle<std::vector<recob::Cluster>>(fClusterModuleLabel);
  art::FindManyP<recob::Hit> clsHitsFind(clsVecHandle, evt, fClusterModuleLabel);
  
  bool prt = true;
  
  // This corresponds to time bins of size 200 ticks * 0.5 us/tick = 0.1 ms
  constexpr float ticksPerHist = 200;
  
  for(unsigned int icl = 0; icl < clsVecHandle->size(); ++icl) {
    art::Ptr<recob::Cluster> cls = art::Ptr<recob::Cluster>(clsVecHandle, icl);
    prt = (icl == 4);
    // only consider the collection plane
    if(cls->Plane().Plane != 2) continue;
    float sTick = cls->StartTick();
    float eTick = cls->EndTick();
    if(sTick > eTick) std::swap(sTick, eTick);
    float dTick = eTick - sTick;
    unsigned short nhist = 1 + dTick / ticksPerHist;
    if(nhist < 5) continue;
    // Get the hits
    std::vector<art::Ptr<recob::Hit> > clsHits;
    clsHitsFind.get(icl, clsHits);
    if(clsHits.size() < 100) continue;
//    if(prt) {
      auto& sht = clsHits[0];
      std::cout<<"Cls "<<icl<<" "<<sht->WireID().TPC<<":"<<sht->WireID().Plane<<":"<<sht->WireID().Wire<<":"<<(int)sht->PeakTime();
      auto& eht = clsHits[clsHits.size()-1];
      std::cout<<"  "<<eht->WireID().TPC<<":"<<eht->WireID().Plane<<":"<<eht->WireID().Wire<<":"<<(int)eht->PeakTime();
      std::cout<<" sTick "<<(int)sTick<<" "<<(int)eTick<<" nhits "<<clsHits.size()<<"\n";
//    }
    // Find the average and maximum charge of these histograms
    std::vector<double> tck(nhist), ave(nhist), cnt(nhist), err(nhist);
    std::vector<float> maxChg(nhist, 10000);
    for(unsigned short nit = 0; nit < 2; ++nit) {
      for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
        tck[ihist] = 0;
        ave[ihist] = 0;
        cnt[ihist] = 0;
        err[ihist] = 0;
      } // ihist
      // Sum to get the (truncated) average
      for(auto& pht : clsHits) {
        unsigned short ihist = (pht->PeakTime() - sTick) / ticksPerHist;
        // require single hits
        if(prt && nit == 0) std::cout<<"Hit "<<pht->WireID().Wire<<":"<<(int)pht->PeakTime()<<" Mult "<<pht->Multiplicity()<<" Q "<<pht->Integral()<<" ihist "<<ihist<<"\n";
        if(pht->Multiplicity() > 2) continue;
        if(ihist > nhist - 1) continue;
        float chg = pht->Integral();
        if(chg > maxChg[ihist]) continue;
        if(prt && nit == 0) fChgHist[ihist]->Fill(chg);
        tck[ihist] += pht->PeakTime() - sTick;
        ave[ihist] += chg;
        err[ihist] += chg * chg;
        ++cnt[ihist];
      } // ii
      for(unsigned short ihist = 0; ihist < nhist; ++ihist) if(cnt[ihist] > 0) {
        tck[ihist] /= cnt[ihist];
        ave[ihist] /= cnt[ihist];
      }
      if(nit == 0) {
        // Calculate a max charge cut on the first iteration
        for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
          maxChg[ihist] = 0;
          if(cnt[ihist] < 5) continue;
          maxChg[ihist] = fLandauCut * ave[ihist];
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
    if(nok < 5) {
      if(prt) std::cout<<" failed nok cut "<<nok<<" Need at least 4 \n";
      continue;
    }

    if(prt) {
      for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
        float xx = (ihist + 0.5) * ticksPerHist + sTick;
        std::cout<<"ihist "<<ihist<<" tick "<<(int)xx<<" ave "<<(int)ave[ihist]<<" +/- "<<(int)err[ihist]<<" cnt "<<(int)cnt[ihist]<<"\n";
      } // ihist
    }

    // fit to find the lifetime
    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;
    double xx, yy, wght, arg;
    double fitcnt = 0;
    
    for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
      if(cnt[ihist] < 3) continue;
      if(err[ihist] == 0) continue;
      xx = tck[ihist];
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
      ++fitcnt;
    } // ihist
    // calculate coefficients
    double delta = sum * sumx2 - sumx * sumx;
    if(delta == 0.) continue;
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    // slope = 1 / lifetime
    double B = (sumxy * sum  - sumx * sumy) / delta;
    if(B > 0) {
//      std::cout<<"Positive lifetime "<<B<<"\n";
      continue;
    }
    B = - 1 / B;
//    std::cout<<"B "<<B;
    double life = B * msPerTick;
    
    // calculate chisq
    double chi = 0;
    for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
      if(cnt[ihist] < 3) continue;
      if(err[ihist] == 0) continue;
      xx = tck[ihist];
      yy = exp(A - xx / B);
      arg = (yy - ave[ihist]) / err[ihist];
      chi += arg * arg;
      if(prt) std::cout<<"chk "<<ihist<<" xx "<<xx<<" yy "<<yy<<" ave "<<ave[ihist]<<" arg "<<arg<<"\n";
    }
    chi /= (double)(fitcnt - 2);
    unsigned short tpc = clsHits[0]->WireID().TPC;
    if(prt) std::cout<<tpc<<" life "<<life<<" chi "<<chi<<"\n";
    fChiDOF->Fill(chi);
    if(chi > fChiCut) continue;
    ++bigLifeCnt[tpc];
    bigLife[tpc] += life;
    bigLifeErr[tpc] += life * life;
    fLife->Fill(life);
    
    double selHits = 0;
    for(auto& hcnt : cnt) selHits += hcnt;
    selHits /= (double)(clsHits.size());
    fFracSelHits->Fill(selHits);
  } // icl

} // analyze


DEFINE_ART_MODULE(nlana::SPLifetime)
