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

#include "larsim/MCCheater/BackTracker.h"


#include "TH1F.h"
#include "TProfile.h"


namespace nlana {
  class SPLifetime;
}


class nlana::SPLifetime : public art::EDAnalyzer {
public:
  explicit SPLifetime(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SPLifetime(SPLifetime const &) = delete;
  SPLifetime(SPLifetime &&) = delete;
  SPLifetime & operator = (SPLifetime const &) = delete;
  SPLifetime & operator = (SPLifetime &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  std::string fClusterModuleLabel;
  double fChiCut;
  std::vector<float> fChgCuts;
  int lastRun;
  int fDebugCluster;
  double bigLifeInv[12];
  double bigLifeInvErr[12];
  double bigLifeInvCnt[12];
  
  TH1F *fLifeInv;
  TH1F *fFracSelHits;
  TH1F *fChiDOF;
  TH1F *fLifeInvCA;
  TH1F *fLifeInvAC;
  
  TProfile *fLifeInv_E;
  TProfile *fLifeInv_Angle;

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
  fLifeInv = tfs->make<TH1F>("LifeInv","LifeInv", 50, 0, 1);
  fFracSelHits = tfs->make<TH1F>("FracSelHits","FracSelHits", 50, 0, 1);
  fChiDOF = tfs->make<TH1F>("ChiDOF","ChiDOF", 80, 0, 80);
  fLifeInvCA = tfs->make<TH1F>("LifeInvCA","LifeInv C to A", 50, 0, 1);
  fLifeInvAC = tfs->make<TH1F>("LifeInvAC","LifeInv A to C", 50, 0, 1);
  
  fLifeInv_E = tfs->make<TProfile>("LifeInv_E","LifeInv_E", 20, 0, 40);
  fLifeInv_Angle = tfs->make<TProfile>("LifeInv_Angle","LifeInv_Angle", 15, -1.5, 1.5);
  
  // initialize an invalid run number
  lastRun = -1;
  for(unsigned short tpc = 0; tpc < 12; ++tpc) {
    bigLifeInv[tpc] = 0;
    bigLifeInvErr[tpc] = 0;
    bigLifeInvCnt[tpc] = 0;
  }

} // beginJob

//--------------------------------------------------------------------
void nlana::SPLifetime::reconfigure(fhicl::ParameterSet const & pset)
{
  // Implementation of optional member function here.
  fClusterModuleLabel  = pset.get<std::string>("ClusterModuleLabel");
  fChgCuts             = pset.get<std::vector<float>>("ChgCuts", {0, 5});
  fChiCut              = pset.get<double>("ChiCut", 3);
  fDebugCluster        = pset.get<int>("DebugCluster");
} // reconfigure

//--------------------------------------------------------------------
void nlana::SPLifetime::endJob()
{
  
  std::cout<<"Run  tpc   lifetime  error count\n";
  for(unsigned short tpc = 0; tpc < 12; ++tpc) {
    if(bigLifeInvCnt[tpc] < 2) continue;
    bigLifeInv[tpc] /= bigLifeInvCnt[tpc];
    bigLifeInvErr[tpc] = bigLifeInvErr[tpc] - bigLifeInvCnt[tpc] * bigLifeInv[tpc] * bigLifeInv[tpc];
    if(bigLifeInvErr[tpc] < 0) continue;
    bigLifeInvErr[tpc] = sqrt(bigLifeInvErr[tpc] / (bigLifeInvCnt[tpc] - 1));
    // convert to error on the mean
    bigLifeInvErr[tpc] /= sqrt(bigLifeInvCnt[tpc]);
    std::cout<<lastRun<<" "<<tpc<<" "<<std::fixed<<std::setprecision(3)<<1/bigLifeInv[tpc]<<" "<<bigLifeInvErr[tpc]<<" "<<(int)bigLifeInvCnt[tpc]<<"\n";
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
  
  art::ServiceHandle<cheat::BackTracker> bt;
  const geo::GeometryCore* geom = lar::providerFrom<geo::Geometry>();  

/*
  for(unsigned short tpc = 0;  tpc < 12; ++tpc) {
    std::cout<<"tpc "<<tpc;
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &thetpc = geom->TPC(tpc, 0);
    thetpc.LocalToWorld(local,world);
    std::cout<<" world "<<(int)world[0]<<" "<<(int)world[1]<<" "<<(int)world[2];
    std::cout<<"\n";
  } // tpc
*/
  
  bool prt = false;
  
  // This corresponds to time bins of size 200 ticks * 0.5 us/tick = 0.1 ms
  constexpr float ticksPerHist = 200;
  
  for(unsigned int icl = 0; icl < clsVecHandle->size(); ++icl) {
    art::Ptr<recob::Cluster> cls = art::Ptr<recob::Cluster>(clsVecHandle, icl);
    // only consider the collection plane
    if(cls->Plane().Plane != 2) continue;
    prt = (fDebugCluster >= 0 && icl == (unsigned int)fDebugCluster);
    float sTick = cls->StartTick();
    float eTick = cls->EndTick();
    if(sTick > eTick) std::swap(sTick, eTick);
    float dTick = eTick - sTick;
    unsigned short nhist = 1 + (unsigned short)(dTick / ticksPerHist);
    if(nhist < 5) continue;
    // Get the hits
    std::vector<art::Ptr<recob::Hit> > clsHits;
    clsHitsFind.get(icl, clsHits);
    if(clsHits.size() < 100) continue;
    if(prt) {
      auto& sht = clsHits[0];
      std::cout<<"Cls "<<icl<<" "<<sht->WireID().TPC<<":"<<sht->WireID().Plane<<":"<<sht->WireID().Wire<<":"<<(int)sht->PeakTime();
      auto& eht = clsHits[clsHits.size()-1];
      std::cout<<"  "<<eht->WireID().TPC<<":"<<eht->WireID().Plane<<":"<<eht->WireID().Wire<<":"<<(int)eht->PeakTime();
      std::cout<<" sTick "<<(int)sTick<<" "<<(int)eTick<<" nhits "<<clsHits.size()<<" nhist "<<nhist<<"\n";
    }
    
    // Find the average and maximum charge of these histograms
    std::vector<double> tck(nhist), ave(nhist), cnt(nhist), err(nhist);
    std::vector<float> minChg(nhist, 0);
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
//        if(prt && nit == 0) std::cout<<"Hit "<<pht->WireID().Wire<<":"<<(int)pht->PeakTime()<<" Mult "<<pht->Multiplicity()<<" Q "<<pht->Integral()<<" ihist "<<ihist<<"\n";
//        if(pht->Multiplicity() > 2) continue;
        if(ihist > nhist - 1) continue;
        float chg = pht->Integral();
        if(chg < minChg[ihist] || chg > maxChg[ihist]) continue;
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
          maxChg[ihist] = fChgCuts[1] * ave[ihist];
          minChg[ihist] = fChgCuts[0] * ave[ihist];
          if(prt) std::cout<<ihist<<" Min "<<(int)minChg[ihist]<<" max "<<(int)maxChg[ihist]<<"\n";
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

    // fit to find the 1/lifetime
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
      xx = (tck[ihist] - tck[0]) * msPerTick;
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
    if(B > 0) continue;
    // calculate the error
    double ndof = fitcnt - 2;
    double varnce = (sumy2 + A*A*sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
    if(varnce == 0) continue;
    double BErr = sqrt(varnce * sum / delta);
    if(prt) std::cout<<"B "<<B<<" "<<BErr;
    
    double lifeInv = -B;
    double lifeInvErr = BErr / (B * B) ;
    
    // calculate chisq
    double chi = 0;
    for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
      if(cnt[ihist] < 3) continue;
      if(err[ihist] == 0) continue;
      xx = (tck[ihist] - tck[0]) * msPerTick;
      yy = exp(A - xx * lifeInv);
      arg = (yy - ave[ihist]) / err[ihist];
      chi += arg * arg;
      if(prt) std::cout<<"chk "<<ihist<<" xx "<<xx<<" yy "<<yy<<" ave "<<ave[ihist]<<" arg "<<arg<<"\n";
    }
    chi /= ndof;
    unsigned short tpc = clsHits[0]->WireID().TPC;
    if(prt) std::cout<<tpc<<" lifeInv "<<lifeInv<<" +/- "<<lifeInvErr<<" chi "<<chi<<"\n";
    fChiDOF->Fill(chi);
    if(chi > fChiCut) continue;
    ++bigLifeInvCnt[tpc];
    bigLifeInv[tpc] += lifeInv;
    bigLifeInvErr[tpc] += lifeInv * lifeInv;
    fLifeInv->Fill(lifeInv);
//    std::cout<<"icl "<<icl<<" tpc "<<tpc<<" lifeInv "<<lifeInv<<" chi "<<chi<<"\n";
    
    double selHits = 0;
    for(auto& hcnt : cnt) selHits += hcnt;
    selHits /= (double)(clsHits.size());
    fFracSelHits->Fill(selHits);

    // temp study MC truth
    art::Ptr<recob::Hit> pht = clsHits[0];
    raw::ChannelID_t channel = geom->PlaneWireToChannel((int)pht->WireID().Plane, (int)pht->WireID().Wire, (int)pht->WireID().TPC, (int)pht->WireID().Cryostat);
    double startTick = pht->PeakTime() - pht->RMS();
    double endTick = pht->PeakTime() + pht->RMS();

    std::vector<sim::TrackIDE> tides;
    bt->ChannelToTrackIDEs(tides, channel, startTick, endTick);
    if(tides.size() != 1) continue;
    int trackID = tides[0].trackID;
    auto part = bt->TrackIDToParticle(trackID);
    if(abs(part->PdgCode()) != 13) std::cout<<"Not a muon "<<part->PdgCode()<<"\n";
//    std::cout<<"TPC "<<tpc<<" start E "<<part->E()<<" Pos "<<(int)part->Vx()<<" "<<(int)part->Vy()<<" "<<(int)part->Vz();
//    std::cout<<" end "<<part->EndE()<<" Pos "<<(int)part->EndX()<<" "<<(int)part->EndY()<<" "<<(int)part->EndZ();
    if(part->Vy() < part->EndY()) std::cout<<" CR going UP!! \n";
    bool goingPosX = (part->Vx() < part->EndX());
    bool posTPC = (tpc == 2 || tpc == 6 || tpc ==10);
    bool CtoA = goingPosX && posTPC;
    if(CtoA) {
//      std::cout<<" C -> A\n";
      fLifeInvCA->Fill(lifeInv);
    } else {
//      std::cout<<" A -> C\n";
      fLifeInvAC->Fill(lifeInv);
    }
    fLifeInv_E->Fill(part->E(), lifeInv);
    fLifeInv_Angle->Fill(cls->StartAngle(), lifeInv);
  } // icl

} // analyze


DEFINE_ART_MODULE(nlana::SPLifetime)
