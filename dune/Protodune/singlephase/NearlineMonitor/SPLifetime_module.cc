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
#include "art/Framework/Core/FileBlock.h"
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
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <array>

//#include "larsim/MCCheater/BackTracker.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"

#define setHistTitles(hist,xtitle,ytitle) hist->GetXaxis()->SetTitle(xtitle); hist->GetYaxis()->SetTitle(ytitle);

//// apa = tpcMapping[tpc]
const std::array<size_t,13> tpcMapping = {{0,4,1,0,0,5,2,0,0,6,3,0,0}};
//// hist->GetXaxis()->SetBinLabel(iBin+1,apaLabels[iBin])
const std::array<std::string,6> apaLabels = {{
                                                "APA-DaS-US/APA5",
                                                "APA-DaS-MS/APA6",
                                                "APA-DaS-DS/APA4",
                                                "APA-RaS-US/APA3", 
                                                "APA-RaS-MS/APA2",
                                                "APA-RaS-DS/APA1",
                                            }};

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
  void reconfigure(fhicl::ParameterSet const & p);
  void respondToOpenInputFile(art::FileBlock const & infileblock) override;

private:

  std::string fClusterModuleLabel;
  double fChiCut;
  std::vector<float> fChgCuts;
  int lastRun;
  int fDebugCluster;
  std::string fInFilename;
  bool fIsRealData;
  double bigLifeInv[12];
  double bigLifeInvErr[12];
  double bigLifeInvCnt[12];
  double signalToNoise[12];
  double signalToNoiseCnt[12];
  unsigned int signalToNoiseClsCnt[12];
  
  TH1F *fLifeInv;
  TH1F *fFracSelHits;
  TH1F *fChiDOF;
  TH1F *fLifeInvCA;
  TH1F *fLifeInvAC;
  
  TProfile *fLifeInv_E;
  TProfile *fLifeInv_Angle;

  TH1F *fDriftTime;
  TH2F *fDriftTimeVTPC;

  TH1F *fSNR;
  TH2F *fSNRVTPC;

  TH1F *fAmplitudes;
  TH1F *fNoise;

};


nlana::SPLifetime::SPLifetime(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset),
  fInFilename("NoInFilenameFound"),
  fIsRealData(true)

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
    signalToNoise[tpc] = 0;
    signalToNoiseCnt[tpc] = 0;
    signalToNoiseClsCnt[tpc] = 0;
  }

  fDriftTime = tfs->make<TH1F>("DriftTime","Drift Time", 500, 0.,10.);
  setHistTitles(fDriftTime,"Cluster Drift Time [ms]", "Clusters / Bin");
  fDriftTimeVTPC = tfs->make<TH2F>("DriftTimeVTPC","Drift Time v. APA", 6,0.5,6.5,500,0.,10.);
  setHistTitles(fDriftTimeVTPC,"","Cluster Drift Time [ms]");

  fSNR = tfs->make<TH1F>("SNR","Signal to Noise Ratio", 1000, 0.,10000.);
  setHistTitles(fSNR,"Signal to Noise Ratio", "Hits / Bin");
  fSNRVTPC = tfs->make<TH2F>("SNRVTPC","Signal to Noise Ratio v. APA", 6,0.5,6.5,1000,0.,10000.);
  setHistTitles(fSNRVTPC,"TPC Number","Signal to Noise Ratio");

  fAmplitudes = tfs->make<TH1F>("Amplitudes","Hit Amplitude", 2000, 0.,4000.);
  setHistTitles(fAmplitudes,"Hit Amplitude [ADC]", "Hits / Bin");
  fNoise = tfs->make<TH1F>("Noise","Wire Noise", 2000, 0.,5.);
  setHistTitles(fNoise,"RMS Noise [ADC]", "Hit Wires / Bin");

  for(size_t apa = 0; apa < 6; ++apa){
    fDriftTimeVTPC->GetXaxis()->SetBinLabel(apa+1,apaLabels.at(apa).c_str());
    fSNRVTPC->GetXaxis()->SetBinLabel(apa+1,apaLabels.at(apa).c_str());
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
  std::ofstream purfile;
  purfile.open("Lifetime_Run" + std::to_string(lastRun) + ".txt");
  //purfile<<"Run, tpc, lifetime, error, count, S/N, num S/N clusters\n";
  purfile<<"Run, tpc, lifetime, error, count, S/N, num S/N clusters, drift time\n";

  for(unsigned short tpc = 0; tpc < 12; ++tpc) {
    if(bigLifeInvCnt[tpc] < 2) continue;
    if(bigLifeInv[tpc] <= 0) continue;
    bigLifeInv[tpc] /= bigLifeInvCnt[tpc];
    bigLifeInvErr[tpc] = bigLifeInvErr[tpc] - bigLifeInvCnt[tpc] * bigLifeInv[tpc] * bigLifeInv[tpc];
    if(bigLifeInvErr[tpc] < 0) continue;
    bigLifeInvErr[tpc] = sqrt(bigLifeInvErr[tpc] / (bigLifeInvCnt[tpc] - 1));
    // convert to error on the mean
    bigLifeInvErr[tpc] /= sqrt(bigLifeInvCnt[tpc]);
    float life = 0;
    if(1/bigLifeInv[tpc] != 0) life = 1 / bigLifeInv[tpc];
    float lifeErr = life * bigLifeInvErr[tpc] / bigLifeInv[tpc];
    float sn = -1;
    if(signalToNoiseCnt[tpc] > 100) sn = signalToNoise[tpc] / signalToNoiseCnt[tpc];
    purfile<<lastRun<<", "<<tpc<<", "<<std::fixed<<std::setprecision(2)<<life<<", "<<lifeErr<<", "<<(int)bigLifeInvCnt[tpc];
    purfile<<", "<<std::setprecision(1)<<sn<<", "<<signalToNoiseClsCnt[tpc];

    // Now do drift time
    size_t apa = tpcMapping.at(tpc);
    if (apa == 0) continue;
    size_t nBinsY = fDriftTimeVTPC->GetNbinsY();
    float driftTime = -1.;
    for (int iBin=nBinsY; iBin >=0; iBin--)
    {
      if (fDriftTimeVTPC->GetBinContent(apa,iBin) > 1)
      {
        driftTime = fDriftTimeVTPC->GetYaxis()->GetBinCenter(iBin);
        break;
      }
    }
    purfile<<", "<<std::setprecision(2)<<driftTime;

    purfile<<"\n";
  }
  purfile.close();

  // Now for images
  TCanvas * canvas = new TCanvas("canvas_SPLifetime");
  canvas->SetLogz();
  canvas->SetBottomMargin(0.13);
  canvas->SetRightMargin(0.12);
  std::string imageFileName;
  // Try to get rid of directory and .root extension
  std::string infilenameStripped = fInFilename;
  if (fIsRealData)
  {
    size_t slashPos = infilenameStripped.find_last_of("/");
    if (slashPos != std::string::npos)
    {
      infilenameStripped = infilenameStripped.substr(slashPos+1,std::string::npos); // get rid of directory
    }
    infilenameStripped = infilenameStripped.substr(0, infilenameStripped.find_last_of(".")); // get rid of .root
  }
  else
  {
    // "run": "run003907_0001_dl05",
    std::stringstream infilenamefake;
    infilenamefake << "run";
    infilenamefake << std::setfill('0') << std::setw(6) << lastRun;
    infilenamefake << std::setw(0) << '_';
    infilenamefake << std::setw(4) << 1;
    infilenamefake << std::setw(0) << "_dl";
    infilenamefake << std::setw(2) << 5;
    infilenameStripped = infilenamefake.str();
  }

  // summary json file
  std::ofstream summaryfile;
  summaryfile.open("summary_purity.json");
  summaryfile << "[\n  {\n    \"run\": \"" << infilenameStripped << "\",\n"
              << "    \"Type\": \"purity\"\n  }\n]";
  summaryfile.close();

  // file list json file
  std::ofstream filelistfile;
  filelistfile.open("purity_FileList.json");
  filelistfile << "[\n  {\n    \"Category\": \"Purity Monitor\",\n    \"Files\": {\n      \"Cluster Drift Time\": \"";

  imageFileName = "driftVTPC_";
  imageFileName += infilenameStripped;
  imageFileName += ".png";
  fDriftTimeVTPC->SetStats(false);
  fDriftTimeVTPC->GetXaxis()->SetLabelSize(0.050);
  fDriftTimeVTPC->Draw("colz");
  canvas->SaveAs(imageFileName.c_str());
  filelistfile << imageFileName<<",";

  imageFileName = "driftVTPC_zoom_";
  imageFileName += infilenameStripped;
  imageFileName += ".png";
  std::string originalTitle = fDriftTimeVTPC->GetTitle();
  fDriftTimeVTPC->SetTitle((originalTitle+" From 0 to 4 ms").c_str());
  fDriftTimeVTPC->GetYaxis()->SetRangeUser(0,4);
  fDriftTimeVTPC->Draw("colz");
  canvas->SaveAs(imageFileName.c_str());
  fDriftTimeVTPC->SetTitle(originalTitle.c_str());
  filelistfile << imageFileName;

  filelistfile << "\"\n    }\n  }\n]";
  filelistfile.close();
  delete canvas;

} // endJob

//--------------------------------------------------------------------
void nlana::SPLifetime::respondToOpenInputFile(art::FileBlock const & infileblock)
{
    fInFilename = infileblock.fileName();
} // respondToOpenInputFile

//--------------------------------------------------------------------
void nlana::SPLifetime::analyze(art::Event const & evt)
{
  
//  int event  = evt.id().event(); 
  int run    = evt.run();
//  int subrun = evt.subRun();
  fIsRealData = evt.isRealData();
  
  if(lastRun < 0) lastRun = run;
  
  if(run != lastRun) mf::LogVerbatim("LIFE")<<"The run number has changed from "<<lastRun<<" to "<<run<<" This might not be good...";
  
  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double msPerTick = 1E-6 * detprop->SamplingRate();


  art::ValidHandle<std::vector<recob::Cluster>> clsVecHandle = evt.getValidHandle<std::vector<recob::Cluster>>(fClusterModuleLabel);
  art::FindManyP<recob::Hit> clsHitsFind(clsVecHandle, evt, fClusterModuleLabel);
  art::ValidHandle<std::vector<recob::Wire>> wireVecHandle = evt.getValidHandle<std::vector<recob::Wire>>("caldata");
  
//  bool prt = false;
  
  // This corresponds to time bins of size 200 ticks * 0.5 us/tick = 0.1 ms
  constexpr float ticksPerHist = 200;
  
  for(unsigned int icl = 0; icl < clsVecHandle->size(); ++icl) {
    art::Ptr<recob::Cluster> cls = art::Ptr<recob::Cluster>(clsVecHandle, icl);
    // only consider the collection plane
    if(cls->Plane().Plane != 2) continue;
//    prt = (fDebugCluster >= 0 && icl == (unsigned int)fDebugCluster);
    float sTick = cls->StartTick();
    float eTick = cls->EndTick();
    if(sTick > eTick) std::swap(sTick, eTick);
    float dTick = eTick - sTick;
    // Get the hits
    std::vector<art::Ptr<recob::Hit> > clsHits;
    clsHitsFind.get(icl, clsHits);
    if(clsHits.size() == 0) continue;
    unsigned short tpc = clsHits[0]->WireID().TPC;
    fDriftTime->Fill(dTick*msPerTick);
    fDriftTimeVTPC->Fill(tpcMapping.at(tpc),dTick*msPerTick);
    unsigned short nhist = 1 + (unsigned short)(dTick / ticksPerHist);
    if(nhist < 5) continue;
    if(clsHits.size() < 100) continue;
/*
    if(prt) {
      auto& sht = clsHits[0];
      std::cout<<"Cls "<<icl<<" "<<sht->WireID().TPC<<":"<<sht->WireID().Plane<<":"<<sht->WireID().Wire<<":"<<(int)sht->PeakTime();
      auto& eht = clsHits[clsHits.size()-1];
      std::cout<<"  "<<eht->WireID().TPC<<":"<<eht->WireID().Plane<<":"<<eht->WireID().Wire<<":"<<(int)eht->PeakTime();
      std::cout<<" sTick "<<(int)sTick<<" "<<(int)eTick<<" nhits "<<clsHits.size()<<" nhist "<<nhist<<"\n";
    }
*/
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
//          if(prt) std::cout<<ihist<<" Min "<<(int)minChg[ihist]<<" max "<<(int)maxChg[ihist]<<"\n";
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
//      if(prt) std::cout<<" failed nok cut "<<nok<<" Need at least 4 \n";
      continue;
    }
/*
    if(prt) {
      for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
        float xx = (ihist + 0.5) * ticksPerHist + sTick;
        std::cout<<"ihist "<<ihist<<" tick "<<(int)xx<<" ave "<<(int)ave[ihist]<<" +/- "<<(int)err[ihist]<<" cnt "<<(int)cnt[ihist]<<"\n";
      } // ihist
    }
*/
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
//    double BErr = sqrt(varnce * sum / delta);
//    if(prt) std::cout<<"B "<<B<<" "<<BErr;
    
    double lifeInv = -B;
//    double lifeInvErr = BErr / (B * B) ;
    
    // calculate chisq
    double chi = 0;
    for(unsigned short ihist = 0; ihist < nhist; ++ihist) {
      if(cnt[ihist] < 3) continue;
      if(err[ihist] == 0) continue;
      xx = (tck[ihist] - tck[0]) * msPerTick;
      yy = exp(A - xx * lifeInv);
      arg = (yy - ave[ihist]) / err[ihist];
      chi += arg * arg;
//      if(prt) std::cout<<"chk "<<ihist<<" xx "<<xx<<" yy "<<yy<<" ave "<<ave[ihist]<<" arg "<<arg<<"\n";
    }
    chi /= ndof;
//    if(prt) std::cout<<tpc<<" lifeInv "<<lifeInv<<" +/- "<<lifeInvErr<<" chi "<<chi<<"\n";
    fChiDOF->Fill(chi);
    if(chi > fChiCut) continue;
    ++bigLifeInvCnt[tpc];
    bigLifeInv[tpc] += lifeInv;
    bigLifeInvErr[tpc] += lifeInv * lifeInv;
    fLifeInv->Fill(lifeInv);
    
    double selHits = 0;
    for(auto& hcnt : cnt) selHits += hcnt;
    selHits /= (double)(clsHits.size());
    fFracSelHits->Fill(selHits);

    fLifeInv_Angle->Fill(cls->StartAngle(), lifeInv);
  } // icl
  
  std::vector<std::pair<unsigned int, float>> chanPedRMS;
//  std::cout << "n wire data: "<<wireVecHandle->size() << std::endl;
  for(unsigned int witer = 0; witer < wireVecHandle->size(); ++witer) {
    art::Ptr<recob::Wire> thisWire(wireVecHandle, witer);
    const recob::Wire::RegionsOfInterest_t& signalROI = thisWire->SignalROI();
    for(const auto& range : signalROI.get_ranges()) {
      const std::vector<float>& signal = range.data();
//      std::cout << "channel: "<<thisWire->Channel()<<" range begin: " << range.begin_index() << " range size: "<<range.size()<<"signal length: " << signal.size() << std::endl;
      if(signal.size() != 1) continue;
      // protect against unrealistic values
      float rms = signal[0];
      if(rms < 0.001) rms = 0.001;
      chanPedRMS.push_back(std::make_pair(thisWire->Channel(), rms));
//      std::cout << "channel: "<<thisWire->Channel()<<" rms: " << rms << std::endl;
    }
  }// witer

  // calculate S/N using shallow angle clusters
  for(unsigned int icl = 0; icl < clsVecHandle->size(); ++icl) {
    art::Ptr<recob::Cluster> cls = art::Ptr<recob::Cluster>(clsVecHandle, icl);
    // only consider the collection plane
    if(cls->Plane().Plane != 2) continue;
    float dTick = std::abs(cls->StartTick() - cls->EndTick());
    if(dTick < 500) continue;
    float dWire = std::abs(cls->StartWire() - cls->EndWire());
    if(dWire < 300) continue;
    // Get the hits
    std::vector<art::Ptr<recob::Hit> > clsHits;
    clsHitsFind.get(icl, clsHits);
    if(clsHits.size() < 300) continue;
    unsigned short tpc = clsHits[0]->WireID().TPC;
    ++signalToNoiseClsCnt[tpc];
//    float aveSN = 0;
//    float cnt = 0;
    for(auto& hit : clsHits) {
      float pedRMS = -1;
      for(auto& chrms : chanPedRMS) {
        if(chrms.first != hit->Channel()) continue;
        pedRMS = chrms.second;
        break;
      } // chrms
      if(pedRMS < 0) continue;
      float snr = hit->PeakAmplitude() / pedRMS;
      //std::cout<<"SNR "<<(int)snr << " S: " << hit->PeakAmplitude() << " N: "<< pedRMS <<"\n";
      signalToNoise[tpc] += snr;
      ++signalToNoiseCnt[tpc];
      fAmplitudes->Fill(hit->PeakAmplitude());
      fNoise->Fill(pedRMS);
      fSNR->Fill(snr);
      fSNRVTPC->Fill(tpcMapping.at(tpc),snr);
//      aveSN += sn;
//      ++cnt;
    } // hit
/*
    if(cnt == 0) continue;
    aveSN /= cnt;
    std::cout<<" aveSN "<<aveSN<<" cnt "<<(int)cnt<<"\n";
*/
  } // icl

} // analyze

DEFINE_ART_MODULE(nlana::SPLifetime)
