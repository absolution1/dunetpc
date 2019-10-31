////////////////////////////////////////////////////////////////////////
// Class:       CheckCRT
// Plugin Type: analyzer (art v3_00_00)
// File:        CheckCRT_module.cc
//
// Generated at Mon Jan  7 12:31:13 2019 by Tingjun Yang using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TTree.h"

namespace pdsp {
  class CheckCRT;
}


class pdsp::CheckCRT : public art::EDAnalyzer {
public:
  explicit CheckCRT(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CheckCRT(CheckCRT const&) = delete;
  CheckCRT(CheckCRT&&) = delete;
  CheckCRT& operator=(CheckCRT const&) = delete;
  CheckCRT& operator=(CheckCRT&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;
  std::vector<uint32_t> crtmask;
  std::vector<double> trkstartx_pandora;
  std::vector<double> trkstarty_pandora;
  std::vector<double> trkstartz_pandora;
  std::vector<double> trkendx_pandora;
  std::vector<double> trkendy_pandora;
  std::vector<double> trkendz_pandora;
  std::vector<double> t0_pandora;

  std::vector<double> trkstartx_pmtrack;
  std::vector<double> trkstarty_pmtrack;
  std::vector<double> trkstartz_pmtrack;
  std::vector<double> trkendx_pmtrack;
  std::vector<double> trkendy_pmtrack;
  std::vector<double> trkendz_pmtrack;
  std::vector<double> t0_pmtrack;

};


pdsp::CheckCRT::CheckCRT(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
}

void pdsp::CheckCRT::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("crt","crt tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("crtmask", &crtmask);
  fTree->Branch("trkstartx_pandora",&trkstartx_pandora);
  fTree->Branch("trkstarty_pandora",&trkstarty_pandora);
  fTree->Branch("trkstartz_pandora",&trkstartz_pandora);
  fTree->Branch("trkendx_pandora",&trkendx_pandora);
  fTree->Branch("trkendy_pandora",&trkendy_pandora);
  fTree->Branch("trkendz_pandora",&trkendz_pandora);
  fTree->Branch("t0_pandora",&t0_pandora);

  fTree->Branch("trkstartx_pmtrack",&trkstartx_pmtrack);
  fTree->Branch("trkstarty_pmtrack",&trkstarty_pmtrack);
  fTree->Branch("trkstartz_pmtrack",&trkstartz_pmtrack);
  fTree->Branch("trkendx_pmtrack",&trkendx_pmtrack);
  fTree->Branch("trkendy_pmtrack",&trkendy_pmtrack);
  fTree->Branch("trkendz_pmtrack",&trkendz_pmtrack);
  fTree->Branch("t0_pmtrack",&t0_pmtrack);

}

void pdsp::CheckCRT::analyze(art::Event const& e)
{

  run = e.run();
  subrun = e.subRun();
  event = e.id().event();

  crtmask.clear();

  trkstartx_pandora.clear();
  trkstarty_pandora.clear();
  trkstartz_pandora.clear();
  trkendx_pandora.clear();
  trkendy_pandora.clear();
  trkendz_pandora.clear();
  t0_pandora.clear();

  trkstartx_pmtrack.clear();
  trkstarty_pmtrack.clear();
  trkstartz_pmtrack.clear();
  trkendx_pmtrack.clear();
  trkendy_pmtrack.clear();
  trkendz_pmtrack.clear();
  t0_pmtrack.clear();

  art::ValidHandle<std::vector<raw::RDTimeStamp>> timeStamps = e.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");

  int trigger = -1;
  // Check that we have good information
  if(timeStamps.isValid() && timeStamps->size() == 1){
    // Access the trigger information. Beam trigger flag = 0xc
    const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
    trigger = timeStamp.GetFlags();
  }

  if (trigger!=13) return;

  auto const& pdspctbs = *e.getValidHandle<std::vector<raw::ctb::pdspctb>>("ctbrawdecoder:daq");
  if (!pdspctbs.empty()){
    const size_t npdspctbs = pdspctbs.size();
    for (size_t i=0; i<npdspctbs; ++i){
//      std::cout << "Number of triggers: " << pdspctbs[i].GetNTriggers() << std::endl;
//      std::cout << "Number of chstatuses: " << pdspctbs[i].GetNChStatuses() << std::endl;
//      std::cout << "Number of feedbacks: " << pdspctbs[i].GetNFeedbacks() << std::endl;
//      std::cout << "Number of miscs: " << pdspctbs[i].GetNMiscs() << std::endl;
//      
//      const std::vector<raw::ctb::Trigger> HLTriggers = pdspctbs[i].GetHLTriggers();
//      std::cout << "Number of HL Triggers: " << HLTriggers.size() << std::endl;
//      for (size_t j=0; j<HLTriggers.size(); ++j){
//        std::cout << "HL Trigger: " << j << " " << HLTriggers[j].word_type << " " << HLTriggers[j].trigger_word << " " << HLTriggers[j].timestamp << std::endl;
//      }
      
      const std::vector<raw::ctb::ChStatus> chs = pdspctbs[i].GetChStatusAfterHLTs();
      std::cout << "Number of CH Status after HLTs: " << chs.size() << std::endl;
      for (size_t j=0; j<chs.size(); ++j){
        //std::cout << " chs after HLT: " << j << " " << chs[j].word_type << " " << chs[j].pds << " " << chs[j].crt << " " << chs[j].beam_hi << " " << chs[j].beam_lo << std::endl;
        crtmask.push_back(chs[j].crt);
      }
//      
//      
//      const std::vector<raw::ctb::Trigger> LLTriggers = pdspctbs[i].GetLLTriggers();
//      std::cout << "Number of LL Triggers: " << LLTriggers.size() << std::endl;
//      for (size_t j=0; j<LLTriggers.size(); ++j){
//        std::cout << "LL Trigger: " << j << " " << LLTriggers[j].word_type << " " << LLTriggers[j].trigger_word << " " << LLTriggers[j].timestamp << std::endl;
//      }		
//      
//      const std::vector<raw::ctb::WordIndex> indexes = pdspctbs[i].GetIndexes();
//      std::cout << "Number of Indexes: " << indexes.size() << std::endl;
//      for (size_t j=0; j<indexes.size(); ++j){
//        std::cout << "Index: " << j << " " << indexes[j].word_type << " " << indexes[j].index << std::endl;
//      }		
      
    }
  }

  art::Handle< std::vector<recob::Track> > pandoratrkHandle;
  std::vector< art::Ptr<recob::Track> > pandoratrks;
  if (e.getByLabel("pandoraTrack", pandoratrkHandle))
    art::fill_ptr_vector(pandoratrks, pandoratrkHandle);
  
  art::Handle< std::vector<recob::PFParticle> > pfpListHandle;
  e.getByLabel("pandora", pfpListHandle);
  art::FindManyP<recob::PFParticle> fmpfp(pandoratrkHandle, e, "pandoraTrack");
  art::FindManyP<anab::T0> fmt0pandora(pfpListHandle, e, "pandora");

  for (size_t i = 0; i<pandoratrks.size(); ++i){
    auto & trk = pandoratrks[i];
    trkstartx_pandora.push_back(trk->Vertex().X());
    trkstarty_pandora.push_back(trk->Vertex().Y());
    trkstartz_pandora.push_back(trk->Vertex().Z());
    trkendx_pandora.push_back(trk->End().X());
    trkendy_pandora.push_back(trk->End().Y());
    trkendz_pandora.push_back(trk->End().Z());
    double t0 = 0;
    auto &pfps = fmpfp.at(trk.key());
    if (!pfps.empty()){
      auto &t0s = fmt0pandora.at(pfps[0].key());
      if (!t0s.empty()){
        t0 = t0s[0]->Time();
      }
    }
    t0_pandora.push_back(t0);
  }

  /*
  art::Handle< std::vector<recob::Track> > pmtracktrkHandle;
  std::vector< art::Ptr<recob::Track> > pmtracktrks;
  if (e.getByLabel("pmtrack", pmtracktrkHandle))
    art::fill_ptr_vector(pmtracktrks, pmtracktrkHandle);
  
  art::FindManyP<anab::T0> fmt0pmtrack(pmtracktrkHandle, e, "pmtrack");

  for (size_t i = 0; i<pmtracktrks.size(); ++i){
    auto & trk = pmtracktrks[i];
    trkstartx_pmtrack.push_back(trk->Vertex().X());
    trkstarty_pmtrack.push_back(trk->Vertex().Y());
    trkstartz_pmtrack.push_back(trk->Vertex().Z());
    trkendx_pmtrack.push_back(trk->End().X());
    trkendy_pmtrack.push_back(trk->End().Y());
    trkendz_pmtrack.push_back(trk->End().Z());
    double t0 = 0;
    auto &t0s = fmt0pmtrack.at(trk.key());
    if (!t0s.empty()){
      t0 = t0s[0]->Time();
    }
    t0_pmtrack.push_back(t0);
  }
  */
  fTree->Fill();
}

DEFINE_ART_MODULE(pdsp::CheckCRT)
