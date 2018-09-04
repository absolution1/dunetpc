////////////////////////////////////////////////////////////////////////
// Class:       CRTSimValidation
// Plugin Type: analyzer (art v2_11_02)
// File:        CRTSimValidation_module.cc
//
// Generated at Thu Jul 26 11:49:48 2018 by Andrew Olivier using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

//LArSoft includes
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"

//Local includes
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//ROOT includes
#include "TH1D.h"

//c++ includes
#include <numeric> //std::accumulate was moved from <algorithm> to <numeric> in c++14

namespace CRT {
  class CRTSimValidation;
}

//Take advantage of "backwards associations" between CRT::Trigger and 
//sim::AuxDetSimChannel to make sure that CRTSim_module is behaving in 
//a believable way.  Produces plots by analyzing CRT::Triggers and the 
//AuxDetSimChannels associated with them (thus only working with MC files).

class CRT::CRTSimValidation : public art::EDAnalyzer {
public:
  explicit CRTSimValidation(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTSimValidation(CRTSimValidation const &) = delete;
  CRTSimValidation(CRTSimValidation &&) = delete;
  CRTSimValidation & operator = (CRTSimValidation const &) = delete;
  CRTSimValidation & operator = (CRTSimValidation &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  //Parameters for reading in CRT::Triggers and associated AuxDetSimChannels.
  art::InputTag fCRTLabel; //Label for the module that produced CRT::Triggers I 
                           //plan to analyze.

  //Histograms I plan to write.  Some TFileService-managed resource owns them, 
  //so all data members here are observer pointers that will not be deleted.  
  TH1D* fDetsWithHits; //Histogram names of AuxDetGeos that have CRT::Triggers
  TH1D* fAllADCs; //Histogram all CRT::Hit ADC values in a single plot.  In 
                  //particle gun, this will show me that DAC threshold parameter 
                  //is having an effect.  With background simulation, could compare with 
                  //data to see how realistic background is.
  TH1D* fDevWithinTrigger; //Similar to sample standard deviation of associated IDEs' times taking
                           //Trigger's timestamp as mean.  Will show me that readout window is working.  
  TH1D* fDeltaTAssoc; //Time difference between each IDE associated with a Trigger and that Trigger's timestamp
  TH1D* fTriggerDeltaT; //Times between Triggers on a given module.  There should be no 
                        //entries below simulated dead time.
};


CRT::CRTSimValidation::CRTSimValidation(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), fCRTLabel(p.get<art::InputTag>("CRTLabel"))
{
  consumes<std::vector<CRT::Trigger>>(fCRTLabel);
  consumes<std::vector<art::Assns<sim::AuxDetSimChannel, CRT::Trigger>>>(fCRTLabel);
}

void CRT::CRTSimValidation::analyze(art::Event const & e)
{
  //Get the data products I plan to work with
  const auto& triggers = e.getValidHandle<std::vector<CRT::Trigger>>(fCRTLabel);

  art::FindManyP<sim::AuxDetSimChannel> trigToSim(triggers, e, fCRTLabel);

  //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
  art::ServiceHandle<geo::Geometry> geom;

  //Mapping from channel to previous Trigger time
  std::unordered_map<size_t, double> prevTimes;
  for(const auto& trigger: *triggers)
  {
    fDetsWithHits->Fill(geom->AuxDet(trigger.Channel()).Name().c_str(), 1.0);
    const auto& hits = trigger.Hits();
    for(const auto& hit: hits) fAllADCs->Fill(hit.ADC());

    const auto found = prevTimes.find(trigger.Channel());
    if(found != prevTimes.end())
    {
      fTriggerDeltaT->Fill(trigger.Timestamp() - found->second);
      found->second = trigger.Timestamp();
    } 
    else prevTimes[trigger.Channel()] = trigger.Timestamp();
  }

  for(size_t trigIt = 0; trigIt < trigToSim.size(); ++trigIt)
  {
    const auto& trig = (*triggers)[trigIt];
    const auto& simChans = trigToSim.at(trigIt);

    size_t nIDEs = 0;
    const auto sumSq = std::accumulate(simChans.begin(), simChans.end(), 0., 
                                       [this, &trig, &nIDEs](double sum, const auto& channel)
                                       {
                                         const auto& ides = channel->AuxDetIDEs();
                                         nIDEs += ides.size();
                                         return sum + std::accumulate(ides.begin(), ides.end(), 0., 
                                                                      [this, &trig](double subSum, const auto& ide)
                                                                      {
                                                                        const auto diff = trig.Timestamp() - 
                                                                                          (ide.exitT + ide.entryT)/2.;
                                                                        fDeltaTAssoc->Fill(diff);
                                                                        if(fabs(diff) > 200) //TODO: Read pset from original module to use readout 
                                                                                             //      time here
                                                                        {
                                                                          mf::LogInfo("LargeDeltaT") << "Found large deltaT: " << diff << "\n"
                                                                                                     << "Timestamp is " << trig.Timestamp() << "\n"
                                                                                                     << "True time is " << (ide.exitT + 
                                                                                                                            ide.entryT)/2. << ".\n";
                                                                        }

                                                                        return subSum + diff*diff;
                                                                      });
                                       });
    fDevWithinTrigger->Fill(std::sqrt(sumSq/(nIDEs-1))); //Root of the mean of squares
  }
}

void CRT::CRTSimValidation::beginJob()
{
  //Register ROOT objects I want to write with the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  fDetsWithHits     = tfs->make<TH1D>("DetsWithHits", "AuxDetGeo Names that Had CRT::Triggers;AuxDetGeo Name;Events", 1, 0, -1); 
  //Trigger automatic binning since am plotting strings.
  fAllADCs          = tfs->make<TH1D>("AllADCs", "ADC Values from CRT::Hits on All Channels;ADC Value;Hits", 500, 0, 2000);
  fDevWithinTrigger = tfs->make<TH1D>("DevWithinTrigger", "Deviation of True Times from Associated Trigger Timestamp;"
                                                          "True Time #sigma;Energy Deposits", 100, 0, 300);
  fDeltaTAssoc      = tfs->make<TH1D>("DeltaTAssoc", "Time Difference Between Associated IDEs and Trigger Timestamp;#DeltaT_{IDE};Energy Deposits", 
                                      200, -100, 100);
  fTriggerDeltaT    = tfs->make<TH1D>("TriggerDeltaT", "#DeltaT Between Consecutive Triggers in a Module;#DeltaT;Trigger Pairs", 20, 0, 30);
}

DEFINE_ART_MODULE(CRT::CRTSimValidation)
