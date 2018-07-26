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

//c++ includes
#include <numeric> //std::accumulate was moved from <algorithm> to <numeric> in c++14

namespace CRT {
  class CRTSimValidation;
}

//Take advantage of "backwards associations" between CRT::Trigger and 
//simb::AuxDetSimChannel to make sure that CRTSim_module is behaving in 
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
  TH1D* fRMSWithinTrigger; //RMS of associated IDEs' times w.r.t. Trigger's timestamp.  
                           //Will show me that readout window is working.  
  TH1D* fTriggerDeltaT; //Times between Triggers on a given module.  There should be no 
                        //entries below simulated dead time.
};


CRT::CRTSimValidation::CRTSimValidation(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), fCRTLabel(p.get<art::InputTag>("CRTLabel"))
{
  consumes<std::vector<CRT::Trigger>>();
  consumes<std::vector<art::Assn<simb::AuxDetSimChannel, CRT::Trigger>>>();
}

void CRT::CRTSimValidation::analyze(art::Event const & e)
{
  //Get the data products I plan to work with
  const auto& triggers = e.getValidHandle<std::vector<CRT::Trigger>>(fCRTLabel);

  art::FindManyP<simb::AuxDetSimChannel> trigToSim(triggers, e, fCRTLabel);

  //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
  art::ServiceHandle<geo::Geometry> geom;

  //Mapping from channel to previous Trigger time
  std::unordered_map<size_t, double> prevTimes;
  for(const auto& trigger: triggers)
  {
    fDetsWithHits->Fill(geom->AuxDet(trigger.Channel()).GetName().c_str(), 1.0);
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
    const auto& trig = triggers[trigIt];
    const auto& simChans = trigToSim.at(trigIt);

    size_t nIDEs = 0;
    const auto sumSq = std::accumulate(simChans.begin(), simChans.end(), 
                                       [&trig, &nIDEs](double sum, const auto& channel)
                                       {
                                         const auto& ides = channel.IDEs();
                                         nIDEs += ides.size();
                                         return sum + std::accumulate(ides.begin(), ides.end(), [&trig](double subSum, const auto& ide)
                                                                                                {
                                                                                                  const auto diff = trig.Timestamp() - 
                                                                                                                    (ide.StopT - ide.StartT);
                                                                                                  return sumSum + diff*diff;
                                                                                                });
                                       });
    fRMSWithinTrigger->Fill(std::sqrt(sumSq/nIDEs)); //Root of the mean of squares
  }
}

void CRT::CRTSimValidation::beginJob()
{
  //Register ROOT objects I want to write with the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  fDetsWithHits     = tfs->make<TH1D>("DetsWithHits", "AuxDetGeo Names that Had CRT::Triggers;AuxDetGeo Name;Events", 1, 0, -1); 
  //Trigger automatic binning since am plotting strings.
  fAllADCs          = tfs->make<TH1D>("AllADCs", "ADC Values from CRT::Hits on All Channels;ADC Value;Hits", 500, 0, 2000);
  fRMSWithinTrigger = tfs->make<TH1D>("RMSWithinTrigger", "RMS of True Times Associated with a CRT::Trigger;True Time RMS;Energy Deposits", 
                                      20, 0, 30);
  fTriggerDeltaT    = tfs->make<TH1D>("TriggerDeltaT", "#DeltaT Between Consecutive Triggers in a Module;#DeltaT;Trigger Pairs", 20, 0, 30);
}

DEFINE_ART_MODULE(CRT::CRTSimValidation)
