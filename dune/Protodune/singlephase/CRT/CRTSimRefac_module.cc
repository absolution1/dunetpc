////////////////////////////////////////////////////////////////////////
// Class:       CRTSimRefac
// Plugin Type: producer (art v2_10_03)
// File:        CRTSimRefac_module.cc
//
// Generated at Wed Jun 27 04:09:39 2018 by Andrew Olivier using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

//     An ART module to simulate how the ProtoDUNE-SP Cosmic Ray Tagger (CRT) system 
//responds to energy deposits. 
//Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
//LArSoft includes

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nug4/ParticleNavigation/ParticleList.h"

//local includes
//#include "CRTTrigger.h"
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//c++ includes
#include <memory>
#include <algorithm>
#include <string>

namespace CRT {
  class CRTSimRefac;
}

class CRT::CRTSimRefac : public art::EDProducer {
public:
  explicit CRTSimRefac(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTSimRefac(CRTSimRefac const &) = delete;
  CRTSimRefac(CRTSimRefac &&) = delete;
  CRTSimRefac & operator = (CRTSimRefac const &) = delete;
  CRTSimRefac & operator = (CRTSimRefac &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;


private:

  //The formats for time and ADC value I will use throughout this module.  If I want to change it later, I only have to change it in one place.  
  typedef long long int time;
  typedef unsigned short adc_t;

  // Member data
  art::InputTag fSimLabel; 


  double fGeVToADC; //Conversion from GeV in detector to ADC counts output by MAROC2.   

  //Parameterization for algorithm that takes voltage entering MAROC2 and converts it to ADC hits.  If I ultimately decide I want to 
  //study impact of more detailed simulation, factor this out into an algorithm for simulating MAROC2 that probably takes continuous 
  //voltage signal as input anyway.  
  time fIntegrationTime; //The time in ns over which charge is integrated in the CRT electronics before being sent to the DAC comparator
  size_t fReadoutWindowSize; //The time in fIntegrationTimes during which activity in a CRT channel is sent to an ADC if the board has 
                             //already decided to trigger
  size_t fDeadtime; //The dead time in fIntegrationTimes after readout during which no energy deposits are processed by CRT boards.
  adc_t fDACThreshold; //DAC threshold for triggering readout for any CRT strip.  
                       //In GeV for now, but needs to become ADC counts one day.  
                       //Should be replaced by either a lookup in a hardware database or 
                       //some constant value one day.  
};


CRT::CRTSimRefac::CRTSimRefac(fhicl::ParameterSet const & p): EDProducer{p}, fSimLabel(p.get<art::InputTag>("SimLabel")), 
                                                              /*fScintillationYield(p.get<double>("ScintillationYield")), 
                                                              fQuantumEff(p.get<double>("QuantumEff")), 
                                                              fDummyGain(p.get<double>("DummyGain")),*/
                                                              fGeVToADC(p.get<double>("GeVToADC")),
                                                              fIntegrationTime(p.get<time>("IntegrationTime")), 
                                                              fReadoutWindowSize(p.get<size_t>("ReadoutWindowSize")), 
                                                              fDeadtime(p.get<size_t>("Deadtime")),
                                                              fDACThreshold(p.get<adc_t>("DACThreshold"))
{
  produces<std::vector<CRT::Trigger>>();
  produces<art::Assns<simb::MCParticle,CRT::Trigger>>(); 

}


void CRT::CRTSimRefac::produce(art::Event & e)
{


  std::vector<art::Handle<sim::AuxDetHitCollection>> allSims;
  e.getManyByType(allSims);

  //const auto & crtHits =  e.getValidHandle<std::vector<sim::AuxDetHit>>("largeant");
  
  art::ServiceHandle < cheat::ParticleInventoryService > partInventory;
  auto trigCol = std::make_unique<std::vector<CRT::Trigger>>();
  //auto simToTrigger = std::make_unique<art::Assns<sim::AuxDetHit, CRT::Trigger>>();
  ///std::unique_ptr< art::Assns<simb::MCParticle, CRT::Trigger>> partToTrigger( new art::Assns<simb::MCParticle, CRT::Trigger>);

  //Utilities to go along with making Assns
  //art::PtrMaker<sim::AuxDetHit> makeSimPtr(e, crtHits.id());
  art::PtrMaker<CRT::Trigger> makeTrigPtr(e);


  art::ServiceHandle<geo::Geometry> geom;
    //<--std::map<int, std::map<time, std::vector<CRT::Hit>>> crtHits;
    std::map<int, std::map<time, std::vector<std::pair<CRT::Hit, int>>>> map_of_crtModule_to_map_of_time_to_vector_of_pair_crtHit_and_trackId;
    for(auto const& auxHits : allSims){
      for(const auto & eDep: * auxHits)
      {
        const size_t tAvg = eDep.GetEntryT();
        //<--crtHits[(eDep.GetID())/64][tAvg/fIntegrationTime].emplace_back(CRT::Hit((eDep.GetID())%64, eDep.GetEnergyDeposited()*0.001f*fGeVToADC));
        map_of_crtModule_to_map_of_time_to_vector_of_pair_crtHit_and_trackId[(eDep.GetID())/64][tAvg/fIntegrationTime].emplace_back(CRT::Hit((eDep.GetID())%64, eDep.GetEnergyDeposited()*0.001f*fGeVToADC),eDep.GetTrackID());
        MF_LOG_DEBUG("TrueTimes") << "Assigned true hit at time " << tAvg << " to bin " << tAvg/fIntegrationTime << ".\n";
      }
    }


  //For each CRT module
  //<--for(const auto & crtModule: crtHits)
  for(const auto & crtModule_to_map_of_time_to_vector_of_pair_crtHit_and_trackId : map_of_crtModule_to_map_of_time_to_vector_of_pair_crtHit_and_trackId)
  {
    int crtChannel = -1;
    int module =   crtModule_to_map_of_time_to_vector_of_pair_crtHit_and_trackId.first;
    
    std::string modStrng="U"+std::to_string(module+1);
    if (module>15) modStrng="D"+std::to_string(module+1-16);
    if ((module+1)==17) crtChannel=22; 
    if ((module+1)==1) crtChannel=24; 
    for (int i=0; i<32; ++i){
      if (crtChannel==22 || crtChannel==24) break;
      const auto& det = geom->AuxDet(i);
      if(det.Name().find(modStrng) != std::string::npos){
        crtChannel=i; break;
      }
    }
    const auto map_of_time_to_vector_of_pair_crtHit_and_trackId = crtModule_to_map_of_time_to_vector_of_pair_crtHit_and_trackId.second;

    MF_LOG_DEBUG("channels") << "Processing channel " << module << "\n";
    
    std::stringstream ss;
    //<--for(const auto& window: timeToHit)
    for(const auto& time_to_vector_of_pair_crtHit_and_trackId : map_of_time_to_vector_of_pair_crtHit_and_trackId)
    {
      ss << "At " << time_to_vector_of_pair_crtHit_and_trackId.first << " ticks, " << time_to_vector_of_pair_crtHit_and_trackId.second.size() << " hits\n";
    } 

    MF_LOG_DEBUG("timeToHitTrackIds") << "Constructed readout windows for module " << crtChannel << ":\n"
                            << ss.str() << "\n";

   
    MF_LOG_DEBUG("timeToHitTrackIds") << "About to loop over " << map_of_time_to_vector_of_pair_crtHit_and_trackId.size() << " time windows that are " << fIntegrationTime << "ns long.\n";
    for(auto time_to_vector_of_pair_crtHit_and_trackId = map_of_time_to_vector_of_pair_crtHit_and_trackId.begin(); time_to_vector_of_pair_crtHit_and_trackId != map_of_time_to_vector_of_pair_crtHit_and_trackId.end(); ) //++time_to_vector_of_pair_crtHit_and_trackId)
    {
      const auto& vector_of_pair_crtHit_and_trackId = time_to_vector_of_pair_crtHit_and_trackId->second;
      const auto aboveThresh = std::find_if(vector_of_pair_crtHit_and_trackId.begin(), vector_of_pair_crtHit_and_trackId.end(), 
                                            [this](const auto& pair_crtHit_and_trackId) { return pair_crtHit_and_trackId.first.ADC() > fDACThreshold; });

      if(aboveThresh != vector_of_pair_crtHit_and_trackId.end()) //If this is true, then I have found a channel above threshold and readout is triggered.  
      {
        //MF_LOG_DEBUG("timeToHitTrackIds") << "Channel " << aboveThresh.Channel() << " has deposit " << aboveThresh.ADC() << " ADC counts that " << "is above threshold.  Triggering readout at time " << time_to_vector_of_pair_crtHit_and_trackId->first << ".\n";
 
        std::vector<CRT::Hit> hits;
        const time timestamp = time_to_vector_of_pair_crtHit_and_trackId->first; //Set timestamp before window is changed.
        const auto end = map_of_time_to_vector_of_pair_crtHit_and_trackId.upper_bound(timestamp+fReadoutWindowSize);
        std::set<uint32_t> channelBusy; //A std::set contains only unique elements.  This set stores channels that have already been read out in 
                                        //this readout window and so are "busy" and cannot contribute any more hits.  
        for(;time_to_vector_of_pair_crtHit_and_trackId != end; ++time_to_vector_of_pair_crtHit_and_trackId)
        {
          for(const auto& pair_crtHit_and_trackId: time_to_vector_of_pair_crtHit_and_trackId->second)
          {
            const auto channel = pair_crtHit_and_trackId.first.Channel();
                                                                  
            if(channelBusy.insert(channel).second) //If this channel hasn't already contributed to this readout window
            {
              hits.push_back(pair_crtHit_and_trackId.first);
             
              /* 
              int trkId = pair_crtHit_and_trackId.second;
              simb::MCParticle &particle = partInventory->TrackIdToParticle(trkId);
              partCol_->push_back(std::move(p));
              art::Ptr<simb::MCParticle> mcp_ptr = art::Ptr<simb::MCParticle>(pid_,partCol_->size()-1,evt->productGetter(pid_));
              tpassn_->addSingle(mct, mcp_ptr, truthInfo);
              */

              ///partToTrigger->addSingle(pair_crtHit_and_trackId.second, makeTrigPtr(trigCol->size()-1));
              //<--//simToTrigger->addSingle(hitPair.second, makeTrigPtr(trigCol->size()-1)); 
              //simToTrigger->addSingle(pair_crtHit_and_trackId.second, makeTrigPtr(trigCol->size()-1)); 
            }
          }
        }

        MF_LOG_DEBUG("CreateTrigger") << "Creating CRT::Trigger...\n";
        //std::cout<< "Timestamp:"<<timestamp*fIntegrationTime<<std::endl;
        trigCol->emplace_back(crtChannel, timestamp*fIntegrationTime, std::move(hits));
        ///partToTrigger->emplace_back(
        ///std::unique_ptr< art::Assns<simb::MCParticle, CRT::Trigger>> partToTrigger( new art::Assns<simb::MCParticle, CRT::Trigger>);

        const auto oldWindow = time_to_vector_of_pair_crtHit_and_trackId;
        if(time_to_vector_of_pair_crtHit_and_trackId != map_of_time_to_vector_of_pair_crtHit_and_trackId.end())
        {
          time_to_vector_of_pair_crtHit_and_trackId = map_of_time_to_vector_of_pair_crtHit_and_trackId.upper_bound(time_to_vector_of_pair_crtHit_and_trackId->first+fDeadtime); 
          MF_LOG_DEBUG("DeadTime") << "Advanced readout window by " << time_to_vector_of_pair_crtHit_and_trackId->first - oldWindow->first
                                                             << " to simulate dead time.\n";
        }
      } //If there was a channel above threshold
      else ++time_to_vector_of_pair_crtHit_and_trackId; //If discriminators did not fire, continue to next time window.
    } //For each time window
  } //For each CRT module

  //Put Triggers and Assns into the event
  MF_LOG_DEBUG("CreateTrigger") << "Putting " << trigCol->size() << " CRT::Triggers into the event at the end of analyze().\n";
  e.put(std::move(trigCol));
  ///e.put(std::move(partToTrigger));
  //e.put(std::move(simToTrigger));
}



DEFINE_ART_MODULE(CRT::CRTSimRefac)
