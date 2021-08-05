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
#include "dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//c++ includes
#include <memory>
#include <algorithm>
#include <string>
#include <map>
#include <unordered_map>

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

  // -- message logger
  mf::LogInfo logInfo_;

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


CRT::CRTSimRefac::CRTSimRefac(fhicl::ParameterSet const & p): EDProducer{p}, 
                                                              logInfo_("CRTSimRefactor"),
                                                              fSimLabel(p.get<art::InputTag>("SimLabel")), 
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


  //const auto & crtHits =  e.getValidHandle<std::vector<sim::AuxDetHit>>("largeant");
  auto const allSims = e.getMany<sim::AuxDetHitCollection>();

  // -- Get all MCParticles to do assns later
  const auto & mcp_handle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant"); // -- TODO: make this an input tag
  art::PtrMaker<simb::MCParticle> makeMCParticlePtr{e,mcp_handle.id()};
  art::ServiceHandle < cheat::ParticleInventoryService > partInventory;
  auto const & mcparticles = *(mcp_handle); //dereference the handle

  // -- Construct map of trackId to MCParticle handle index to do assns later
  std::unordered_map<int, int> map_trackID_to_handle_index;
  for (size_t idx = 0; idx < mcparticles.size(); ++idx){
    int tid = mcparticles[idx].TrackId();
    map_trackID_to_handle_index.insert(std::make_pair(tid,idx));
  }
  
  auto trigCol = std::make_unique<std::vector<CRT::Trigger>>();

  std::unique_ptr< art::Assns<simb::MCParticle, CRT::Trigger>> partToTrigger( new art::Assns<simb::MCParticle, CRT::Trigger>);

  art::PtrMaker<CRT::Trigger> makeTrigPtr(e);


  art::ServiceHandle<geo::Geometry> geom;
  std::map<int, std::map<time, std::vector<std::pair<CRT::Hit, int>>>> crtHitsModuleMap;
  for(auto const& auxHits : allSims){
    for(const auto & eDep: * auxHits)
    {
      const size_t tAvg = eDep.GetEntryT();
      crtHitsModuleMap[(eDep.GetID())/64][tAvg/fIntegrationTime].emplace_back(CRT::Hit((eDep.GetID())%64, eDep.GetEnergyDeposited()*0.001f*fGeVToADC),eDep.GetTrackID());
      mf::LogDebug("TrueTimes") << "Assigned true hit at time " << tAvg << " to bin " << tAvg/fIntegrationTime << ".\n";
    }
  }


  // -- For each CRT module
  for(const auto & crtHitsMappedByModule : crtHitsModuleMap)
  {
    int crtChannel = -1;
    int module = crtHitsMappedByModule.first;
    
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
    const auto crtHitsMappedByTime = crtHitsMappedByModule.second;

    mf::LogDebug("channels") << "Processing channel " << module << "\n";
    
    std::stringstream ss;

    mf::LogDebug("timeToHitTrackIds") << "Constructed readout windows for module " << crtChannel << ":\n"
                            << ss.str() << "\n";



     auto lastTimeStamp=time(0);
    int i=0;
  for(auto window : crtHitsMappedByTime) 
    {
	if (i!=0 && (time(fDeadtime)+lastTimeStamp)>window.first && lastTimeStamp<window.first) continue;
	i++;
      const auto& hitsInWindow = window.second;
      const auto aboveThresh = std::find_if(hitsInWindow.begin(), hitsInWindow.end(), 
                                            [this](const auto& hitPair) { return hitPair.first.ADC() > fDACThreshold; });


if(aboveThresh != hitsInWindow.end()){

       std::vector<CRT::Hit> hits;
	std::set<int> trkIDCheck;
        const time timestamp = window.first; //Set timestamp before window is changed.
        const time end = (timestamp+fReadoutWindowSize);
        std::set<uint32_t> channelBusy; //A std::set contains only unique elements.  This set stores channels that have already been read out in 
                                        //this readout window and so are "busy" and cannot contribute any more hits. 
    for(auto busyCheckWindow : crtHitsMappedByTime ){
	if (time(busyCheckWindow.first)<timestamp || time(busyCheckWindow.first)>end) continue;
          for(const auto& hitPair: window.second){
            const auto channel = hitPair.first.Channel(); //TODO: Get channel number back without needing art::Ptr here.  
                                                                    //      Maybe store crt::Hits.
            if(channelBusy.insert(channel).second){ 
              hits.push_back(hitPair.first);


	      if (hitPair.first.ADC()>fDACThreshold){

              int tid = hitPair.second;



	      trkIDCheck.insert(tid);

	  
	    }
            }
          }
       }

	for (int tid : trkIDCheck){
	      // -- safe index retrieval
	      int index = 0;
              auto search = map_trackID_to_handle_index.find(tid);
              if (search != map_trackID_to_handle_index.end()){
                index = search->second;
                mf::LogDebug("GetAssns") << "Found index : " << index;
              } else {                                                                                        
                mf::LogDebug("GetAssns") << "No matching index... strange";
                continue;
              }
              
              // -- Sanity check, not needed
              simb::MCParticle particle = mcparticles[index];
              mf::LogDebug("GetAssns") << "TrackId from particle obtained with index " << index
                << " is : " << particle.TrackId() << " , expected: " << tid;
              mf::LogDebug("GetMCParticle") << particle;
              
              auto const mcptr = makeMCParticlePtr(index);
              partToTrigger->addSingle(mcptr, makeTrigPtr(trigCol->size()-1));

	}
	//std::cout<<"Hits Generated:"<<hits.size()<<std::endl;
        lastTimeStamp=window.first;

        MF_LOG_DEBUG("CreateTrigger") << "Creating CRT::Trigger...\n";
        trigCol->emplace_back(crtChannel, timestamp*fIntegrationTime, std::move(hits)); 
	} // For each readout with a triggerable hit
    }  // For each time window

  } //For each CRT module

  // -- Put Triggers and Assns into the event
  mf::LogDebug("CreateTrigger") << "Putting " << trigCol->size() << " CRT::Triggers into the event at the end of analyze().\n";
  e.put(std::move(trigCol));
  e.put(std::move(partToTrigger));

}



DEFINE_ART_MODULE(CRT::CRTSimRefac)
