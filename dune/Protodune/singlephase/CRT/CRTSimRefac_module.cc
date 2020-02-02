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
  std::vector<art::Handle<sim::AuxDetHitCollection>> allSims;
  e.getManyByType(allSims);

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

   

    for(auto time_to_vector_of_pair_crtHit_and_trackId = crtHitsMappedByTime.begin(); time_to_vector_of_pair_crtHit_and_trackId != crtHitsMappedByTime.end(); ) //++time_to_vector_of_pair_crtHit_and_trackId)
    {
      const auto& vector_of_pair_crtHit_and_trackId = time_to_vector_of_pair_crtHit_and_trackId->second;
      const auto aboveThresh = std::find_if(vector_of_pair_crtHit_and_trackId.begin(), vector_of_pair_crtHit_and_trackId.end(), 
                                            [this](const auto& pair_crtHit_and_trackId) { return pair_crtHit_and_trackId.first.ADC() > fDACThreshold; });

      if(aboveThresh != vector_of_pair_crtHit_and_trackId.end()) //If this is true, then I have found a channel above threshold and readout is triggered.  
      {
        mf::LogDebug("timeToHitTrackIds") << "Channel " << (aboveThresh->first).Channel() << " has deposit " << (aboveThresh->first).ADC() << " ADC counts that " << "is above threshold.  Triggering readout at time " << time_to_vector_of_pair_crtHit_and_trackId->first << ".\n";
 
        std::vector<CRT::Hit> hits;
        const time timestamp = time_to_vector_of_pair_crtHit_and_trackId->first; //Set timestamp before window is changed.
        const auto end = crtHitsMappedByTime.upper_bound(timestamp+fReadoutWindowSize);

	    std::set<int> trkIDCheck; //Backtracking set to get rid of duplicate trkIDs that are above threshold

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


	      if (pair_crtHit_and_trackId.first.ADC()>fDACThreshold){

              int tid = pair_crtHit_and_trackId.second;



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

        mf::LogDebug("CreateTrigger") << "Creating CRT::Trigger...\n";

        trigCol->emplace_back(crtChannel, timestamp*fIntegrationTime, std::move(hits));

        const auto oldWindow = time_to_vector_of_pair_crtHit_and_trackId;
        if(time_to_vector_of_pair_crtHit_and_trackId != crtHitsMappedByTime.end())
        {
          time_to_vector_of_pair_crtHit_and_trackId = crtHitsMappedByTime.upper_bound(time_to_vector_of_pair_crtHit_and_trackId->first+fDeadtime); 
          mf::LogDebug("DeadTime") << "Advanced readout window by " << time_to_vector_of_pair_crtHit_and_trackId->first - oldWindow->first
                                                             << " to simulate dead time.\n";
        }
      } //If there was a channel above threshold
      else ++time_to_vector_of_pair_crtHit_and_trackId; //If discriminators did not fire, continue to next time window.
    } //For each time window
  } //For each CRT module

  // -- Put Triggers and Assns into the event
  mf::LogDebug("CreateTrigger") << "Putting " << trigCol->size() << " CRT::Triggers into the event at the end of analyze().\n";
  e.put(std::move(trigCol));
  e.put(std::move(partToTrigger));

}



DEFINE_ART_MODULE(CRT::CRTSimRefac)
