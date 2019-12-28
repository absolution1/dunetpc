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

//LArSoft includes

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/AuxDetHit.h"

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

}


void CRT::CRTSimRefac::produce(art::Event & e)
{


  std::vector<art::Handle<sim::AuxDetHitCollection>> allSims;
  e.getManyByType(allSims);


  //const auto & crtHits =  e.getValidHandle<std::vector<sim::AuxDetHit>>("largeant");
  

  auto trigCol = std::make_unique<std::vector<CRT::Trigger>>();
  //auto simToTrigger = std::make_unique<art::Assns<sim::AuxDetHit, CRT::Trigger>>();

  //Utilities to go along with making Assns
  //art::PtrMaker<sim::AuxDetHit> makeSimPtr(e, crtHits.id());
  art::PtrMaker<CRT::Trigger> makeTrigPtr(e);


  art::ServiceHandle<geo::Geometry> geom;
    std::map<int, std::map<time, std::vector<CRT::Hit>>> crtHits;
    for(auto const& auxHits : allSims){
    for(const auto & eDep: * auxHits)
    {
	
        const size_t tAvg = (eDep.GetExitT()+eDep.GetEntryT())/2.;
        crtHits[(eDep.GetID())/64][tAvg/fIntegrationTime].emplace_back(CRT::Hit((eDep.GetID())%64, eDep.GetEnergyDeposited()*0.001f*fGeVToADC));
        MF_LOG_DEBUG("TrueTimes") << "Assigned true hit at time " << tAvg << " to bin " << tAvg/fIntegrationTime << ".\n";
      
    }
    }


  //For each CRT module
  for(const auto & crtModule: crtHits)
  {
     int crtChannel=-1;
     int module=crtModule.first;
    
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
    const auto timeToHits=crtModule.second;

    MF_LOG_DEBUG("channels") << "Processing channel " << module << "\n";
    
    std::stringstream ss;
    for(const auto& window: timeToHits)
    {
      ss << "At " << window.first << " ticks, " << window.second.size() << " hits\n";
    } 

    MF_LOG_DEBUG("timeToHits") << "Constructed readout windows for module " << crtChannel << ":\n"
                            << ss.str() << "\n";

   
    MF_LOG_DEBUG("timeToHits") << "About to loop over " << timeToHits.size() << " time windows that are " << fIntegrationTime << "ns long.\n";
    for(auto window = timeToHits.begin(); window != timeToHits.end(); ) //++window)
    {
      const auto& hitsInWindow = window->second;
      const auto aboveThresh = std::find_if(hitsInWindow.begin(), hitsInWindow.end(), 
                                            [this](const auto& hitPair) { return hitPair.ADC() > fDACThreshold; });

      if(aboveThresh != hitsInWindow.end()) //If this is true, then I have found a channel above threshold and readout is triggered.  
      {
        //MF_LOG_DEBUG("timeToHits") << "Channel " << aboveThresh.Channel() << " has deposit " << aboveThresh.ADC() << " ADC counts that " << "is above threshold.  Triggering readout at time " << window->first << ".\n";
 
        std::vector<CRT::Hit> hits;
        const time timestamp = window->first; //Set timestamp before window is changed.
        const auto end = timeToHits.upper_bound(timestamp+fReadoutWindowSize);
        std::set<uint32_t> channelBusy; //A std::set contains only unique elements.  This set stores channels that have already been read out in 
                                        //this readout window and so are "busy" and cannot contribute any more hits.  
        for(;window != end; ++window)
        {
          for(const auto& hitPair: window->second)
          {
            const auto channel = hitPair.Channel();  
                                                                  
            if(channelBusy.insert(channel).second) //If this channel hasn't already contributed to this readout window
            {
              hits.push_back(hitPair);
 
 
              //simToTrigger->addSingle(hitPair.second, makeTrigPtr(trigCol->size()-1)); 
            }
          }
        }


        MF_LOG_DEBUG("CreateTrigger") << "Creating CRT::Trigger...\n";
	//std::cout<< "Timestamp:"<<timestamp*fIntegrationTime<<std::endl;
        trigCol->emplace_back(crtChannel, timestamp*fIntegrationTime, std::move(hits)); 

        const auto oldWindow = window;
        if(window != timeToHits.end()) window = timeToHits.upper_bound(window->first+fDeadtime); 
        if(window != timeToHits.end()) MF_LOG_DEBUG("DeadTime") << "Advanced readout window by " << window->first - oldWindow->first
                                                             << " to simulate dead time.\n";
      } //If there was a channel above threshold
      else ++window; //If discriminators did not fire, continue to next time window.
    } //For each time window
  } //For each CRT module

  //Put Triggers and Assns into the event
  MF_LOG_DEBUG("CreateTrigger") << "Putting " << trigCol->size() << " CRT::Triggers into the event at the end of analyze().\n";
  e.put(std::move(trigCol));
  //e.put(std::move(simToTrigger));
}



DEFINE_ART_MODULE(CRT::CRTSimRefac)
