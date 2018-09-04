////////////////////////////////////////////////////////////////////////
// Class:       CRTSim
// Plugin Type: producer (art v2_10_03)
// File:        CRTSim_module.cc
//
// Generated at Wed Jun 27 04:09:39 2018 by Andrew Olivier using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

//     An ART module to simulate how the ProtoDUNE-SP Cosmic Ray Tagger (CRT) system 
//responds to energy deposits.  Starts with sim::AuxDetSimChannel, loops over SimChannels 
//produced for CRT volumes, groups energy deposits across SimChannels into "readout packets", 
//and produces one CRT::Hit for each (group?  Need more information about board) sim::AuxDetIDE
//within a "readout packet".  
//
//     According to Matt Strait and Vishvas, the CRT boards already do some hit reconstruction for us.  
//This module should apply any energy deposit-level physics (like Birks' Law and Poisson fluctuations 
//in number of photons?), apply detector response (from an as-yet-unavailable MySQL database), then 
//do the same reconstruction steps as the CRT boards.  Each CRT board is responsible for an entire 
//module's channels.  A board only reads out when it finds a "hit" with an energy deposit above some 
//threshold.  Then, it reads out the "hit"s on ALL of its' channels in one "packet" that becomes a 
//single artdaq::Fragment.  

//TODO: art::Assns<sim::AuxDetSimChannel, std::vector<CRT::Hit>> for backtracking?
//TODO: Associate CRT::Hits from a readout window?  Seems possible to do this in "real data".   
//TODO: Produce and associate to raw::ExternalTriggers?  

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
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"

//local includes
//#include "CRTTrigger.h"
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//c++ includes
#include <memory>
#include <algorithm>

namespace CRT {
  class CRTSim;
}

class CRT::CRTSim : public art::EDProducer {
public:
  explicit CRTSim(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTSim(CRTSim const &) = delete;
  CRTSim(CRTSim &&) = delete;
  CRTSim & operator = (CRTSim const &) = delete;
  CRTSim & operator = (CRTSim &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  //void beginJob() override;
  /*void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;*/

private:

  //The formats for time and ADC value I will use throughout this module.  If I want to change it later, I only have to change it in one place.  
  typedef long long int time;
  typedef unsigned short adc_t;

  // Member data
  art::InputTag fSimLabel; //The label of the module that produced sim::AuxDetSimChannels

  //Parameters I could probably get from calibration + electronics documentation to make this simulation more realistic.  
  //Decided that I am more interested in speed for now.
  /*double fScintillationYield; //The converion from GeV to energy to photons
  double fQuantumEff; //The conversion from photons to photo-electrons leaving the photocathode
  double fDummyGain; //Dummy conversion of photo-electrons leaving the photocathode to voltage entering the MAROC2 chip.   
  double fResistanceBeforeMAROC; //Resistance converting current after gain to voltage entering MAROC2 */

  //Simple parameterization in lieu of complicated information above.  I can probably still get this from the calibration output.
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


CRT::CRTSim::CRTSim(fhicl::ParameterSet const & p): fSimLabel(p.get<art::InputTag>("SimLabel")), 
                                                              /*fScintillationYield(p.get<double>("ScintillationYield")), 
                                                              fQuantumEff(p.get<double>("QuantumEff")), 
                                                              fDummyGain(p.get<double>("DummyGain")),*/
                                                              fGeVToADC(p.get<double>("GeVToADC")),
                                                              fIntegrationTime(p.get<time>("IntegrationTime")), 
                                                              fReadoutWindowSize(p.get<size_t>("ReadoutWindowSize")), 
                                                              fDeadtime(p.get<size_t>("Deadtime")),
                                                              fDACThreshold(p.get<adc_t>("DACThreshold"))
{
  //Tell ART that I convert std::vector<AuxDetSimChannel> to CRT::Hits associated with raw::ExternalTriggers
  produces<std::vector<CRT::Trigger>>();
  produces<art::Assns<sim::AuxDetSimChannel, CRT::Trigger>>(); 
  consumes<std::vector<sim::AuxDetSimChannel>>(fSimLabel);
}

//Turn sim::AuxDetSimChannels into CRT::Hits. 
void CRT::CRTSim::produce(art::Event & e)
{
  //Get collection of AuxDetSimChannel as input
  auto channels = e.getValidHandle<std::vector<sim::AuxDetSimChannel>>(fSimLabel);

  //Collections of CRT::Trigger and AuxDetSimChannel-CRT::Trigger associations that I will std::move() into e at the end of produce
  auto trigCol = std::make_unique<std::vector<CRT::Trigger>>();
  auto simToTrigger = std::make_unique<art::Assns<sim::AuxDetSimChannel, CRT::Trigger>>();

  //Utilities to go along with making Assns
  art::PtrMaker<sim::AuxDetSimChannel> makeSimPtr(e, channels.id());
  art::PtrMaker<CRT::Trigger> makeTrigPtr(e, *this);

  //Get access to geometry for each event (TODO: -> subrun?) in case CRTs move later
  art::ServiceHandle<geo::Geometry> geom;

  //Group AuxDetSimChannels by module they came from for easier access later.  Storing AuxDetSimChannels as art::Ptrs so I can make Assns later.
  //TODO: I'm spending a lot of effort and confusing code on keeping track of sim::AuxDetSimChannels for Assns at the end of this module.  
  //      Are raw -> simulation Assns worth this much trouble?  Can I restructure my logic somehow to make this association more natural?  
  std::map<uint32_t, art::PtrVector<sim::AuxDetSimChannel>> moduleToChannels;
  for(size_t channelPos = 0; channelPos < channels->size(); ++channelPos)
  {
    const auto id = (*channels)[channelPos].AuxDetID();
    const auto& det = geom->AuxDet(id);
    if(det.Name().find("CRT") != std::string::npos) //If this is a CRT AuxDet
    {
      moduleToChannels[id].push_back(makeSimPtr(channelPos)); 
    } //End if this is a CRT AuxDet
  } //End for each simulated sensitive volume -> CRT strip

  LOG_DEBUG("moduleToChannels") << "Got " << moduleToChannels.size() << " modules to process in detector simulation.\n";

  //For each CRT module
  for(auto& pair: moduleToChannels)
  {
    auto& module = pair.second;

    //Probably a good idea to write a function/algorithm class for each of these
    //TODO: Any leftover physics like Birks' Law?  I can handle Birks' Law more accurately if I 
    //      can get access to individual Geant steps.  
    //TODO: Read detector response from MariaDB database on DAQ machine?
    //TODO: Simulate detector response with quantum efficiency and detection efficiency?
    //TODO: simulate time -> timestamp.  Check out how 35t example in dunetpc/dune/DetSim does this.

    //Integrate "energy deposited" over time, then form a time-ordered sparse-vector (=std::map) of CRT::Hits
    //associated with the art::Ptr that produced each hit. Associating each hit with a time value gives me the 
    //option to skip over dead time later.  Each time bin contains up to one hit per channel.
    LOG_DEBUG("moduleToChannels") << "Processing " << module.size() << " channels in module " << pair.first << "\n";
    std::map<time, std::vector<std::pair<CRT::Hit, art::Ptr<sim::AuxDetSimChannel>>>> timeToHits; //Mapping from integration time bin to list of hit-Ptr pairs
    for(const auto& channel: module)
    {
      LOG_DEBUG("channels") << "Processing channel " << channel->AuxDetSensitiveID() << "\n";
      const auto& ides = channel->AuxDetIDEs();
      for(const auto& eDep: ides)
      {
        const size_t tAvg = (eDep.exitT+eDep.entryT)/2.;
        timeToHits[tAvg/fIntegrationTime].emplace_back(CRT::Hit(channel->AuxDetSensitiveID(), eDep.energyDeposited*fGeVToADC), channel);
        LOG_DEBUG("TrueTimes") << "Assigned true hit at time " << tAvg << " to bin " << tAvg/fIntegrationTime << ".\n";
      }
    }

    std::stringstream ss;
    for(const auto& window: timeToHits)
    {
      ss << "At " << window.first << " ticks, " << window.second.size() << " hits\n";
    } 

    LOG_DEBUG("timeToHits") << "Constructed readout windows for module " << pair.first  << ":\n"
                            << ss.str() << "\n";

    //Group AuxDetIDEs into CRT board "readout packets".  First, look for an energy 
    //deposit that triggers board readout.  Then, find all energy deposits within 
    //fReadoutWindowSize of each triggering hit and create a "readout packet".  
    //Once an energy deposit has been "read out", it cannot trigger any more 
    //"readout packet"s. 

    //A real CRT module triggers as soon as it sees a channel above threshold. 
    //Time-order the IDEs on each strip, then set triggerTime to the earliest
    //IDE from any channel that is above threshold.  Channel number breaks ties.
   
    //TODO: Eventually read out CRT modules only until one module triggers to simulate the "traffic jam" situation I think could happen?
    LOG_DEBUG("timeToHits") << "About to loop over " << timeToHits.size() << " time windows that are " << fIntegrationTime << "ns long.\n";
    for(auto window = timeToHits.begin(); window != timeToHits.end(); ) //++window)
    {
      const auto& hitsInWindow = window->second;
      const auto aboveThresh = std::find_if(hitsInWindow.begin(), hitsInWindow.end(), 
                                            [this](const auto& hitPair) { return hitPair.first.ADC() > fDACThreshold; });

      if(aboveThresh != hitsInWindow.end()) //If this is true, then I have found a channel above threshold and readout is triggered.  
      {
        LOG_DEBUG("timeToHits") << "Channel " << aboveThresh->first.Channel() << " has deposit " << aboveThresh->first.ADC() << " ADC counts that "
                                << "is above threshold.  Triggering readout at time " << window->first << ".\n";
        //TODO: Integrate all channels over the readout window?  What happens if a channel had energy deposits twice during the readout window?
        //Write all channels with activity in the readout window to a CRT::Trigger.  Ignore repeated hits in a channel for now.   
        std::vector<CRT::Hit> hits;
        const time timestamp = window->first; //Set timestamp before window is changed.
        const auto end = timeToHits.upper_bound(timestamp+fReadoutWindowSize);
        std::set<uint32_t> channelBusy; //A std::set contains only unique elements.  This set stores channels that have already been read out in 
                                        //this readout window and so are "busy" and cannot contribute any more hits.  
        for(;window != end; ++window)
        {
          //TODO: Read out channels, but make sure the same channel is not read out twice because ADC is busy?  Requires nested map?
          for(const auto& hitPair: window->second)
          {
            const auto channel = hitPair.first.Channel(); //TODO: Get channel number back without needing art::Ptr here.  
                                                                    //      Maybe store crt::Hits.
            if(channelBusy.insert(channel).second) //If this channel hasn't already contributed to this readout window
            {
              hits.push_back(hitPair.first);
 
              //Create Assns to AuxDetSimChannels now to avoid having to store art::Ptrs again.  If art::Ptr 
              //actually had to hold a raw pointer to a concrete instance of a class or else be nullptr, then 
              //the following would be gibberish.  However, the comments about art::PtrMaker and what little I 
              //still remember of the "old" method of creating art::Assns in LArSoft lead me to believe that 
              //I am actually working with art::Ptrs that don't exist yet when creating an Assn.  After all, 
              //the art::Event doesn't even own my vector<CRT::Trigger> yet.  So, I can refer to a hypothetical 
              //CRT::Trigger when creating Assns and thus avoid the need to loop over used AuxDetSimChannels 
              //again.  This saves me a horrible STL template instantiation.  
              simToTrigger->addSingle(hitPair.second, makeTrigPtr(trigCol->size()-1)); //TODO: I could add a "data" object that 
                                                                                                   //      gives the indices of the IDEs used 
                                                                                                   //      to create this Trigger
            }
          }
        }

        //Create a CRT::Trigger to represent this module's readout window.  Associate the AuxDetSimChannels used to make this CRT::Trigger.
        //This is the hypothetical CRT::Trigger I was talking about when creating Assns.
        LOG_DEBUG("CreateTrigger") << "Creating CRT::Trigger...\n";
        trigCol->emplace_back(pair.first, timestamp*fIntegrationTime, std::move(hits)); 
                                                                         //TODO: Convert trigger time into timestamp.  From line 211 of 
                                                                         //      DetSim/Modules/SimCounter35t_module.cc, it looks like I 
                                                                         //      need to know more about how CRT timestamps are constructed 
                                                                         //      in data before proceeding.  

        //Advance window past dead time so that no energy deposits in dead time are read out.
        const auto oldWindow = window;
        if(window != timeToHits.end()) window = timeToHits.upper_bound(window->first+fDeadtime); 
        if(window != timeToHits.end()) LOG_DEBUG("DeadTime") << "Advanced readout window by " << window->first - oldWindow->first
                                                             << " to simulate dead time.\n";
      } //If there was a channel above threshold
      else ++window; //If discriminators did not fire, continue to next time window.
    } //For each time window
  } //For each CRT module

  //Put Triggers and Assns into the event
  LOG_DEBUG("CreateTrigger") << "Putting " << trigCol->size() << " CRT::Triggers into the event at the end of analyze().\n";
  e.put(std::move(trigCol));
  e.put(std::move(simToTrigger));
}

//Tell ART what data products this modules works with and retrieve any resources that don't change throughout the job here.
/*void CRT::CRTSim::beginJob()
{
}*/

/*void CRT::CRTSim::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void CRT::CRTSim::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void CRT::CRTSim::endJob()
{
  // Implementation of optional member function here.
}

void CRT::CRTSim::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void CRT::CRTSim::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void CRT::CRTSim::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void CRT::CRTSim::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void CRT::CRTSim::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void CRT::CRTSim::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}*/

DEFINE_ART_MODULE(CRT::CRTSim)
