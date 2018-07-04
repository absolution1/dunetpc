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

//dunetpc includes
#include "dunetpc/dune/Protodune/singlephase/CRTSim/CRTTrigger.h"

//c++ includes
#include <memory>

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

  // Member data
  art::InputTag fSimLabel; //The label of the module that produced sim::AuxDetSimChannels
  size_t fReadoutWindowSize; //The size of a CRT readout window in ns
  const double fDummyTriggerThreshold; //Dummy threshold for triggering readout for any CRT strip.  
                                       //In GeV for now, but needs to become ADC counts one day.  
                                       //Should be replaced by either a lookup in a hardware database or 
                                       //some constant value one day.  
};


CRT::CRTSim::CRTSim(fhicl::ParameterSet const & p): fSimLabel(p.get<art::InputTag>("SimLabel")), 
                                                              fReadoutWindowSize(p.get<double>("ReadoutWindowSize")), 
                                                              fDummyTriggerThreshold(p.get<double>("DummyTriggerThreshold"))
{
  //Tell ART that I convert std::vector<AuxDetSimChannel> to CRT::Hits associated with raw::ExternalTriggers
  produces<std::vector<CRT::Trigger>>();
  produces<art::Assns<sim::AuxDetSimChannel, CRT::Trigger>>(); 
  consumes<std::vector<sim::AuxDetSimChannel>>();
}

//Turn sim::AuxDetSimChannels into CRT::Hits. 
void CRT::CRTSim::produce(art::Event & e)
{
  //Get collection of AuxDetSimChannel as input
  const auto channels = e.getValidHandle<std::vector<sim::AuxDetSimChannel>>(fSimLabel);

  //Collections of raw::ExternalTrigger and CRT::Hits that I will std::move() into e at the end of produce
  std::unique_ptr<std::vector<CRT::Trigger>> trigCol;
  std::unique_ptr<art::Assns<sim::AuxDetSimChannel, CRT::Trigger>> simToTrigger;

  //Utilities to go along with making Assns
  art::PtrMaker<sim::AuxDetSimChannel> makeSimPtr(e, channels.id());
  art::PtrMaker<CRT::Trigger> makeTrigPtr(e, *this);

  //Get access to geometry for each event (TODO: -> subrun?) in case CRTs move later
  art::ServiceHandle<geo::Geometry> geom;

  //Group AuxDetSimChannels by module they came from for easier access later
  std::unordered_map<uint32_t, std::list<sim::AuxDetSimChannel>> moduleToChannel;
  for(const auto& channel: *channels)
  {
    const auto id = channel.AuxDetID();
    const auto& det = geom->AuxDet(id);
    if(det.Name().find("CRT") != std::string::npos) //If this is a CRT AuxDet
    {
      moduleToChannel[id].push_back(moduleToChannel); //At this point, I have made my own copy of each AuxDetSimChannel to modify as I please
    } //If this is a CRT AuxDet
  } //For each simulated sensitive volume -> CRT strip

  //For each CRT module
  for(auto& pair: moduleToChannel)
  {
    auto& module = pair.second;

    //Probably a good idea to write a function/algorithm class for each of these
    //TODO: Any leftover physics like Birks' Law?  
    //TODO: Read detector response from MariaDB database on DAQ machine?
    //TODO: Simulate detector response with quantum efficiency and detection efficiency?
    //TODO: simulate time -> timestamp.  Check out how 35t example in dunetpc/dune/DetSim does this.

    const double threshold = fDummyTriggerThreshold; //Threshold for an energy deposit to trigger readout in this strip.  In GeV for now, but 
                                                     //ultimately becomes a calibrated ADC value that might be different for each strip.

    //Group AuxDetIDEs into CRT board "readout packets".  First, look for an energy 
    //deposit that triggers board readout.  Then, find all energy deposits within 
    //fReadoutWindowSize of each triggering hit and create a "readout packet".  
    //Once an energy deposit has been "read out", it cannot trigger any more 
    //"readout packet"s. 

    //A real CRT module triggers as soon as it sees a channel above threshold. 
    //Time-order the IDEs on each strip, then set triggerTime to the earliest
    //IDE from any channel that is above threshold.  Channel number breaks ties.
    double triggerTime = std::numeric_limits<double>::max();
    bool triggered = true; //TODO: Once I know the readout window size, I shouldn't need this variable. 

    for(auto& strip: module) //For each strip in this module
    {
      auto& ides = strip.AuxDetIDEs();

      //Find first IDE that is above trigger threshold
      auto earliest = ides.end();
      for(auto ide = ides.begin(); ide != ides.end(); ++ide)
      {
        if(ide->energyDeposit > threshold) //If above threshold
        {
          if(earliest == ides.end() || ide->StartT < earliest->StartT) earliest = ide; //If there is no above threshold IDE or this IDE is 
                                                                                       //earlier, make it the earliest IDE.  
        }
      }

      //Record this strip's earliest time if it is earliest energy deposit of all strips so far
      if(earliest->StartT < triggerTime) 
      {
        triggerTime = earliest->StartT;
        triggered = true;
      }
    }

    if(triggered) //TODO: If found a trigger time.  Probably easy to do once I know the length of the readout window
    {
      //Find all energy deposits that happened within a CRT readout window of this hit.  Create a CRT::Hit for each.  Remove read out IDEs 
      //from my local list of IDEs so that they do not appear in other hits.   
      std::vector<CRT::Hit> hits;
      for(auto& strip: module)
      { 
        strip.remove_if([*this, &triggerTime, &hits, &strip](const auto& ide)
                        {
                          if(strip.StartT - triggerTime >= fReadoutWindowSize) return false; //If this IDE is not in the trigger window, move on
                          
                          //Otherwise, remove this IDE from the list of IDEs for future triggers and create a hit for it
                          hits.emplace_back(strip.AuxDetSensitiveID(), ide.energyDeposit); //TODO: turn energyDeposit into ADC hits
                          return true;
                        });
      }
       
      //Create a CRT::Trigger to represent this module's readout window.  Associate the AuxDetSimChannels used to make this CRT::Trigger.
      trigCol->emplace_back(module.AuxDetID(), triggerTime, std::move(hits)); //TODO: Convert trigger time into timestamp.  From line 211 of 
                                                                              //      DetSim/Modules/SimCounter35t_module.cc, it looks like I 
                                                                              //      need to know more about how CRT timestamps are constructed 
                                                                              //      in data before proceeding.  
      simToTrigger->AddSingle(makeSimPtr(module), makeTrigPtr(trigCol->size()-1));
    }
  } //For each CRT module

  //Put Triggers and Assns into the event
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
