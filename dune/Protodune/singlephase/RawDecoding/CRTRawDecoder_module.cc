////////////////////////////////////////////////////////////////////////
// Class:       CRTRawDecoder
// Plugin Type: producer (art v2_10_03)
// File:        CRTRawDecoder_module.cc
//
// Generated at Wed Jul  4 08:14:43 2018 by Andrew Olivier using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//dune-artdaq includes
#include "artdaq-core/Data/ContainerFragment.hh"

//dune-raw-data includes
#include "dune-raw-data/Overlays/CRTFragment.hh"

//dunetpc includes
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//ROOT includes
#include "TGraph.h"

//c++ includes
#include <memory>

using namespace CRT;

//A CRTRawDecoder takes artdaq::Fragments made from Cosmic Ray Tagger input 
//and produces a CRT::Trigger for each time a CRT module triggered in this Event.  
namespace CRT
{
  class CRTRawDecoder : public art::EDProducer 
  {
    public:
      explicit CRTRawDecoder(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.
    
      // Plugins should not be copied or assigned.
      CRTRawDecoder(CRTRawDecoder const &) = delete;
      CRTRawDecoder(CRTRawDecoder &&) = delete;
      CRTRawDecoder & operator = (CRTRawDecoder const &) = delete;
      CRTRawDecoder & operator = (CRTRawDecoder &&) = delete;
    
      // Required functions.
      void produce(art::Event & e) override;
    
      // Selected optional functions.  I might get services that I expect to change here.  So far, seems like only service I might need is timing.
      void beginJob() override;
      //void beginRun(art::Run & r) override;
      //void beginSubRun(art::SubRun & sr) override;
    
    private:
    
      // Declare member data here.
      const art::InputTag fFragTag; //Label and product instance name (try eventdump.fcl on your input file) of the module that produced 
                                      //artdaq::Fragments from CRT data that I will turn into CRT::Triggers. Usually seems to be "daq" from artdaq.  

      const bool fLookForContainer; //Should this module look for "container" artdaq::Fragments instead of "single" artdaq::Fragments?  
                                    //If you ran the CRT board reader with the "request_mode" FHICL parameter set to "window", you need 
                                    //to look for "container" Fragments and probably won't get any useful Fragments otherwise.  If you 
                                    //ran the CRT board reader with "request_mode" set to "ignore", then set fLookForContainer to "false" 
                                    //which is the default.  
                                    //You probably want this set to "true" because the CRT should be run in "window" request_mode for 
                                    //production.  You might set it to "false" if you took debugging data in which you wanted to ignore 
                                    //timestamp matching with other detectors.  

      // Compartmentalize internal functionality so that I can reuse it with both regular Fragments and "container" Fragments
      void FragmentToTriggers(const artdaq::Fragment& artFrag, std::unique_ptr<std::vector<CRT::Trigger>>& triggers);

      // For the first Event of every job, I want to set fEarliestTime to the earliest time in that Event
      void SetEarliestTime(const artdaq::Fragment& frag);

      //Sync diagnostic plots
      struct PerModule
      {
        PerModule(art::TFileDirectory& parent, const size_t module): fModuleDir(parent.mkdir("Module"+std::to_string(module)))
        {
          fLowerTimeVersusTime = fModuleDir.makeAndRegister<TGraph>("LowerTimeVersusTime", "Raw 32 Bit Timestamp versus Elapsed Time in Seconds;"
                                                                                           "Time [s];Raw Timestamp [ticks]");
        }

        art::TFileDirectory fModuleDir; //Directory for plots from this module
        TGraph* fLowerTimeVersusTime; //Graph of lower 32 bits of raw timestamp versus processed timestamp in seconds.  
      }; 

      void createSyncPlots();

      std::vector<PerModule> fSyncPlots; //Mapping from module number to sync diagnostic plots in a directory 
      uint64_t fEarliestTime; //Earliest time in clock ticks
  };
  
  
  CRTRawDecoder::CRTRawDecoder(fhicl::ParameterSet const & p): fFragTag(p.get<std::string>("RawDataTag")), 
                                                               fLookForContainer(p.get<bool>("LookForContainer", false)),
                                                               fEarliestTime(std::numeric_limits<decltype(fEarliestTime)>::max())
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<CRT::Trigger>>();
    consumes<std::vector<artdaq::Fragment>>(fFragTag);
 
    //Register callback to make new plots on every file
    art::ServiceHandle<art::TFileService> tfs;
    tfs->registerFileSwitchCallback(this, &CRTRawDecoder::createSyncPlots);
  }

  void CRTRawDecoder::FragmentToTriggers(const artdaq::Fragment& artFrag, std::unique_ptr<std::vector<CRT::Trigger>>& triggers)
  {
    CRT::Fragment frag(artFrag);
                                                                                                                                                   
    MF_LOG_DEBUG("CRTGoodEvent") << "Is this Fragment good?  " << ((frag.good_event())?"true":"false") << "\n";
    /*frag.print_header();
    frag.print_hits();*/
                                                                                                                                                   
    std::vector<CRT::Hit> hits;
                                                                                                                                                   
    //Make a CRT::Hit from each non-zero ADC value in this Fragment
    for(size_t hitNum = 0; hitNum < frag.num_hits(); ++hitNum)
    {
      const auto hit = *(frag.hit(hitNum));
      MF_LOG_DEBUG("CRTRaw") << "Channel: " << (int)(hit.channel) << "\n"
                          << "ADC: " << hit.adc << "\n";
                                                                                                                                                   
      hits.emplace_back(hit.channel, hit.adc);
      //MF_LOG_DEBUG("CRT Hits") CRT::operator << hits.back() << "\n"; //TODO: Some function template from the message service interferes with my  
                                                                    //      function template from namespace CRT.  using namespace CRT seems like 
                                                                    //      it should solve this, but it doesn't seem to.
    }
                                                                                                                                                   
    MF_LOG_DEBUG("CRTFragments") << "Module: " << frag.module_num() << "\n"
                              << "Number of hits: " << frag.num_hits() << "\n"
                              << "Fifty MHz time: " << frag.fifty_mhz_time() << "\n";
                                                                                                                                                   
    triggers->emplace_back(frag.module_num(), frag.fifty_mhz_time(), std::move(hits)); //TODO: Get AuxDet index from channel map
    
    //Make diagnostic plots for sync pulses
    const auto& plots = fSyncPlots[frag.module_num()];
    const double deltaT = (frag.fifty_mhz_time() - fEarliestTime)*1.6e-8; //TODO: Get size of clock ticks from a service
    if(deltaT > 0 && deltaT < 1e6) //Ignore time differences less than 1s and greater than 1 day
                                   //TODO: Understand why these cases come up
    {
      plots.fLowerTimeVersusTime->SetPoint(plots.fLowerTimeVersusTime->GetN(), deltaT, frag.raw_backend_time());
    }
    else
    {
      mf::LogWarning("SyncPlots") << "Got time difference " << deltaT << " that was not included in sync plots.\n"
                                  << "lhs is " << frag.fifty_mhz_time() << ", and rhs is " << fEarliestTime << ".\n"
                                  << "Hardware raw time is " << frag.raw_backend_time() << ".\n";
    }
  }

  void CRTRawDecoder::SetEarliestTime(const artdaq::Fragment& frag)
  {
    CRT::Fragment crt(frag);
    if(crt.fifty_mhz_time() < fEarliestTime) fEarliestTime = crt.fifty_mhz_time();
  }
  
  //Read artdaq::Fragments produced by fFragTag, and use CRT::Fragment to convert them to CRT::Triggers.  
  void CRTRawDecoder::produce(art::Event & e)
  {
    //Create an empty container of CRT::Triggers.  Any Triggers in this container will be put into the 
    //event at the end of produce.  I will try to fill this container, but just not produce any CRT::Triggers 
    //if there are no input artdaq::Fragments.  
    auto triggers = std::make_unique<std::vector<CRT::Trigger>>();

    try
    {
      //Try to get artdaq::Fragments produced from CRT data.  The following line is the reason for 
      //this try-catch block.  I don't expect anything else to throw a cet::Exception.
      const auto& fragHandle = e.getValidHandle<std::vector<artdaq::Fragment>>(fFragTag);

      if(fLookForContainer)
      {
        //If this is the first event, set fEarliestTime
        if(fEarliestTime == std::numeric_limits<decltype(fEarliestTime)>::max())
        {
          for(const auto& frag: *fragHandle) 
          {
            artdaq::ContainerFragment container(frag);
            for(size_t pos = 0; pos < container.block_count(); ++pos) SetEarliestTime(*container[pos]);
          }
        }

        for(const auto& artFrag: *fragHandle)
        {
          artdaq::ContainerFragment container(artFrag);
          for(size_t pos = 0; pos < container.block_count(); ++pos) FragmentToTriggers(*container[pos], triggers); 
        }
      }
      else
      {
        //If this is the first event, set fEarliestTime
        if(fEarliestTime == std::numeric_limits<decltype(fEarliestTime)>::max())
        {
          for(const auto& frag: *fragHandle) SetEarliestTime(frag);
        }

        //Convert each fragment into a CRT::Trigger.
        for(const auto& artFrag: *fragHandle) FragmentToTriggers(artFrag, triggers);
      }
    }
    catch(const cet::exception& exc) //If there are no artdaq::Fragments in this Event, just add an empty container of CRT::Triggers.
    {
      mf::LogWarning("MissingData") << "No artdaq::Fragments produced by " << fFragTag << " in this event, so "
                                    << "not doing anything.\n";
    }

    //Put a vector of CRT::Triggers into this Event for other modules to read.
    e.put(std::move(triggers));
  }
  
  void CRT::CRTRawDecoder::beginJob()
  {
    createSyncPlots();
  }

  void CRT::CRTRawDecoder::createSyncPlots()
  {
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fSyncPlots.clear(); //Clear any previous pointers to histograms.  Some TFile owned them.  
    for(size_t module = 0; module < 32; ++module) fSyncPlots.emplace_back(*tfs, module); //32 modules in the ProtoDUNE-SP CRT
  }
  
  /*void CRT::CRTRawDecoder::beginRun(art::Run & r)
  {
    // Implementation of optional member function here.
  }
  
  void CRT::CRTRawDecoder::beginSubRun(art::SubRun & sr)
  {
    // Implementation of optional member function here.
  }*/
}

DEFINE_ART_MODULE(CRT::CRTRawDecoder)
