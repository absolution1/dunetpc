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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//dune-raw-data includes
#include "dune-raw-data/Overlays/CRTFragment.hh"

//dunetpc includes
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

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
      //void beginJob() override;
      //void beginRun(art::Run & r) override;
      //void beginSubRun(art::SubRun & sr) override;
    
    private:
    
      // Declare member data here.
      art::InputTag fFragLabel; //Label of the module that produced artdaq::Fragments 
                                //from CRT data that I will turn into CRT::Triggers.
  
      //TODO: Sync diagnostic plots
      /*struct PerModule
      {
        PerModule(util::Directory parent, const size_t module): fModuleDir(parent.mkdir("Module"+std::to_string(module)))
        {
          fLowerTimeVersusTime = fModuleDir.makeAndRegister<TGraph>("LowerTimeVersusTime", "Raw 32 Bit Timestamp versus Elapsed Time in Seconds;"
                                                                                           "Time [s];Raw Timestamp [ticks]");
        }

        util::Directory fModuleDir; //Directory for plots from this module
        TGraph* fLowerTimeVersusTime; //Graph of lower 32 bits of raw timestamp versus processed timestamp in seconds.  
      }; 

      std::vector<PerModule> fSyncPlots; //Mapping from module number to sync diagnostic plots in a directory 
      uint64_t fEarliestTime; //Earliest time in clock ticks*/
  };
  
  
  CRTRawDecoder::CRTRawDecoder(fhicl::ParameterSet const & p): fFragLabel(p.get<std::string>("RawDataLabel"), p.get<std::string>("RawDataInstance"))
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<CRT::Trigger>>();
  }
  
  //Read artdaq::Fragments produced by fFragLabel, and use CRT::Fragment to convert them to CRT::Triggers.  
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
      const auto& fragHandle = e.getValidHandle<std::vector<artdaq::Fragment>>(fFragLabel);
      
      //Convert each fragment into a CRT::Trigger.
      for(const auto& artFrag: *fragHandle)
      {
        CRT::Fragment frag(artFrag);

        //TODO: Remove me
        LOG_DEBUG("CRTGoodEvent") << "Is this Fragment good?  " << ((frag.good_event())?"true":"false") << "\n";
        frag.print_header();
        frag.hexdump();

        std::vector<CRT::Hit> hits;
  
        //Make a CRT::Hit from each non-zero ADC value in this Fragment
        for(size_t hitNum = 0; hitNum < frag.num_hits(); ++hitNum)
        {
          const auto hit = *(frag.hit(hitNum));
          LOG_DEBUG("CRTRaw") << "Channel: " << (int)(hit.channel) << "\n"
                              << "ADC: " << hit.adc << "\n";
  
          hits.emplace_back(hit.channel, hit.adc);
          //LOG_DEBUG("CRT Hits") CRT::operator << hits.back() << "\n"; //TODO: Some function template from the message service interferes with my  
                                                                        //      function template from namespace CRT.  using namespace CRT seems like 
                                                                        //      it should solve this, but it doesn't seem to.
        }
  
        LOG_DEBUG("CRTFragments") << "Module: " << frag.module_num() << "\n"
                                  << "Number of hits: " << frag.num_hits() << "\n"
                                  << "Fifty MHz time: " << frag.fifty_mhz_time() << "\n";
  
        triggers->emplace_back(frag.module_num(), frag.fifty_mhz_time(), std::move(hits)); //TODO: Get AuxDet index from channel map
        
        /* //TODO: Get raw lower timestamp for this plot
        //Make diagnostic plots for sync pulses
        const auto& plots = fSyncPlots[trigger.Channel()];
        const double deltaT = (triggers.back().Timestamp() - fEarliestTime)*1.6e-8;
        if(deltaT > 0 && deltaT < 1e6) //Ignore time differences less than 1s and greater than 1 day
                                       //TODO: Understand why these cases come up
        {
          plots.fLowerTimeVersusTime->SetPoint(plots.fLowerTimeVersusTime->GetN(), deltaT, pair.second*16);
        }
        else
        {
          //TODO: Print something out here
        }*/
      } 
    }
    catch(const cet::exception& exc) //If there are no artdaq::Fragments in this Event, just add an empty container of CRT::Triggers.
    {
      mf::LogWarning("MissingData") << "No artdaq::Fragments produced by " << fFragLabel << " in this event, so not doing anything.\n";
    }
  
    //Put a vector of CRT::Triggers into this Event for other modules to read.
    e.put(std::move(triggers));
  }
  
  /*void CRT::CRTRawDecoder::beginJob()
  {
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    for(size_t module = 0; module < 32; ++module) fSyncPlots.emplace_back(tfs, module); //32 modules in the ProtoDUNE-SP CRT
  }
  
  void CRT::CRTRawDecoder::beginRun(art::Run & r)
  {
    // Implementation of optional member function here.
  }
  
  void CRT::CRTRawDecoder::beginSubRun(art::SubRun & sr)
  {
    // Implementation of optional member function here.
  }*/
}

DEFINE_ART_MODULE(CRT::CRTRawDecoder)
