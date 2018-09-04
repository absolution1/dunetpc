//File: DQMPlotter.cpp
//Brief: An "algorithm" class for creating online monitoring plots in a way that works 
//       either with or without LArSoft.  Should start out by plotting diagnostics from 
//       perl plotting code, then branch out as we identify more pathologies.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRT_ONLINEPLOTTER_CPP
#define CRT_ONLINEPLOTTER_CPP

//crt-core includes
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"
#include "../srcs/dunetpc/dune/Protodune/singlephase/CRT/alg/util/Timer.h"

//ROOT includes
#include "TDirectory.h"
#include "TH1D.h"
#include "TGraph.h"

//c++ includes
#include <map>
#include <string>

namespace CRT
{
  template <class TFS> 
  //TFS is an object that defines a function like:
  //
  //template <class HIST, class ARGS...> HIST* make(ARGS... args)
  //
  //Looks like ART's TFileService.  
  class DQMPlotter
  {
    public:
      
      DQMPlotter(TFS& tfs): fFileService(tfs), fFirstTime(std::numeric_limits<uint64_t>::max()), fElapsed(nullptr)
      {
        fTimestamps = fFileService->template make<TH1D>("timestamps", "Timestamps for All Modules;Timestamp;Triggers", 1000, 6.58923730e18, 6.58923750e18); //Trigger automatic binning?
        fTriggerDeltaT = fFileService->template make<TH1D>("triggerDeltaT", "Time Difference Between Two Triggers on Any Channel;"
                                                                            "#Delta T;Triggrs", 75, -50, 200);
        for(size_t module = 0; module < 32; ++module) //TODO: 32 is the number of modules in the CRT.  Replace it with  
                                                      //      some more descriptive constant.
        {
          const auto name = "module"+std::to_string(module);
          fModules.emplace_back(fFileService->mkdir(name));
        }
      }

      virtual ~DQMPlotter() = default; //All pointers kept by this class are owned 
                                          //by another object.  

      
      void AnalyzeEvent(const std::vector<CRT::Trigger>& triggers) //Make plots from CRT::Triggers
      {
        const auto end = triggers.cend();
        auto prev = end;

        //If first timestamp hasn't been set to a suitable value yet, try to set it now to the 
        //earliest timestamp in this event.
        if(fFirstTime == std::numeric_limits<uint64_t>::max())
        {
          for(const auto& trigger: triggers) 
          {
            if(trigger.Timestamp() > 1e18 && trigger.Timestamp() < fFirstTime) 
            {
              fFirstTime =  trigger.Timestamp();
            }
          }
          fElapsed.reset(new CRT::Timer(fFirstTime));
        }

        for(auto iter = triggers.cbegin(); iter != end; ++iter)
        { 
          const auto& trigger = *iter;
          fTimestamps->Fill(trigger.Timestamp());
          if(prev != end)
          {
            fTriggerDeltaT->Fill(trigger.Timestamp()-prev->Timestamp());
          }

          const auto id = trigger.Channel();
          auto& module = fModules[id]; 
          
          if(trigger.Timestamp() > 1e16) //TODO: What do the very small timestamps represent?  They have UNIX timestamp of 0
          {
            //Plot lower 32 timestamp bits versus time since beginning of event.
            //To get time since beginning of event, I have to extract each part of the timestamp and convert the lower half to seconds.  
            const double elapsed = fElapsed->elapsed(trigger.Timestamp());
            module.fLowerTimeVersusTime->SetPoint(module.fLowerTimeVersusTime->GetN(), elapsed, (uint32_t)trigger.Timestamp()); 
          }

          if(module.fPrevTimestamp < std::numeric_limits<unsigned long long>::max()) //Don't fill fPrevTimestamp if there WAS no 
                                                                                     //fPrevTimestamp
          {
            module.fTriggerDeltaT->Fill(trigger.Timestamp()-module.fPrevTimestamp);
          }
          module.fPrevTimestamp = trigger.Timestamp(); //This also sets the first fPrevTimestamp value

          const auto& hits = trigger.Hits();
          module.fHitsPerTrigger->Fill(hits.size());
          for(const auto& hit: hits)
          {
            const auto channel = hit.Channel();
            module.fHits->Fill(channel);

            const auto adcFound = module.fChannelToADC.find(channel);
            if(adcFound == module.fChannelToADC.end())
            {
              const auto name = "channel"+std::to_string(channel);
              module.fChannelToADC[channel] = module.fDir.template make<TH1D>(name.c_str(), ("Hits on Channel "
                                                                              +std::to_string(channel)).c_str(), 
                                                                              500, 0, 2000);
            } 

            module.fChannelToADC[channel]->Fill(hit.ADC());
          } //For each Hit
          prev = iter;
        } //For each Trigger
      }

    private:
      TFS fFileService; //Handle to tool for writing out files

      //Plots that will be filled during the lifetime of this object
      //Group plots for a module into one object
      template <class DIRECTORY>
      struct ModulePlots
      {
        ModulePlots(DIRECTORY&& dir): fDir(dir), fPrevTimestamp(std::numeric_limits<decltype(fPrevTimestamp)>::max()) 
        {
          fHits = fDir.template make<TH1D>("Hits", "Hits for Each Channel;Channel;Hits", 65, 0, 64);
          fTriggerDeltaT = fDir.template make<TH1D>("TriggerDeltaT", "Time Between Triggers;#Delta Timestamp;Pairs of Triggers", 
                                                    100, 0., 100.);
          fHitsPerTrigger = fDir.template make<TH1D>("HitsPerTrigger", "Hits in each Trigger;Hits;Triggers", 65, 0, 64);

          fLowerTimeVersusTime = fDir.template makeAndRegister<TGraph>("LowerTimeVersusTime", "Lower 32 Timestamp Bits versus Total Time in Job;"
                                                                       "Total Time in Job [s];Lower 32 Bits [ns]");
        }

        DIRECTORY fDir; //Directory to which this module's plots will be written
        TH1D* fHits; //Plot of number of hits in each channel for this module
        TH1D* fTriggerDeltaT; //Plot time between Triggers on this module
        TH1D* fHitsPerTrigger; //Number of Hits per Trigger
        std::map<size_t, TH1D*> fChannelToADC; //Plot of ADC counts for each channel
      
        unsigned long long fPrevTimestamp; //Timestamp of the previous Trigger
        
        //Plots made over time
        TGraph* fLowerTimeVersusTime; //Plot lower 32 bits of timestamp as a function of total timestamp to see lower 32 bits increase to 
                                      //timing of sync pulses before resetting.  
      };

      //TODO: It seems like I only need fFileService for this decltype declaration.  Otherwise, I could move from a class template to 
      //      just a function template for the constructor.
      std::vector<ModulePlots<decltype(fFileService->mkdir("test"))>> fModules; //A group of plots for each module
      TH1D* fTimestamps; //Plot of timestampts for all Triggers
      TH1D* fTriggerDeltaT; //Time difference between any Trigger and its' neighbor within an event.  

      //Keep track of elapsed time since earliest Trigger
      uint64_t fFirstTime; //First timestamp of any Trigger in first event.
      std::unique_ptr<CRT::Timer> fElapsed; //TODO: Hold a pointer to fElapsed so I can delay its' initialization to the first event
  };
}

#endif //CRT_ONLINEPLOTTER_CPP
