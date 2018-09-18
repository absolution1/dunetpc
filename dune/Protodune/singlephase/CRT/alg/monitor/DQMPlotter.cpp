//File: DQMPlotter.cpp
//Brief: An "algorithm" class for creating online monitoring plots in a way that works 
//       either with or without LArSoft.  Should start out by plotting diagnostics from 
//       perl plotting code, then branch out as we identify more pathologies.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRT_ONLINEPLOTTER_CPP
#define CRT_ONLINEPLOTTER_CPP

//crt-core includes
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//ROOT includes
#include "TDirectory.h"
#include "TH1D.h"
#include "TGraph.h"

//c++ includes
#include <map>
#include <string>

//TODO: Remove me
#include <iostream>

//Hash function I want to use below
namespace
{
  struct hash
  {
    size_t operator ()(const CRT::Trigger& trigger) const
    {
      return trigger.Channel()/4ul; //Return Trigger's frame.  4 modules per frame.  This has nothing to do with frame numbers on cables!
    }
  };
}

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
      
      DQMPlotter(TFS& tfs, const double timeTickSize): fFileService(tfs), fFirstTime(std::numeric_limits<uint64_t>::max()), 
                                                       fTickSize(timeTickSize)
      {
        fTimestamps = fFileService->template make<TH1D>("timestamps", "Timestamps for All Modules;Timestamp;Triggers", 300, 6.592522e18, 6.592525e18); //Trigger automatic binning?
        fTriggerDeltaT = fFileService->template make<TH1D>("triggerDeltaT", "Time Difference Between Two Triggers on Any Channel;"
                                                                            "#Delta T[ns];Triggrs", 1000, -50, 100);

        fTriggerTime = fFileService->template make<TH1D>("triggerTime", "Time of each Trigger in Run;Time [s];Triggers", 100, 0, 3);

        fDeltaTSeconds =  fFileService->template make<TH1D>("deltaTSeconds", "Difference in Elapsed Time Between Adjacent Triggers;#Delta T [s];"
                                                                             "Triggers", 400, -0.3, 0.3);

        fEventLength = fFileService->template make<TH1D>("EventLength", "Time Ticks Between Earliest and Latest Triggers;#DeltaT [ticks];Events", 
                                                         100, 0, 100);
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
        std::cout << "New event with " << triggers.size() << " Triggers.\n";
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
        }

        //Map from Trigger to times in frame along each axis.  First element is first 2 Triggers (alternates x and y), and second element is 
        //second 2 triggers (alternates y and x).  What the "first" element means is different between frames, but it is consistent within a 
        //frame.  
        ::hash toFrame;
        std::map<size_t, std::pair<std::vector<uint64_t>, std::vector<uint64_t>>> frameToOverlapTimes; 
        uint64_t earliest = std::numeric_limits<uint64_t>::max(), latest = 1e16;
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

          //Plot times of Triggers that overlap with this one, then add this Trigger's time to its frame
          const auto frameNumber = toFrame(trigger);
          const bool isFirst = (trigger.Channel()%4)/2ul; //Integer division determines whether this Trigger is in "first pair" of modules in its frame
          auto& frame = frameToOverlapTimes[frameNumber];
          const auto& overlaps = isFirst?frame.second:frame.first;
          auto& ownAxis = isFirst?frame.first:frame.second;
          std::cout << "Comparing module " << trigger.Channel() << " from frame " << frameNumber << " to " 
                    << overlaps.size() << " overlaps at time " << trigger.Timestamp() << ".  Its own axis has " << ownAxis.size() 
                    << " Triggers already.\n";
          for(const auto& time: overlaps) 
          {
            module.fDeltaTOverlap->Fill(trigger.Timestamp()-time); 
          }
          ownAxis.push_back(trigger.Timestamp());

          if(trigger.Timestamp() > 1e16 && trigger.Timestamp() < earliest) earliest = trigger.Timestamp();
          if(trigger.Timestamp() > latest) latest = trigger.Timestamp();
          
          if(trigger.Timestamp() > 1e16) //TODO: What do the very small timestamps represent?  They have UNIX timestamp of 0
          {
            //Plot lower 32 timestamp bits versus time since beginning of event.
            //To get time since beginning of event, I have to extract each part of the timestamp and convert the lower half to seconds.  
            const double elapsed = (trigger.Timestamp() - fFirstTime)*fTickSize; //In seconds.  
            fTriggerTime->Fill(elapsed);
            fDeltaTSeconds->Fill(elapsed-module.fPrevSeconds); 
            module.fPrevSeconds = elapsed; 

            if(module.fPrevTimestamp < std::numeric_limits<unsigned long long>::max()) //Don't fill fPrevTimestamp if there WAS no 
                                                                                       //fPrevTimestamp
            {
              module.fTriggerDeltaT->Fill(trigger.Timestamp()-module.fPrevTimestamp);
            }
            module.fPrevTimestamp = trigger.Timestamp(); //This also sets the first fPrevTimestamp value
          }

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
            module.fAllADC->Fill(hit.ADC());
          } //For each Hit
          prev = iter;
        } //For each Trigger

        if(triggers.size() > 1) fEventLength->Fill(latest-earliest);
      }

    private:
      TFS fFileService; //Handle to tool for writing out files

      //Plots that will be filled during the lifetime of this object
      //Group plots for a module into one object
      template <class DIRECTORY>
      struct ModulePlots
      {
        ModulePlots(DIRECTORY&& dir): fDir(dir), fPrevTimestamp(std::numeric_limits<decltype(fPrevTimestamp)>::max()), fPrevSeconds(0)
        {
          fHits = fDir.template make<TH1D>("Hits", "Hits for Each Channel;Channel;Hits", 65, 0, 64);
          fTriggerDeltaT = fDir.template make<TH1D>("TriggerDeltaT", "Time Between Triggers;#Delta Timestamp [ticks];Pairs of Triggers", 
                                                    1e3, 0., 1e6);
          fHitsPerTrigger = fDir.template make<TH1D>("HitsPerTrigger", "Hits in each Trigger;Hits;Triggers", 65, 0, 64);

          fAllADC = fDir.template make<TH1D>("AllADC", "ADCs for All Channels in this Module;ADC;Hits", 1024, 0, 4096);

          fDeltaTOverlap = fDir.template make<TH1D>("DeltaTOverlap", "Time Ticks Between Overlapping Modules' Triggers in an Event;"
                                                                     "#DeltaT [ticks];Overlaps", 200, 0, 200);
        }

        DIRECTORY fDir; //Directory to which this module's plots will be written
        TH1D* fHits; //Plot of number of hits in each channel for this module
        TH1D* fTriggerDeltaT; //Plot time between Triggers on this module
        TH1D* fHitsPerTrigger; //Number of Hits per Trigger
        TH1D* fAllADC; //ADC values for all channels in this module
        TH1D* fDeltaTOverlap; //Time difference between pairs of overlapping modules
        std::map<size_t, TH1D*> fChannelToADC; //Plot of ADC counts for each channel
      
        unsigned long long fPrevTimestamp; //Timestamp of the previous Trigger
        double fPrevSeconds; //Previous elapsed time in seconds
                             //timing of sync pulses before resetting.  
      };

      //TODO: It seems like I only need fFileService for this decltype declaration.  Otherwise, I could move from a class template to 
      //      just a function template for the constructor.
      std::vector<ModulePlots<decltype(fFileService->mkdir("test"))>> fModules; //A group of plots for each module
      TH1D* fTimestamps; //Plot of timestampts for all Triggers
      TH1D* fTriggerDeltaT; //Time difference between any Trigger and its' neighbor within an event.  
      TH1D* fTriggerTime; //Elapsed time before each Trigger was detected
      TH1D* fDeltaTSeconds; //Difference between elpased time between Trigger and previous Trigger
      TH1D* fEventLength; //Time between earliest and latest timestamp in each AnalyzeEvent() call in clock ticks

      //Keep track of elapsed time since earliest Trigger
      uint64_t fFirstTime; //First timestamp of any Trigger in first event.
      const double fTickSize; //Size of a time tick in seconds
  };
}

#endif //CRT_ONLINEPLOTTER_CPP
