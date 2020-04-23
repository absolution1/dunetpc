//File: ValidationPlotter.cpp
//Brief: An "algorithm" class for creating online monitoring plots in a way that works 
//       either with or without LArSoft.  Should start out by plotting diagnostics from 
//       perl plotting code, then branch out as we identify more pathologies.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRT_ONLINEPLOTTER_CPP
#define CRT_ONLINEPLOTTER_CPP

//crt-core includes
#include "dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//crt-alg includes
#include "dune/Protodune/singlephase/CRT/alg/geom/Geometry.h"

//ROOT includes
#include "TDirectory.h"
#include "TH1D.h"

//c++ includes
#include <map>
#include <string>
#include <algorithm> //For std::max_element()

//TODO: Remove me
#include <iostream>

namespace CRT
{
  template <class TFS> 
  //TFS is an object that defines a function like:
  //
  //template <class HIST, class ARGS...> HIST* make(ARGS... args)
  //
  //Looks like ART's TFileService.  
  class ValidationPlotter
  {
    public:
      
      ValidationPlotter(TFS& tfs, std::unique_ptr<CRT::Geometry>&& geom): fFileService(tfs), 
                                                                          fGeom(std::move(geom))
      {
      }

      virtual ~ValidationPlotter() = default; //All pointers kept by this class are owned 
                                          //by another object.  

      
      void AnalyzeEvent(const std::vector<CRT::Trigger>& triggers) //Make plots from CRT::Triggers
      {
        std::cout << "New event\n";
        //First, fill a nested map to ADC values to find overlaps.  
        CRT::geoMap<int16_t> stripToHits;

        for(const auto& trigger: triggers) 
        {
          std::cout << "On module " << trigger.Channel() << ":\n";
          const auto& hits = trigger.Hits();
          for(const auto& hit: hits)
          {
            std::cout << "\tGot a hit on strip " << hit.Channel() << "\n";
            const auto id = fGeom->StripID(trigger.Channel(), hit.Channel());
            
            //Create a place to put plots from this (strip, module) pair in case it will be needed later.  
            auto foundModule = fModuleToDir.find(id);
            if(foundModule == fModuleToDir.end()) fModuleToDir.emplace(id, fFileService->mkdir("Module"+std::to_string(trigger.Channel())));

            auto foundStrip = fStripToChannel.find(id);
            if(foundStrip == fStripToChannel.end()) fStripToChannel.emplace(id, hit.Channel());

            if(hit.ADC() > stripToHits[id]) stripToHits[id] = hit.ADC(); //TODO: Assuming for now that there is no more 
                                                                         //      than 1 hit for each strip in an event.  
                                                                         //      If there were more than one hit, I'd 
                                                                         //      keep the largest here because I am 
                                                                         //      looking for the largest ADC on each 
                                                                         //      plane for matching later.  
          }
        }

        //Next, find the highest ADC hit in each layer and look for an overlap with the other layer in its module.
        //Fill a map with all such overlaps sorted by plane so that I can find 2D pixels made of pairs of highest-energy 
        //hits.  
        CRT::map<std::vector<std::pair<CRT::StripID, int64_t>>, CRT::FrameID, CRT::PlaneID> stripToOverlappingHits;
        for(const auto& framePair: stripToHits)
        {
          //If there are 2 planes in this frame with hits, then 
          //there could be 2D overlaps in this frame.
          if(framePair.second.size() > 1)
          {
            for(const auto& planePair: framePair.second)
            {
              for(const auto& modulePair: planePair.second)
              {
                //There have to be 2 layers with hits in this plane for overlaps to be possible.
                if(modulePair.second.size() > 1)
                {
                  auto layerPair = modulePair.second.begin();
                  auto& stripPairs = layerPair->second;
                  const auto firstMax = std::max_element(stripPairs.begin(), stripPairs.end(), [](const auto& first, const auto& second)
                                                                                               {
                                                                                                 return first.second < second.second; 
                                                                                               });

                  ++layerPair; //Go to the other layer in this plane
                  auto& secondStripPairs = layerPair->second;
                  const auto secondMax = std::max_element(secondStripPairs.begin(), secondStripPairs.end(), [](const auto& first, const auto& second)
                                                                                                {
                                                                                                  return first.second < second.second; 
                                                                                                });

                  if(firstMax->first.Overlaps(secondMax->first)) 
                  {
                    const auto& firstID = firstMax->first;
                    stripToOverlappingHits[firstID].push_back(*firstMax);

                    const auto& secondID = secondMax->first;
                    stripToOverlappingHits[secondID].push_back(*secondMax);
                  } //If firstMax->first overlaps secondMax->first
                } //If more than one layer has hits
              } //For each module
            } //For each plane in this frame
          } //If more than one plane has hits
        } //For each frame with hits

        //Finally, plot the ADCs of strips that are part of 2D pixels
        for(const auto& framePair: stripToOverlappingHits)
        {
          if(framePair.second.size() > 1) //If there is more than one plane with overlapping hit pairs within a frame, 
                                          //that plane has 2D pixels.
          {
            for(const auto& planePair: framePair.second)
            {
              const std::vector<std::pair<CRT::StripID, int64_t>>& strips = planePair.second; 
              for(const auto& stripPair: strips) //TODO: Why doesn't planePair.second work here?! 
              {
                //Look for plots for this StripID and create them if necessary.  This way, 
                //I'll only get plots for modules that were part of overlapping pairs.  
                const auto& id = stripPair.first;
                auto foundStrip = fStrips.find(stripPair.first);
                if(foundStrip == fStrips.end())
                {
                  foundStrip = fStrips.emplace(stripPair.first, StripPlots(fModuleToDir.at(id), fStripToChannel.at(id))).first;
                }

                auto& strip = foundStrip->second; 
                strip.fOverlappingADCs->Fill(stripPair.second);
              } //For each strip with an overlap in this plane
            } //For each plane in this frame
          } //If more than one plane
        } //For each frame with pairs of overlapping strips
      } //AnalyzeEvent()

    private:
      TFS fFileService; //Handle to tool for writing out files
      using DIRECTORY = decltype(fFileService->mkdir("null"));

      //Group plots for a module into one object
      struct StripPlots
      {
        StripPlots(DIRECTORY& dir, const size_t channel)
        {
          fOverlappingADCs = dir.template make<TH1D>(("channel"+std::to_string(channel)+"OverlappingADCs").c_str(), 
                                                     "ADC Values for 2D-Matched Hits;ADC;Hits", 
                                                     1024, 0, 4095);
        }

        TH1D* fOverlappingADCs; //ADC values for this strip when it overlaps with something
      };

      //Plots produced
      std::map<CRT::StripID, StripPlots> fStrips; //Plots for each strip

      //Configuation state
      std::unique_ptr<CRT::Geometry> fGeom; //Handle to geometrical description of CRT

      //Hack to get around lack of backwards mapping in Geometry interface
      std::map<CRT::ModuleID, DIRECTORY> fModuleToDir; //DIRECTORY-channel number pair for each StripID processed so far.  
      std::map<CRT::StripID, size_t> fStripToChannel; //Mapping from strip number to electronics channel number within its board
  };
}

#endif //CRT_ONLINEPLOTTER_CPP
