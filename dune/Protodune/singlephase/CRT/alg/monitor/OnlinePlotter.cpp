//File: OnlinePlotter.cpp
//Brief: An "algorithm" class for creating online monitoring plots in a way that works 
//       either with or without LArSoft.  Should start out by plotting diagnostics from 
//       perl plotting code, then branch out as we identify more pathologies.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRT_ONLINEPLOTTER_CPP
#define CRT_ONLINEPLOTTER_CPP

//crt-core includes
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//ROOT includes
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"

//c++ includes
#include <unordered_map> 
#include <string>
#include <iostream> //TODO: Remove me

namespace CRT
{
  template <class TFS> 
  //TFS is an object that defines a function like:
  //
  //template <class HIST, class ARGS...> HIST* make(ARGS... args)
  //
  //Looks like ART's TFileService.  
  class OnlinePlotter
  {
    public:
      
      OnlinePlotter(TFS& tfs): fRunNum(0), fFileService(tfs), fRunStartTime(std::numeric_limits<double>::max()), fRunStopTime(0), 
                               fStartTotalTime(std::numeric_limits<double>::max()), fStopTotalTime(0),
                                                              fModuleToUSB ({{0, 13}, 
                                                                             {1, 13},
                                                                             {2, 13},
                                                                             {3, 13},
                                                                             {4, 13},
                                                                             {5, 13},
                                                                             {6, 13},
                                                                             {7, 13},
                                                                             {8, 14},
                                                                             {9, 14},
                                                                             {10, 14},
                                                                             {11, 14},
                                                                             {12, 14}, 
                                                                             {13, 14},
                                                                             {14, 14},
                                                                             {15, 14},
                                                                             {16, 3},
                                                                             {17, 3},
                                                                             {18, 3},
                                                                             {19, 3},
                                                                             {20, 22},
                                                                             {21, 22},
                                                                             {22, 22},
                                                                             {23, 22},
                                                                             {24, 22},
                                                                             {25, 22},
                                                                             {26, 22},
                                                                             {27, 22},
                                                                             {28, 3},
                                                                             {29, 3},
                                                                             {30, 3},
                                                                             {31, 3}})
      {
         //TODO: Get unordered_mapping from module to USB from some parameter passed to constructor
         //TODO: Get tick length from parameter passed to constructor
      }

      virtual ~OnlinePlotter() = default; //All pointers kept by this class are owned 
                                          //by another object.  

      void ReactEndRun(const std::string& /*fileName*/)
      {
        //Scale all rate histograms here with total elapsed time in run
        const auto deltaT = fRunStopTime-fRunStartTime;
        if(deltaT > 0)
        {
          std::cout << "Elapsed time for this run is " << deltaT << ".  Total elapsed time is " << fStopTotalTime-fStartTotalTime << ".\n";
          const auto timeInv = 1./deltaT;
          fCurrentRunPlots->fMeanRate->Scale(timeInv); 
          fCurrentRunPlots->fMeanRatePerBoard->Scale(timeInv);

          //Fill histograms that profile over board or USB number
          const auto nChannels = fCurrentRunPlots->fMeanRate->GetXaxis()->GetNbins();
          const auto nModules = fCurrentRunPlots->fMeanRate->GetYaxis()->GetNbins();
          for(auto channel = 0; channel < nChannels; ++channel)
          {
            for(auto module = 0; module < nModules; ++module)
            {
              const auto count = fCurrentRunPlots->fMeanRate->GetBinContent(fCurrentRunPlots->fMeanRate->GetBin(channel, module));
              fCurrentRunPlots->fMeanRatePerBoard->Fill(module, count);

              const auto foundUSB = fModuleToUSB.find(module);
              if(foundUSB != fModuleToUSB.end()) 
              {
                const auto usb = foundUSB->second;
                auto found = fUSBToForeverPlots.find(usb);
                if(found == fUSBToForeverPlots.end())
                {
                  found = fUSBToForeverPlots.emplace(usb, fFileService->mkdir("USB"+std::to_string(usb))).first;
                }

                auto& forever = found->second;
                forever.fMeanRate->Fill(fStopTotalTime - fStartTotalTime, count);
                forever.fMeanADC->Fill(fStopTotalTime - fStartTotalTime, 
                                       fCurrentRunPlots->fMeanADC->GetBinContent(fCurrentRunPlots->fMeanADC->GetBin(channel, module)));
              } //If found a USB for this module
            } //For each module in mean rate/ADC plots
          } //For each channel in mean rate/ADC plots
        } //If deltaT > 0

        //Create new directory for next run.  Keep track of run number to do this.  
        ++fRunNum;
      }

      void ReactBeginRun(const std::string& /*fileName*/)
      {
        fCurrentRunPlots.reset(new PerRunPlots(fFileService->mkdir("Run"+std::to_string(fRunNum))));

        fRunStartTime = std::numeric_limits<double>::max();
        fRunStopTime = 0;
      }

      void AnalyzeEvent(const std::vector<CRT::Trigger>& triggers) //Make plots from CRT::Triggers
      {
        for(const auto& trigger: triggers)
        {
          const auto timestamp = trigger.Timestamp();
          if(timestamp > 1e16)
          {
            //Update time bounds based on this timestamp
            const auto unixTime = timestamp >> 32;
            if(unixTime < fRunStartTime) fRunStartTime = unixTime;
            if(unixTime > fRunStopTime) fRunStopTime = unixTime;
            if(unixTime < fStartTotalTime) fStartTotalTime = unixTime;
            if(unixTime > fStopTotalTime) fStopTotalTime = unixTime;

            const auto module = trigger.Channel();
            const auto& hits = trigger.Hits(); 
            double sumADC = 0;
            for(const auto& hit: hits)
            {
              const auto channel = hit.Channel();
              fCurrentRunPlots->fMeanRate->Fill(channel, module); 
              fCurrentRunPlots->fMeanADC->Fill(channel, module, hit.ADC());
              sumADC += hit.ADC();
            }
            fCurrentRunPlots->fMeanADCPerBoard->Fill(module, sumADC/hits.size());
          } //If UNIX timestamp is not 0
        }
      }

    private:
      //Keep track of the current run number
      size_t fRunNum;

      TFS fFileService; //Handle to tool for writing out files
      using DIRECTORY=decltype(fFileService->mkdir("test"));

      //Plots I want for each USB board for each Run
      struct PerRunPlots
      {
        PerRunPlots(DIRECTORY&& dir): fDir(dir)
        {
          fMeanRate = dir.template make<TH2D>("MeanRate", "Mean Rates;Channel;Module;Rate", 64, 0, 64, 32, 0, 32);

          fMeanADC = dir.template make<TProfile2D>("MeanADC", "Mean ADC Values;Channel;Module;Rate [Hz]", 64, 0, 64, 32, 0, 32, 0., 4096);

          fMeanRatePerBoard = dir.template make<TProfile>("MeanRateBoard", "Mean Rate per Board;Board;Rate [Hz]", 
                                                          32, 0, 32);
          fMeanADCPerBoard = dir.template make<TProfile>("MeanADCBoard", "Mean ADC Value per Board;Board;ADC",
                                                         32, 0, 32);
        }

        DIRECTORY fDir; //Directory for the current run where these plots are kept
        TH2D* fMeanRate; //Plot of mean hit rate per channel 
        TProfile2D* fMeanADC; //Plot of mean ADC per channel
        TProfile* fMeanRatePerBoard; //Plot of mean rate for each board
        TProfile* fMeanADCPerBoard; //Plot of mean ADC per board
      };

      std::unique_ptr<PerRunPlots> fCurrentRunPlots; //Per run plots for the current run
      double fRunStartTime; //Earliest UNIX time in run
      double fRunStopTime; //Latest UNIX time in run
      double fStartTotalTime; //Earliest UNIX time stamp in entire job
      double fStopTotalTime; //Lateset UNIX time stamp in entire job

      //Plots I want for each USB board integrated over all(?) Runs
      struct ForeverPlots
      {
        ForeverPlots(DIRECTORY&& dir): fDir(dir)
        {
          fMeanRate = dir.template make<TProfile>("MeanRateHistory", "Mean Rate History;Time [s];Rate [Hz]",
                                                  1000, 0, 3600);
          fMeanADC  = dir.template make<TProfile>("MeanADCHistory", "Mean ADC History;Time [s];ADC", 
                                                  1000, 0, 3600);
        }

        DIRECTORY fDir; //Directory for this USB's overall plots
        TProfile* fMeanRate; //Mean hit rate in time
        TProfile* fMeanADC; //Mean hit ADC in time
      };

      std::unordered_map<unsigned int, ForeverPlots> fUSBToForeverPlots; //Mapping from USB number to ForeverPlots.  

      //Configuration parameters
      std::unordered_map<unsigned int, unsigned int> fModuleToUSB; //Mapping from module number to USB
  };
}

#endif //CRT_ONLINEPLOTTER_CPP
