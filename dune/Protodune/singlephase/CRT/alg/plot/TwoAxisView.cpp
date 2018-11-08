//File: TwoAxisView.cpp
//Brief: A histogramming interface that draws module number on one axis and channel number on the other.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//Include header
#include "plot/TwoAxisView.h"

//ROOT includes
#include "TStyle.h"

namespace 
{
  struct StyleSentry
  {
    StyleSentry(): fOldStyle(gStyle) 
    {
      if(gStyle) gStyle = new TStyle(*gStyle);
      else gStyle = new TStyle();
    }

    ~StyleSentry() { gStyle = fOldStyle; }

    private:
     TStyle* fOldStyle; //The old gStyle that will be restored when this object is destroyed
  };
}

namespace CRT
{
  TwoAxisView::TwoAxisView(const std::string& name, const std::string& title, const std::string& xTitle, const std::string& yTitle, 
                           const std::string& zTitle): ChannelView()
  {
    fHist = TH2D(name.c_str(), (title+";"+xTitle+";"+yTitle+";"+zTitle).c_str(), ChannelsPerModule, 0, ChannelsPerModule, NModules, 0, NModules);
    fHist.SetStats(false);
    fHist.SetMaximum(4096); //Set range to hard-coded maximum ADC value in CRT hardware
  }
   
  TwoAxisView::TwoAxisView(const std::string& name, const std::string& title, const std::string& xTitle, const std::string& yTitle,
                           const std::string& zTitle, TPad* pad): ChannelView(pad)
  {
    fHist = TH2D(name.c_str(), (title+";"+xTitle+";"+yTitle+";"+zTitle).c_str(), ChannelsPerModule, 0, ChannelsPerModule, NModules, 0, NModules);
  }

  TwoAxisView::TwoAxisView(const std::string& name, const std::string& title, const std::string& zTitle): TwoAxisView(name, title, "channel", 
                                                                                                                      "module", zTitle)
  {
  }

  TwoAxisView::TwoAxisView(TPad* pad, const std::string& name, const std::string& title, const std::string& zTitle): 
                           TwoAxisView(name, title, "channel", "module", zTitle, pad)
  {
  }

  void TwoAxisView::doFill(const size_t module, const size_t channel, const double weight)
  {
    fHist.Fill(channel, module, weight);
  }

  void TwoAxisView::doSetValue(const size_t module, const size_t channel, const double value)
  {
    fHist.SetBinContent(fHist.GetBin(channel, module), value);
  }

  void TwoAxisView::doDraw(const char* option)
  {
    StyleSentry old; //Save the old style during this function so I don't change it

    gStyle->SetOptStat(0);
    fHist.UseCurrentStyle();
    fHist.Draw(option);
  }

  void TwoAxisView::doReset(const char* option)
  {
    fHist.Reset(option);
  }
}
