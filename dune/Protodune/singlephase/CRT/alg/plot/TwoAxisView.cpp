//File: TwoAxisView.cpp
//Brief: A histogramming interface that draws module number on one axis and channel number on the other.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//Include header
#include "plot/TwoAxisView.h"

//Include CRTID definitions for doFill() and doSetValue()
#include "util/CRTID.h"

namespace CRT
{
  TwoAxisView::TwoAxisView(const std::string& name, const std::string& title, const std::string& xTitle, const std::string& yTitle, 
                           const std::string& zTitle): ChannelView()
  {
    fHist = TH2D(name.c_str(), (title+";"+xTitle+";"+yTitle+";"+zTitle).c_str(), ChannelsPerModule, 0, ChannelsPerModule, NModules, 0, NModules);
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

  TwoAxisView::TwoAxisView(const std::string& name, const std::string& title, const std::string& zTitle, TPad* pad): 
                           TwoAxisView(name, title, "channel", "module", zTitle, pad)
  {
  }

  void TwoAxisView::doFill(const ChannelID& channel)
  {
    fHist.Fill(channel.ChannelNum, channel.ModuleNum);
  }

  void TwoAxisView::doSetValue(const ChannelID& channel, const double value)
  {
    fHist.SetBinContent(fHist.GetBin(channel.ChannelNum, channel.ModuleNum), value);
  }

  void TwoAxisView::doDraw(const char* option)
  {
    fHist.Draw(option);
  }

  void TwoAxisView::doReset(const char* option)
  {
    fHist.Reset(option);
  }
}
