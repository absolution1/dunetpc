//File: SpaceView.cpp
//Brief: A histogramming interface that draws module number on one axis and channel number on the other.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//Include header
#include "plot/SpaceView.h"

//ROOT includes
#include "TStyle.h"

//Local includes
#include "StyleSentry.cpp"

//c++ includes
#include <iostream>

namespace CRT
{
  void SpaceView::ConfigHistogram(TH2& hist)
  {
    hist.SetStats(false);
    //auto zAxis = hist.GetZaxis();
    //TODO: x and y labels based on approximate module positions?
  }

  void SpaceView::SetupPads()
  {
    auto main = GetMainPad(); //Observer pointer
    main->Divide(2);

    //Set up upstream pad
    main->cd(1);
    gPad->SetLogz();

    //Set up downstream pad
    main->cd(2);
    gPad->SetLogz();
  }

  SpaceView::SpaceView(const std::string& name, const std::string& title, const std::string& zTitle): ChannelView()
  {
    const auto binLimit = 4*ChannelsPerModule*2; //Frames are 4 modules high.  2 strips per module to simulate offsets of "bottom" strip layer
    fUpstream = TH2D(("upstream"+name).c_str(), ("Upstream;cartoon x;cartoon y;"+zTitle).c_str(), binLimit, 0, binLimit, binLimit, 0, binLimit);
    ConfigHistogram(fUpstream);

    fDownstream = TH2D(("downstream"+name).c_str(), ("Downstream;cartoon x;cartoon y;"+zTitle).c_str(), binLimit, 0, binLimit, binLimit, 0, binLimit);
    ConfigHistogram(fDownstream);

    SetupPads();
  }
   
  SpaceView::SpaceView(TPad* pad, const std::string& name, const std::string& title, const std::string& zTitle): ChannelView(pad)
  {
    const auto binLimit = 4*ChannelsPerModule*2;
    std::cout << "binLimit is " << binLimit << " = 4*" << ChannelsPerModule << "*2\n";
    fUpstream = TH2D(("upstream"+name).c_str(), ("Upstream;cartoon x;cartoon y;"+zTitle).c_str(), binLimit, 0, binLimit, binLimit, 0, binLimit);
    fDownstream = TH2D(("downstream"+name).c_str(), ("Downstream;cartoon x;cartoon y;"+zTitle).c_str(), binLimit, 0, binLimit, binLimit, 0, binLimit);

    ConfigHistogram(fUpstream);
    ConfigHistogram(fDownstream);

    SetupPads();
  }

  SpaceView::SpaceView(const double zMax, const std::string& name, const std::string& title, const std::string& zTitle): SpaceView(name, title, zTitle)
  {
    fUpstream.SetMaximum(zMax);
    fDownstream.SetMaximum(zMax);
  }

  SpaceView::SpaceView(TPad* pad, const double zMax, const std::string& name, const std::string& title, const std::string& zTitle): SpaceView(pad, name, title, zTitle)
  {
    fUpstream.SetMaximum(zMax);
    fDownstream.SetMaximum(zMax);
  }

  void SpaceView::doFill(const size_t module, const size_t channel, const double weight)
  {
    doSomething((module < NModules/2l)?fUpstream:fDownstream, module, channel, 
                [&weight](auto& hist, const auto x, const auto y) { hist.Fill(x, y, weight); });
  }

  template <class FUNC>
  void SpaceView::doSomething(TH2& hist, const size_t module, const size_t channel, FUNC&& func)
  {
    //When doing things like this in the future, a frame would be a good logical abstraction layer.  
    const size_t local = module%16;
    const int nXBins = hist.GetXaxis()->GetNbins(), nYBins = hist.GetYaxis()->GetNbins(), endOfFrame = 2l*ChannelsPerModule*2l; 
    //Frames are 2 modules x 2 modules.  2 columns per strip
    const bool secondLayer = (channel > 32);

    std::cout << "module is " << module << ", channel is " << channel << ", local is " << local << ".\n";
    if(local < 2) //Top beam-left
    {
      const auto xbin = endOfFrame-(local*ChannelsPerModule+channel)*2+secondLayer; 
      std::cout << "Filling beam-left with subtraction at xbin=" << xbin << "\n";
      for(int ybin = endOfFrame; ybin < nYBins; ++ybin) 
      {
        //Fill both strips that this channel overlaps
        func(hist, xbin, ybin);
        func(hist, xbin-1, ybin);
      }
    }
    else if(local > 13) //Top beam-right
    {
      const auto xbin = nXBins-((local-14)*ChannelsPerModule+channel)*2+secondLayer; 
      std::cout << "Filling beam-right with subtraction at xbin=" << xbin << "\n";
      for(int ybin = endOfFrame; ybin < nYBins; ++ybin) 
      {
        //Fill both strips that this channel overlaps
        func(hist, xbin, ybin);
        func(hist, xbin-1, ybin);
      }
    }
    else if(local > 5 && local < 10) //Bottom vertical module
    {
      const auto xbin = ((local-6)*ChannelsPerModule+channel)*2-secondLayer; //Entire face is 4 modules x 4 modules
      std::cout << "Filling at xbin=" << xbin << "\n";
      for(int ybin = 0; ybin < endOfFrame; ++ybin) 
      {
        //Fill both strips that this channel overlaps
        func(hist, xbin, ybin);
        func(hist, xbin+1, ybin);
      }
    }
    else if(local < 6) //Beam-left horizontal module
    {
      const auto ybin = nYBins-((local-2)*ChannelsPerModule+channel)*2+secondLayer; //Entire face is 4 modules x 4 modules
      std::cout << "Filling at ybin=" << ybin << "\n";
      for(int xbin = 0; xbin < endOfFrame; ++xbin) 
      {
        //Fill both strips that this channel overlaps
        func(hist, xbin, ybin);
        func(hist, xbin, ybin-1);
      }
    }
    else //Beam-right horizontal module
    {
      const auto ybin = ((local-10)*ChannelsPerModule+channel)*2-secondLayer; //Entire face is 4 modules x 4 modules
      std::cout << "Filling at ybin=" << ybin << "\n";
      for(int xbin = endOfFrame; xbin < nXBins; ++xbin) 
      {
        //Fill both strips that this channel overlaps
        func(hist, xbin, ybin);
        func(hist, xbin, ybin+1);
      }
    }
  }

  void SpaceView::doSetValue(const size_t module, const size_t channel, const double value)
  {
    doSomething((channel < 32)?fUpstream:fDownstream, module, channel, 
                [&value](auto& hist, const auto x, const auto y) { hist.SetBinContent(x, y, value); });
  }

  void SpaceView::doDraw(const char* option)
  {
    util::StyleSentry old; //Save the old style during this function so I don't change it

    gStyle->SetOptStat(0);
    fUpstream.UseCurrentStyle();
    fDownstream.UseCurrentStyle();

    auto main = GetMainPad();
    main->cd(1);
    fUpstream.Draw(option);

    main->cd(2);
    fDownstream.Draw(option);
  }

  void SpaceView::doReset(const char* option)
  {
    fUpstream.Reset(option);
    fDownstream.Reset(option);
  }
}
