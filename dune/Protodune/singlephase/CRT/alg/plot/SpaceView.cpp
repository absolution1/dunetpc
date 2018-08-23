//File: SpaceView.cpp
//Brief: Implementation of a histogram-like interface for plotting 
//       CRT information in the relative spatial positions of the CRT 
//       modules.  Implements ChannelView. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TBox.h"

//Include header
#include "SpaceView.h"

namespace CRT
{
  SpaceView::SpaceView(): ChannelView(), fModuleToBox()
  {
    BuildGeometry();
  }

  SpaceView::SpaceView(TPad* pad): ChannelView(pad), fModuleToBox()
  {
    BuildGeometry();
  }

  //Create the TBoxes in a CRT module given that module's center x and y position
  void SpaceView::BuildModule(const double x, const double y, const bool vertical)
  {
    std::vector<Strip> module;
    double moduleWidth, moduleHeight;
    if(vertical) 
    {
      moduleWidth = ChannelsPerModule*fStripWidth;
      moduleHeight = ChannelsPerModule*fStripLength;
    }
    else
    {
      moduleWidth = ChannelsPerModule*fStripLength;
      moduleHeight = ChannelsPerModule*fStripWidth;
    }
     
    bool upperLayer = true; //Each module consists of two layers of strips.  The lower layer of strips is offset 
                            //half of a strip width away from strip 32.  
    for(size_t strip = 0; strip < ChannelsPerModule; ++strip)
    {
      if(strip > ChannelsPerModule/2) upperLayer = true;
      Strip sens;
      sens.fValue = 0;
      sens.fBox = TBox((strip-0.5)*fStripWidth - moduleWidth/2. - x - upperLayer*fStripWidth/2., moduleHeight/2. - y, 
                       (strip+0.5)*fStripWidth - x - fStripWidth/2., -moduleHeight/2. - y);
      //TODO: Box outlines, probably in black
      module.push_back(sens);
    }

    fModuleToBox.push_back(module);
  }

  //Create TBoxes grouped by channel number according to how the CRT is layed out.  
  //TODO: If I ever write a geometry interface for the CRT, replace most of what this 
  //      function does with that geometry interface.  
  void SpaceView::BuildGeometry()
  {
    //Create TBoxes for the upsteam modules according to https://www.dropbox.com/home/Protodune/CRT%20cables?preview=Protodune_CRT_Mapping.pptx#
    const double moduleWidth = ChannelsPerModule*fStripWidth, moduleHeight = ChannelsPerModule*fStripLength;
    //TODO: Use FrameWidth/Height to make things line up better.  
    
    const double offset = (2.+0.1)*moduleWidth; //Offset upstream by downstream by this amount so that there will be space between them
    //Build upstream beam-left (Jura)
    //TODO: Upstream Label
    //Top
    BuildModule(-moduleWidth*1./2.-offset, moduleHeight*1./2., true);
    BuildModule(-moduleWidth*3./2.-offset, moduleHeight*1./2., true);
    BuildModule(-moduleHeight*1./2.-offset, moduleWidth*3./2., false);
    BuildModule(-moduleHeight*1./2.-offset, moduleWidth*1./2., false);

    //Bottom
    BuildModule(-moduleHeight*1./2.-offset, -moduleWidth*1./2., false);
    BuildModule(-moduleHeight*1./2.-offset, -moduleWidth*3./2., false);
    BuildModule(-moduleWidth*3./2.-offset, -moduleHeight*1./2., true);
    BuildModule(-moduleWidth*1./2.-offset, -moduleHeight*1./2., true);

    //TODO: Make room for the beam pipe
    //Build upstream beam-right (Saleve)
    //Bottom
    BuildModule(moduleWidth*1./2.-offset, -moduleHeight*1./2., true);
    BuildModule(moduleWidth*3./2.-offset, -moduleHeight*1./2., true);
    BuildModule(moduleHeight*1./2.-offset, -moduleWidth*3./2., false);
    BuildModule(moduleheight*1./2.-offset, -moduleWidth*1./2., false);
    
    //Top
    BuildModule(moduleHeight*1./2.-offset, moduleWidth*1./2., false);
    BuildModule(moduleHeight*1./2.-offset, moduleWidth*3./2., false);
    BuildModule(moduleWidth*3./2.-offset, moduleHeight*1./2., true);
    BuildModule(moduleWidth*1./2.-offset, moduleHeight*1./2., true);

    //Build downstream beam-left (Jura)
    //TODO: Downstream label
    //Top
    BuildModule(-moduleWidth*1./2.+offset, moduleHeight*1./2., true);
    BuildModule(-moduleWidth*3./2.+offset, moduleHeight*1./2., true);
    BuildModule(-moduleHeight*1./2.+offset, moduleWidth*3./2., false);
    BuildModule(-moduleHeight*1./2.+offset, moduleWidth*1./2., false);
                                                                 
    //Bottom                                                        
    BuildModule(-moduleHeight*1./2.+offset, -moduleWidth*1./2., false);
    BuildModule(-moduleHeight*1./2.+offset, -moduleWidth*3./2., false);
    BuildModule(-moduleWidth*3./2.+offset, -moduleHeight*1./2., true);
    BuildModule(-moduleWidth*1./2.+offset, -moduleHeight*1./2., true);
                                                                 
    //Build downstream beam-right (Saleve)
    //Bottom
    BuildModule(moduleWidth*1./2.+offset, -moduleHeight*1./2., true);
    BuildModule(moduleWidth*3./2.+offset, -moduleHeight*1./2., true);
    BuildModule(moduleHeight*1./2.+offset, -moduleWidth*3./2., false);
    BuildModule(moduleheight*1./2.+offset, -moduleWidth*1./2., false);
    
    //Top
    BuildModule(moduleHeight*1./2.+offset, moduleWidth*1./2., false);
    BuildModule(moduleHeight*1./2.+offset, moduleWidth*3./2., false);
    BuildModule(moduleWidth*3./2.+offset, moduleHeight*1./2., true);
    BuildModule(moduleWidth*1./2.+offset, moduleHeight*1./2., true);
  }

  void SpaceView::doFill(const ChannelID& channel)
  {
    ++fModuleToBox[channel.ModuleNum][channel.ChannelNum].fValue;
  }

  void SpaceView::doSetValue(const ChannelID& channel, const double value)
  {
    fModuleToBox[channel.ModuleNum][channel.ChannelNum].fValue = value;
  }

  void SpaceView::doDraw(const char* /*option*/)
  {
    for(size_t module = 0; module < fModuleToBox.size()/2.; ++module)
    {
      for(auto& channel: module)
      {
        auto& box = channel.fBox;
        box.SetFillColorAlpha(fPalette.GetValueColor(box.fValue), 0.25); //Hard-code alpha value here so that 4 overlapping 
                                                                         //channels give an alpha value of 1.
        box.Draw();
      }
    }

    fPalette.Draw();
  }

  void SpaceView::doReset(const char* /*option*/)
  {
    for(auto& module: fModuleToBox)
    {
      for(auto& channel: module) channel.fValue = 0.;
    }
  }
} 
