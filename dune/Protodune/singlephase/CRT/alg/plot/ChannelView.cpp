//File: ChannelView.cpp
//Brief: Implementation of common interface utilities for classes that map (module, channel) 
//       pairs to a 2D histogram on a TPad.  Derive from ChannelView to use it.
//Author: Andrew Olivier aolivier@ur.rochester.edu  

//Include header
#include "plot/ChannelView.h"

//Include some utilities for dealing with CRT::Triggers from Core
#include "util/CRTID.h"

namespace CRT
{
  //I own fPad.  Delete it in destructor.
  ChannelView::ChannelView(): fPad(new TPad("ChannelViewDefault", "CRT::ChannelView Default Pad", 0, 0, 1, 1), MaybeDeleter<TPad>(true))
  {
  }

  //Someone else owns fPad.  Don't delete it in destructor.  
  ChannelView::ChannelView(TPad* pad): fPad(pad, MaybeDeleter<TPad>(false))
  {
  }

  //Implement default destructor so that class templates that want class to only be declared in header will be happy.  
  ChannelView::~ChannelView()
  {
  }

  void ChannelView::Fill(const ChannelID& channel)
  {
    //I could add "utility" processing here so that all derived classes benefit from it
    doFill(channel);
  }

  void ChannelView::SetValue(const ChannelID& channel, const double value)
  {
    //I could add "utility" processing here so that all derived classes benefit from it
    doSetValue(channel, value);
  }

  void ChannelView::Draw(const char* option)
  {
    fPad->cd();
    doDraw(option);
  }
}
