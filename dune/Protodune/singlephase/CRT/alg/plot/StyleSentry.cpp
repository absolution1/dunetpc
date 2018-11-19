//File: StyleSentry.cpp
//Brief: Allows for changes to TStyle that are local to some scope by making a copy 
//       of the old gStyle for changes and restoring back to the original on destruction.
//       Inspired by how ART's TFileService handles gFile.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//ROOT includes
#include "TStyle.h"

namespace util
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
