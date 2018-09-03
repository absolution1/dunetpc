//File: ChannelView.h
//Brief: Interface for mapping some value in (module, channel) space to 
//       a concrete graphical object that is drawn in a TPad.  Example implementation (derived) 
//       classes might include a simple wrapper over a TH2D and an overlaid 
//       pair of TH2Ds that draw channels as they relate in space.  Derive 
//       from this abstract base class to use it.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRT_CHANNELVIEW_H
#define CRT_CHANNELVIEW_H

//ROOT includes
#include "TPad.h"

//c++ includes
#include <memory>

namespace CRT
{
  class ChannelID;

  class ChannelView
  {
    public:
      ChannelView(); //Create my own TPad for this ChannelView
      ChannelView(TPad* pad); //Use a TPad the user provides for all drawing.  Might be useful if you want to divide a 
                              //TPad for some GUI.

      virtual ~ChannelView();

      //Public interface.  Private implementation below.
      void Fill(const ChannelID& channel);
      void SetValue(const ChannelID& channel, const double value);
      void Draw(const char* option);
      void Reset(const char* option);

    protected:
      //Implement the following methods in a base class to use this interface.
      virtual void doFill(const ChannelID& channel) = 0; //Add a module-channel pair to histogram
      virtual void doSetValue(const ChannelID& channel, const double value) = 0; //Set histogram value for a module-channel 
                                                                                 //pair.  Lets histogram be used as a graph.
      virtual void doDraw(const char* option) = 0; //Queue whatever graphical object(s) are doing the drawing to draw on next update.  
      virtual void doReset(const char* option) = 0; //Reset statistical contents of model that keeps track of bin values

      //Since I am planning to hard-code number of CRTs for now, put it all in one place so that I can maintain that decision in one place
      static constexpr size_t NModules = 32;
      static constexpr size_t ChannelsPerModule = 64;

    private:
      template <class T>
      struct MaybeDeleter
      {
        MaybeDeleter(const bool toDelete): fDeleteMe(toDelete) {}

        void operator ()(T* obj) const
        {
          if(fDeleteMe) delete obj;
        }

        const bool fDeleteMe;
      };

      //WARNING: I might not really own fPad below even though it is a unique_ptr.  In one constructor, it will be 
      //         constructed with the default deleter (and so owned).  In the other constructor, it will be created 
      //         with a no-op deleter (and so not owned).  
      std::unique_ptr<TPad, MaybeDeleter<TPad>> fPad;
  };
}

#endif //CRT_CHANNELVIEW_H
