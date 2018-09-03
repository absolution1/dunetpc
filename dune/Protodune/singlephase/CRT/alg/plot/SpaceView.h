//File: SpaceView.h
//Brief: Histograms CRT modules as a cartoon of their relative 
//       positions in space as seen from the incoming beam's point of view. 
//       Implements ChannelView. 
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRT_SPACEVIEW_H
#define CRT_SPACEVIEW_H

//Base class header
#include "plot/ChannelView.h"

namespace CRT
{
  class SpaceView: public ChannelView
  {
    public:
      SpaceView(); //Default constructor provides its' own TPad for drawing
      SpaceView(TPad* pad); //Let the user provide a TPad for drawing.  Could be useful with 
                            //something like TFileService to write the TPad out to a file later.

      virtual ~SpaceView(); //Trivial destructor, but defer it to the implementation file in 
                            //case I want to use something like Factory later.

      //The rest of the public interface is defined by the ChannelView base class

    protected:
      //Implement required member functions from ChannelView
      virtual void doFill(const ChannelID& channel) override; //Add a module-channel pair to histogram
      virtual void doSetValue(const ChannelID& channel, const double value) override; //Set histogram value for a module-channel 
                                                                                      //pair.  Lets histogram be used as a graph.
      virtual void doDraw(const char* option) override; //Queue whatever graphical object(s) are doing the drawing to draw on next update.
      virtual void doReset(const char* option) override; //Reset contents of bin value model to 0

    private:
      //Group a graphics object together with the bin value it represents
      struct Strip
      {
        TBox fBox; //Graphics object that visualizes data
        double fBinValue; //Value in this bin.  Converted to fBox's color.
      };

      std::vector<std::vector<Strip>> fModuleToBox; //Mapping from (module identifier, strip identifier) to graphics object-bin value pair
      TPaletteAxis fPalette; //Visualization and model for mapping a value to a color in ROOT's current palette.
      //TODO: Some circle primitive for the beam pipe?

      //Geometry data used to put channels in the right position
      static constexpr double fStripWidth = 5.; //Width of a CRT scintillator strip in cm
      static constexpr double fStripLength = 320.; //Length of longest dimension of a CRT strip in cm
      
      void BuildGeometry(); //Build the map from ChannelID to TBox using geometry information about CRT detector
      void BuildModule(const double x, const double y, const bool vertical); //Create the Strips in one CRT module using its' center x and y coordinates
  };
}

#endif //CRT_SPACEVIEW_H  
