//File: SpaceView.h
//Brief: A SpaceView is a ChannelView (mapping from (module, channel) pairs to histogram bins)
//       that draws x coordinate on one axis and y coordinate on another.  
//Author: Andrew Olivier aolivier@ur.rochester.edu  

//Include header
#include "plot/ChannelView.h"

//ROOT includes
#include "TH2D.h"

namespace CRT
{
  class SpaceView: public ChannelView
  {
    public:
      //Axis labels are automatically generated for x and y axes
      SpaceView(const std::string& name="CRTEvd", const std::string& title="CRT Event Display", const std::string& zTitle="Hits");
      SpaceView(TPad* pad, const std::string& name="CRTEvd", const std::string& title="CRT Event Display", const std::string& zTitle="Hits");

      //Automatic axis labels with maximum z value.  Useful for ADCs
      SpaceView(const double zMax, const std::string& name="CRTEvd", const std::string& title="CRT Event Display", const std::string& zTitle="Hits");
      SpaceView(TPad* pad, const double zMax, const std::string& name="CRTEvd", const std::string& title="CRT Event Display", const std::string& zTitle="Hits");

      //Public interface is defined in base class
    protected:
      virtual void doFill(const size_t module, const size_t channel, const double weight) override; 
      virtual void doSetValue(const size_t module, const size_t channel, const double value) override;
      virtual void doDraw(const char* option) override;
      virtual void doReset(const char* option) override;

    private:
      TH2D fUpstream; //Cartoon of upstream hit positions
      TH2D fDownstream; //Cartoon of downstream hit positions
      
      //Concentrate code for deciding which histogram to fill in one place
      template <class FUNC> //FUNC is a callable object that takes a TH2& and two bin numbers and uses them to 
                            //update the values in a histogram.  
      void doSomething(TH2& hist, const size_t channel, const size_t module, FUNC&& func); 
      void ConfigHistogram(TH2& hist); 
      void SetupPads();  
  };
}
