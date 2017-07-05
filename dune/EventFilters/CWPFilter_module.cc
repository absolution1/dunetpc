#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "TMath.h"
//  Filter that uses counter information to select tracks parallel to the 
//    wires planes and passing close to them.  
#include "lardataobj/RawData/ExternalTrigger.h"
#include <fstream>

namespace filt{

 class CWPFilter : public art::EDFilter {
   public:
    explicit      CWPFilter(fhicl::ParameterSet const& pset);
    virtual      ~CWPFilter() { }
    virtual bool  filter(art::Event& e);
    void          reconfigure(fhicl::ParameterSet const& pset);

  private:

   //ofstream         which_Counters;
   int              coincidenceWindow = 5;
   std::string      fCounterModuleLabel;
   std::vector<int> fFirstSet;
   std::vector<int> fSecondSet;
  };

   CWPFilter::CWPFilter(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  void CWPFilter::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCounterModuleLabel = pset.get< std::string >("CounterModuleLabel");
    fFirstSet           = pset.get<std::vector<int>>("firstset");
    fSecondSet          = pset.get<std::vector<int>>("secondset");
  }
  
  bool CWPFilter::filter(art::Event& e)
  {
    bool keepFlag=1;

    int event = e.id().event();
    //int run   = e.run();

    std::cout << " event " << event << std::endl;

    // 2016Feb11 - MStancari
    // Look for counter pairs auxdetid=(9,8,7) (29,30,31) 
    //   93,94,95 can be used in addition to 8,7, but may not yet be part 
    //   of the EL-WU trigger

    art::Handle< std::vector< raw::ExternalTrigger> > externalTriggerListHandle;
    e.getByLabel(fCounterModuleLabel, externalTriggerListHandle);
    std::vector< art::Ptr< raw::ExternalTrigger> > trigs;
    art::fill_ptr_vector(trigs,externalTriggerListHandle);

    std::vector<std::vector<std::vector<long long>>> kept_hits(fFirstSet.size());

    for(unsigned int i = 0; i < fFirstSet.size(); i++)
    {
      kept_hits.at(i).resize(2);
    }
   
    unsigned int nchits = trigs.size();

    for(unsigned int i = 0; i < fFirstSet.size(); i++)
    {
      for(unsigned int j = 0; j < nchits; j++)
      {
        int       auxdetid  = trigs.at(j)->GetTrigID();
        long long trig_time = trigs.at(j)->GetTrigTime();

        if(auxdetid==fFirstSet.at(i))
        {
          kept_hits.at(i).at(0).push_back(trig_time);
          std::cout << "Counter hit ID: " << auxdetid << " Trig time: " << trig_time << std::endl; 
        }
        else if(auxdetid==fSecondSet.at(i))
        {
          kept_hits.at(i).at(1).push_back(trig_time);
          std::cout << "Counter hit ID: " << auxdetid << " Trig time: " << trig_time << std::endl; 
        }
      }
    }

    for(unsigned int i = 0; i < fFirstSet.size(); i++)
    {
      bool break_signal = 0;
      if(kept_hits.at(i).at(0).size()>0&&kept_hits.at(i).at(1).size()>0)
      {
        for(unsigned int j = 0; j < kept_hits.at(i).at(0).size(); j++)
        {
          for(unsigned int k = 0; k < kept_hits.at(i).at(1).size(); k++)
          {
            if(kept_hits.at(i).at(0).size()>=kept_hits.at(i).at(1).size())
            {
              if(k<=j)
              {
                if(std::abs(kept_hits.at(i).at(0).at(j)-kept_hits.at(i).at(1).at(k))<coincidenceWindow)
                {
                  std::cout << "Coincident IDs : "    << fFirstSet.at(i) << " & " << fSecondSet.at(i)
                            << " Concidence window: " << std::abs(kept_hits.at(i).at(0).at(j)-kept_hits.at(i).at(1).at(k)) << std::endl;

                  /*
                  which_Counters.open("/dune/data/users/abooth/whichCounters.txt", std::ios_base::app);
                  which_Counters << "Event: "       << event  << " | Run: "       << run << " | ";
                  which_Counters << fFirstSet.at(i) << "\t\t" << fSecondSet.at(i) << std::endl;
                  which_Counters.close();
                  */

                  keepFlag     = 1;
                  break_signal = 1;
                  break;
                }
                else
                {
                  keepFlag = 0;
                }
              }
            }
            else
            {
              if(j<=k)
              {
                if(std::abs(kept_hits.at(i).at(0).at(j)-kept_hits.at(i).at(1).at(k))<coincidenceWindow)
                {
                  std::cout << "Coincident IDs : "     << fFirstSet.at(i) << " & " << fSecondSet.at(i)
                            << " Coincidence Window: " << std::abs(kept_hits.at(i).at(0).at(j)-kept_hits.at(i).at(1).at(k)) << std::endl; 

                  /*
                  which_Counters.open("/dune/data/users/abooth/whichCounters.txt", std::ios_base::app);
                  which_Counters << "Event: "       << event  << " | Run: "       << run << " | ";
                  which_Counters << fFirstSet.at(i) << "\t\t" << fSecondSet.at(i) << std::endl;
                  which_Counters.close();
                  */

                  keepFlag     = 1;
                  break_signal = 1;
                  break;
                }
                else
                {
                  keepFlag = 0;
                }
              }
            }
          }
          if(break_signal == 1)
          {
            break;
          }
        }
      }
      else
      {
        keepFlag = 0;
      }

      if(break_signal == 1)
      {
        break; 
      }
    }

    std::cout << "event=" << event << " kept=" << keepFlag << std::endl;

    return keepFlag;
  }
  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(CWPFilter)

} // namespace filt


