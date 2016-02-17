#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "TMath.h"
//  Filter that uses counter information to select tracks parallel to the 
//    wires planes and passing close to them.  
#include "RawData/ExternalTrigger.h"

namespace filt{

 class CWPFilter : public art::EDFilter {
   public:
     explicit CWPFilter(fhicl::ParameterSet const& pset);
    virtual ~CWPFilter() { }
    virtual bool filter(art::Event& e);
    void    reconfigure(fhicl::ParameterSet const& pset);

  private:

   int coincidenceWindow = 5;
   std::string fCounterModuleLabel;

  };

   CWPFilter::CWPFilter(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

  }

  void CWPFilter::reconfigure(fhicl::ParameterSet const& pset)
  {
  fCounterModuleLabel = pset.get< std::string >("CounterModuleLabel");
      
  }
  
  bool CWPFilter::filter(art::Event& e)
  {


    // set keepFlag=0 to reject the event
    bool keepFlag=1;

    // Get event number from the event.
    int event = e.id().event();

    std::cout << " event " << event << std::endl;

    // 2016Feb11 - MStancari
    // Look for counter pairs auxdetid=(9,8,7) (29,30,31) 
    //   93,94,95 can be used in addition to 8,7, but may not yet be part 
    //   of the EL-WU trigger

    art::Handle< std::vector< raw::ExternalTrigger> > externalTriggerListHandle;
    e.getByLabel(fCounterModuleLabel, externalTriggerListHandle);
    std::vector< art::Ptr< raw::ExternalTrigger> > trigs;
    art::fill_ptr_vector(trigs,externalTriggerListHandle);

    std::vector<std::vector<int>> C6_31(2);
    std::vector<std::vector<int>> C7_30(2);
    std::vector<std::vector<int>> C8_29(2);
    std::vector<std::vector<int>> C9_28(2);
    std::vector<std::vector<int>> C6_32(2);
    std::vector<std::vector<int>> C10_28(2);
    std::vector<std::vector<int>> C16_43(2);
    std::vector<std::vector<int>> C21_38(2);
    std::vector<std::vector<int>> C0_27(2);
    std::vector<std::vector<int>> C5_22(2);

    unsigned int nchits = trigs.size();

    std::cout << "Number of counter hits: " << nchits << std::endl;

    for(unsigned int i = 0; i < nchits; i++)
    {
      int auxdetid = trigs.at(i)->GetTrigID();
      
      if(auxdetid==6)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C6_31.at(0).push_back(tick);
      }
      else if(auxdetid==31)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C6_31.at(1).push_back(tick);
      }
      else if(auxdetid==7)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C7_30.at(0).push_back(tick);
      }
      else if(auxdetid==30)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C7_30.at(1).push_back(tick);
      }
      else if(auxdetid==8)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C8_29.at(0).push_back(tick);
      }
      else if(auxdetid==29)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C8_29.at(1).push_back(tick);
      }
      else if(auxdetid==9)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C9_28.at(0).push_back(tick);
      }
      else if(auxdetid==28)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C9_28.at(1).push_back(tick);
      }
      else if(auxdetid==6)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C6_32.at(0).push_back(tick);
      }
      else if(auxdetid==32)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C6_32.at(1).push_back(tick);
      }
      else if(auxdetid==10)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C10_28.at(0).push_back(tick);
      }
      else if(auxdetid==28)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C10_28.at(1).push_back(tick);
      }
      else if(auxdetid==16)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C16_43.at(0).push_back(tick);
      }
      else if(auxdetid==43)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C16_43.at(1).push_back(tick);
      }
      else if(auxdetid==21)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C21_38.at(0).push_back(tick);
      }
      else if(auxdetid==38)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C21_38.at(1).push_back(tick);
      }
      else if(auxdetid==0)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C0_27.at(0).push_back(tick);
      }
      else if(auxdetid==27)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C0_27.at(1).push_back(tick);
      }
      else if(auxdetid==5)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C5_22.at(0).push_back(tick);
      }
      else if(auxdetid==22)
      {
        int tick = trigs.at(i)->GetTrigTime();
        C5_22.at(1).push_back(tick);
      }
    }

    std::cout << "EL10_1: " << C6_31.at(0).size()  << std::endl;
    std::cout << "WU4: "    << C6_31.at(1).size()  << std::endl;
    std::cout << "EL9: "    << C7_30.at(0).size()  << std::endl;
    std::cout << "WU3: "    << C7_30.at(1).size()  << std::endl;
    std::cout << "EL8: "    << C8_29.at(0).size()  << std::endl;
    std::cout << "WU2: "    << C8_29.at(1).size()  << std::endl;
    std::cout << "EL7: "    << C9_28.at(0).size() << std::endl;
    std::cout << "WU1_1: "  << C9_28.at(1).size() << std::endl;
    std::cout << "EL10_2: " << C6_32.at(0).size() << std::endl;
    std::cout << "WU5: "    << C6_32.at(1).size() << std::endl;
    std::cout << "EL6: "    << C10_28.at(0).size()  << std::endl;
    std::cout << "WU1_2: "  << C10_28.at(1).size()  << std::endl;
    std::cout << "NL1: "    << C16_43.at(0).size()  << std::endl;
    std::cout << "SU1: "    << C16_43.at(1).size()  << std::endl;
    std::cout << "NL6: "    << C21_38.at(0).size()  << std::endl;
    std::cout << "SU6: "    << C21_38.at(1).size()  << std::endl;
    std::cout << "SL1: "    << C0_27.at(0).size()  << std::endl;
    std::cout << "NU6: "    << C0_27.at(1).size()  << std::endl;
    std::cout << "SL6: "    << C5_22.at(0).size()  << std::endl;
    std::cout << "NU1: "    << C5_22.at(1).size()  << std::endl;

    if(C6_31.at(0).size()>0&&C6_31.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C6_31.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C6_31.at(1).size(); j++)
        {
          if(C6_31.at(0).size()>=C6_31.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C6_31.at(0).at(i)-C6_31.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL10 & WU4, difference in trigger times: " << std::abs(C6_31.at(0).at(i)-C6_31.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C6_31.at(0).at(i)-C6_31.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL10 & WU4, difference in trigger times: " << std::abs(C6_31.at(0).at(i)-C6_31.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C7_30.at(0).size()>0&&C7_30.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C7_30.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C7_30.at(1).size(); j++)
        {
          if(C7_30.at(0).size()>=C7_30.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C7_30.at(0).at(i)-C7_30.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL9 & WU3, difference in trigger times: " << std::abs(C7_30.at(0).at(i)-C7_30.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C7_30.at(0).at(i)-C7_30.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL9 & WU3, difference in trigger times: " << std::abs(C7_30.at(0).at(i)-C7_30.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C8_29.at(0).size()>0&&C8_29.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C8_29.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C8_29.at(1).size(); j++)
        {
          if(C8_29.at(0).size()>=C8_29.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C8_29.at(0).at(i)-C8_29.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL8 & WU2, difference in trigger times: " << std::abs(C8_29.at(0).at(i)-C8_29.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C8_29.at(0).at(i)-C8_29.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL8 & WU2, difference in trigger times: " << std::abs(C8_29.at(0).at(i)-C8_29.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C9_28.at(0).size()>0&&C9_28.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C9_28.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C9_28.at(1).size(); j++)
        {
          if(C9_28.at(0).size()>=C9_28.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C9_28.at(0).at(i)-C9_28.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL7 & WU1, difference in trigger times: " << std::abs(C9_28.at(0).at(i)-C9_28.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C9_28.at(0).at(i)-C9_28.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL7 & WU1, difference in trigger times: " << std::abs(C9_28.at(0).at(i)-C9_28.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C6_32.at(0).size()>0&&C6_32.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C6_32.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C6_32.at(1).size(); j++)
        {
          if(C6_32.at(0).size()>=C6_32.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C6_32.at(0).at(i)-C6_32.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL10 & WU5, difference in trigger times: " << std::abs(C6_32.at(0).at(i)-C6_32.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C6_32.at(0).at(i)-C6_32.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL10 & WU5, difference in trigger times: " << std::abs(C6_32.at(0).at(i)-C6_32.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C10_28.at(0).size()>0&&C10_28.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C10_28.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C10_28.at(1).size(); j++)
        {
          if(C10_28.at(0).size()>=C10_28.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C10_28.at(0).at(i)-C10_28.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL6 & WU1, difference in trigger times: " << std::abs(C10_28.at(0).at(i)-C10_28.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C10_28.at(0).at(i)-C10_28.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between EL6 & WU1, difference in trigger times: " << std::abs(C10_28.at(0).at(i)-C10_28.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C16_43.at(0).size()>0&&C16_43.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C16_43.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C16_43.at(1).size(); j++)
        {
          if(C16_43.at(0).size()>=C16_43.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C16_43.at(0).at(i)-C16_43.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between NL1 & SU1, difference in trigger times: " << std::abs(C16_43.at(0).at(i)-C16_43.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C16_43.at(0).at(i)-C16_43.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between NL1 & SU1, difference in trigger times: " << std::abs(C16_43.at(0).at(i)-C16_43.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C21_38.at(0).size()>0&&C21_38.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C21_38.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C21_38.at(1).size(); j++)
        {
          if(C21_38.at(0).size()>=C21_38.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C21_38.at(0).at(i)-C21_38.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between NL6 & SU6, difference in trigger times: " << std::abs(C21_38.at(0).at(i)-C21_38.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C21_38.at(0).at(i)-C21_38.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between NL6 & SU6, difference in trigger times: " << std::abs(C21_38.at(0).at(i)-C21_38.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C0_27.at(0).size()>0&&C0_27.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C0_27.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C0_27.at(1).size(); j++)
        {
          if(C0_27.at(0).size()>=C0_27.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C0_27.at(0).at(i)-C0_27.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between SL1 & NU6, difference in trigger times: " << std::abs(C0_27.at(0).at(i)-C0_27.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C0_27.at(0).at(i)-C0_27.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between SL1 & NU6, difference in trigger times: " << std::abs(C0_27.at(0).at(i)-C0_27.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else if(C5_22.at(0).size()>0&&C5_22.at(1).size()>0)
    {
      for(unsigned int i = 0; i < C5_22.at(0).size(); i++)
      {
        bool break_Signal = 0;
        for(unsigned int j = 0; j < C5_22.at(1).size(); j++)
        {
          if(C5_22.at(0).size()>=C5_22.at(1).size())
          {
            if(j<=i)
            {
              if(std::abs(C5_22.at(0).at(i)-C5_22.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between SL6 & NU1, difference in trigger times: " << std::abs(C5_22.at(0).at(i)-C5_22.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
          else
          {
            if(i<=j)
            {
              if(std::abs(C5_22.at(0).at(i)-C5_22.at(1).at(j))<coincidenceWindow)
              {
                std::cout << "Coincidence between SL6 & NU1, difference in trigger times: " << std::abs(C5_22.at(0).at(i)-C5_22.at(1).at(j)) << std::endl;
                break_Signal = 1;
                break;
              }
              else
              {
                keepFlag=0;
              }
            }
          }
        }
        if(break_Signal==1)
        {
          break;
        }
      }
    }
    else
    {
      keepFlag=0;
    }

    std::cout << "event=" << event << " Kept=" << keepFlag << std::endl;

    /*
    std::vector<int> HitsEast;
    std::vector<int> HitsWest;

    for (unsigned int ih=0; ih<nchits; ++ih) #
    {
      int auxdetid = (trigs.at(ih))->GetTrigID();
      //      std::cout << ih << "  counter " << auxdetid << std::endl;
      if (auxdetid==9 || auxdetid==8 || auxdetid==7 ||
          auxdetid==93 || auxdetid==94 || auxdetid==95 ) 
      {  
        int tick = trigs.at(ih)->GetTrigTime();
        HitsEast.push_back(tick);
      }
      else if (auxdetid==29 || auxdetid==30 || auxdetid==31 ) 
      {
        int tick = trigs.at(ih)->GetTrigTime();
        HitsWest.push_back(tick);
      }
    }   
    std::cout << " EAST:  "  << HitsEast.size() << 
      " WEST: " << HitsWest.size() << std::endl;

    if (HitsEast.size()>0 && HitsWest.size()>0 )
    {
      for (unsigned int ih=0;ih<HitsEast.size();ih++) 
      {
        std::cout << "east hit " << (int)ih << "  " << HitsEast.at(ih) << std::endl;
      }
      for (unsigned int ih2=0;ih2<HitsWest.size();ih2++) 
      {
        std::cout << "west hit " << (int)ih2 << "  " << HitsWest.at(ih2) << std::endl;
      }
    }

    else {keepFlag=0;}
    */

    return keepFlag;
  }
  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(CWPFilter)

} // namespace filt


