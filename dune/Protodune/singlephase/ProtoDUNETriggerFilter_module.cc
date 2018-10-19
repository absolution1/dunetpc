// A very simple module to filter out the events with beam triggers
// leigh.howard.whitehead@cern.ch

#include <iostream>
#include <utility>
#include <set>

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 

#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"

namespace filt{

  class ProtoDUNETriggerFilter : public art::EDFilter {
  public:
    explicit ProtoDUNETriggerFilter(fhicl::ParameterSet const & pset);
    virtual ~ProtoDUNETriggerFilter() {};
    virtual bool filter(art::Event& e);

  private:

    std::vector<unsigned int> fTimingFlagSel;
    std::vector<unsigned int> fTimingFlagDeSel;
    std::string fTimingLabel;
    std::string fTimingInstance;
    std::string fTriggerLabel;
    std::string fTriggerInstance;
    bool fBeamTrigBool;

  };

  ProtoDUNETriggerFilter::ProtoDUNETriggerFilter::ProtoDUNETriggerFilter(fhicl::ParameterSet const & pset)
  {
    fBeamTrigBool = pset.get<bool>("BeamTrigBool",false);  // just use Leigh's beam selector
    std::vector<unsigned int> defaulttriglist;
    defaulttriglist.push_back(0xc);
    fTimingFlagSel = pset.get<std::vector<unsigned int> >("TimingFlagSelectList",defaulttriglist);
    std::vector<unsigned int> emptylist;
    fTimingFlagDeSel = pset.get<std::vector<unsigned int> >("TimingFlagDeselectList",emptylist);

    fTimingLabel = pset.get<std::string>("TimingLabel","timingrawdecoder");
    fTimingInstance = pset.get<std::string>("TimingInstance","daq");
    fTriggerLabel = pset.get<std::string>("TriggerLabel","ctbrawdecoder");
    fTriggerInstance = pset.get<std::string>("TriggerInstance","daq");
  }

  bool ProtoDUNETriggerFilter::filter(art::Event & evt){

    if (fBeamTrigBool)
      {
	// The ProtoDUNE data utility tells us if we have a beam trigger
	protoana::ProtoDUNEDataUtils dataUtil;
	return dataUtil.IsBeamTrigger(evt);
      }

    bool result = false;

    // Fetch the trigger and timing clock.

    std::string myname = "PDSPTriggerFilter_module.cc: ";

    art::Handle<std::vector<raw::RDTimeStamp>> htims;
    evt.getByLabel(fTimingLabel, fTimingInstance, htims);

    art::Handle<std::vector<raw::ctb::pdspctb> > hctb;
    evt.getByLabel(fTriggerLabel, fTriggerInstance, hctb);

    // using the ctb triggers is not yet implemented

    if ( ! htims.isValid() ) {
      std::cout << myname << "WARNING: Timing clocks product not found." << std::endl;
    } else if (  htims->size() != 1 ) {
      std::cout << myname << "WARNING: Unexpected timing clocks size: " << htims->size() << std::endl;
      for ( unsigned int itim=0; itim<htims->size() && itim<50; ++itim ) {
	std::cout << myname << "  " << htims->at(itim).GetTimeStamp() << std::endl;
      }
    } else {
      const raw::RDTimeStamp& tim = htims->at(0);

      //timingClock = tim.GetTimeStamp();
      // See https://twiki.cern.ch/twiki/bin/view/CENF/TimingSystemAdvancedOp#Reference_info
      //trigFlag = tim.GetFlags();

      bool selectflagresult = false;
      if (fTimingFlagSel.size())
	{
	  for (size_t i=0; i<fTimingFlagSel.size(); ++i)
	    {
	      if ( (tim.GetFlags() & fTimingFlagSel.at(i)) == fTimingFlagSel.at(i))  // require exact match of all bits in selection list 
		{
		  selectflagresult = true;
		  break;
		}
	    }
	}
      else
	{
	  selectflagresult = true;
	}

      bool deselectflagresult = false;
      if (fTimingFlagDeSel.size())
	{
	  for (size_t i=0; i<fTimingFlagDeSel.size(); ++i)
	    {
	      if ( (tim.GetFlags() & fTimingFlagDeSel.at(i)) ==  fTimingFlagDeSel.at(i))  // require exact match of all bits in selection list 
		{
		  deselectflagresult = true;
		  break;
		}
	    }
	}

      result = selectflagresult && (! deselectflagresult);

      // select everything if we have empty input vectors

      if (fTimingFlagSel.size() == 0 && fTimingFlagDeSel.size() == 0)
	{
	  result = true;
	}

      // tim.GetFlags is 0xc if beam.  
      //std::cout << myname << "Trigger flag: " << trigFlag << " (";
      //bool isBeam = trigFlag == 0xc;
      //bool isFake = trigFlag >= 0x8 && trigFlag <= 0xb;
      //if ( isBeam ) std::cout << "Beam";
      //else if ( isFake ) std::cout << "Fake";
      //else std::cout << "Unexpected";
      //std::cout << ")" << std::endl;
    }

    return result;


  }

  DEFINE_ART_MODULE(ProtoDUNETriggerFilter)

}
