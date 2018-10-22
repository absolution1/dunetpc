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

    unsigned int fLogLevel;
    std::vector<unsigned int> fTimingFlagSelectList;
    std::vector<unsigned int> fTimingFlagDeselectList;
    std::string fTimingLabel;
    std::string fTimingInstance;
    std::string fTriggerLabel;
    std::string fTriggerInstance;
    bool fBeamTrigBool;

  };

  ProtoDUNETriggerFilter::ProtoDUNETriggerFilter::ProtoDUNETriggerFilter(fhicl::ParameterSet const & pset) {
    using std::cout;
    using std::endl;
    const std::string myname = "ProtoDUNETriggerFilter::ctor: ";
    fLogLevel = pset.get<unsigned int>("LogLevel");
    fBeamTrigBool = pset.get<bool>("BeamTrigBool",false);  // just use Leigh's beam selector
    std::vector<unsigned int> defaulttriglist;
    defaulttriglist.push_back(0xc);
    fTimingFlagSelectList = pset.get<std::vector<unsigned int> >("TimingFlagSelectList",defaulttriglist);
    std::vector<unsigned int> emptylist;
    fTimingFlagDeselectList = pset.get<std::vector<unsigned int> >("TimingFlagDeselectList",emptylist);

    fTimingLabel = pset.get<std::string>("TimingLabel","timingrawdecoder");
    fTimingInstance = pset.get<std::string>("TimingInstance","daq");
    fTriggerLabel = pset.get<std::string>("TriggerLabel","ctbrawdecoder");
    fTriggerInstance = pset.get<std::string>("TriggerInstance","daq");
    if ( fLogLevel >= 1 ) {
      cout << myname << "         LogLevel: " << fLogLevel << endl;
      cout << "       TimingFlagSelectList: [";
      bool first = true;
      for ( unsigned int flg : fTimingFlagSelectList ) {
        if ( first ) first = false;
        else cout << ", ";
        cout << flg;
      }
      cout << "]" << endl;
      cout << "     TimingFlagDeselectList: [";
      first = true;
      for ( unsigned int flg : fTimingFlagDeselectList ) {
        if ( first ) first = false;
        else cout << ", ";
        cout << flg;
      }
      cout << "]" << endl;
      cout << myname << "      TimingLabel: " << fTimingLabel << endl;
      cout << myname << "   TimingInstance: " << fTimingLabel << endl;
      cout << myname << "     TriggerLabel: " << fTriggerLabel << endl;
      cout << myname << "  TriggerInstance: " << fTriggerLabel << endl;
    }
  }

  bool ProtoDUNETriggerFilter::filter(art::Event & evt) {
    using std::cout;
    using std::endl;
    const std::string myname = "ProtoDUNETriggerFilter::filter: ";

    if (fBeamTrigBool)
      {
	// The ProtoDUNE data utility tells us if we have a beam trigger
	protoana::ProtoDUNEDataUtils dataUtil;
	return dataUtil.IsBeamTrigger(evt);
      }

    bool result = false;

    // Fetch the trigger and timing clock.

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

      // See https://twiki.cern.ch/twiki/bin/view/CENF/TimingSystemAdvancedOp#Reference_info
      unsigned int trigFlag = tim.GetFlags();

      bool selectflagresult = false;
      if ( fTimingFlagSelectList.size() )
	{
	  for (size_t i=0; i<fTimingFlagSelectList.size(); ++i)
	    {
	      if ( trigFlag == fTimingFlagSelectList.at(i))  // require exact match of the value (not trigger bits but an enum)
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
      if (fTimingFlagDeselectList.size())
	{
	  for (size_t i=0; i<fTimingFlagDeselectList.size(); ++i)
	    {
	      if ( trigFlag == fTimingFlagDeselectList.at(i))  // require exact match of the value (not trigger bits but an enum)
		{
		  deselectflagresult = true;
		  break;
		}
	    }
	}

      result = selectflagresult && (! deselectflagresult);

      // select everything if we have empty input vectors

      if (fTimingFlagSelectList.size() == 0 && fTimingFlagDeselectList.size() == 0)
	{
	  result = true;
	}
    }

    if ( fLogLevel >=2 ) std::cout << myname << "Returning " << (result ? "true" : "false") << endl;
    return result;

  }

  DEFINE_ART_MODULE(ProtoDUNETriggerFilter)

}
