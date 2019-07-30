#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
//#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/Protodune/singlephase/RawDecoding/data/RDStatus.h"

 class PDSPTPCDataInterfaceParent {
  public:
    virtual ~PDSPTPCDataInterfaceParent() noexcept = default;

    // method to get raw digits, RDTimeStamps, RDStatuses and associations from all input fragments specified by an input label (like "daq:ContainerTPC") 

    virtual int retrieveData(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
		   art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, std::vector<raw::RDStatus> &rdstatuses, std::string outputLabel) = 0;

    // method to get raw digits, RDTimeStamps, RDStatuses and associations from all input fragments specified by an input label (like "daq:ContainerTPC") but ony for
    // APA's (== crate numbers) on a list.  If the list contains a -1 in it, it returns all APA data found in the input label.

    virtual int retrieveDataAPAListWithLabels(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
					      art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, std::vector<raw::RDStatus> &rdstatuses, std::string outputLabel,
					      std::vector<int> &apalist) = 0;

    // method to get raw digits, RDTimeStamps, RDStatuses and associations for a specified list of APA's.  The list of possible labels on which to find
    // APA data is proved by fcl configuration.

    virtual int retrieveDataForSpecifiedAPAs(art::Event &evt, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
					     art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, std::vector<raw::RDStatus> &rdstatuses, std::string outputLabel, 
					     std::vector<int> &apalist) = 0;

  };
