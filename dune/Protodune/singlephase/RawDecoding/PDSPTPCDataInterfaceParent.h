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
    virtual int retrieveData(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
		   art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, std::vector<raw::RDStatus> &rdstatuses, std::string outputLabel) = 0;

  };
