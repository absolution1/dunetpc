#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"

#include "lardataobj/RawData/RDTimeStamp.h"

protoana::ProtoDUNEDataUtils::ProtoDUNEDataUtils(){

}

protoana::ProtoDUNEDataUtils::~ProtoDUNEDataUtils(){

}

// Access the trigger information to see if this is a beam trigger
bool protoana::ProtoDUNEDataUtils::IsBeamTrigger(art::Event const & evt) const{

  bool isBeam = false;

  // Accessing the trigger information as done in DataPrepModule
  // The information is stored in the time stamps
  art::Handle<std::vector<raw::RDTimeStamp>> timeStamps;
  evt.getByLabel("timingrawdecoder","daq",timeStamps);

  // Return false if we have no time stamps
  if(!timeStamps.isValid()) return isBeam;
  // We should only have one RDTimeStamp
  if(timeStamps->size() > 1) return isBeam;
  
  // Access the trigger information. Beam trigger flag = 0xc
  const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
  isBeam = (timeStamp.GetFlags() == 0xc);
  
  return isBeam;
}

