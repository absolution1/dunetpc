#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/RawDigit.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

protoana::ProtoDUNEDataUtils::ProtoDUNEDataUtils(fhicl::ParameterSet const& p){
  this->reconfigure(p);
}

protoana::ProtoDUNEDataUtils::~ProtoDUNEDataUtils(){

}

void protoana::ProtoDUNEDataUtils::reconfigure(fhicl::ParameterSet const& p){
  fTimingTag            = p.get<art::InputTag>("TimingTag");
  fRawDigitTag          = p.get<art::InputTag>("RawDigitTag");
  fRawDigitTimeStampTag = p.get<art::InputTag>("RawDigitTimeStampTag");
}

// Access the trigger information to see if this is a beam trigger
bool protoana::ProtoDUNEDataUtils::IsBeamTrigger(art::Event const & evt) const{

  bool isBeam = false;

  // Accessing the trigger information as done in DataPrepModule
  // The information is stored in the time stamps
  art::Handle<std::vector<raw::RDTimeStamp>> timeStamps;
  evt.getByLabel(fTimingTag,timeStamps);

  // Return false if we have no time stamps
  if(!timeStamps.isValid()) return isBeam;
  // We should only have one RDTimeStamp
  if(timeStamps->size() > 1) return isBeam;
  
  // Access the trigger information. Beam trigger flag = 0xc
  const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
  isBeam = (timeStamp.GetFlags() == 0xc);
  
  return isBeam;
}

// ----------------------------------------------------------------------------
int protoana::ProtoDUNEDataUtils::GetNActiveFembsForAPA(art::Event const & evt, int apa) const {


// Get pd channel map
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;

  // set only saves unique elements
  std::set<int> apaset;

  // Get raw digits time stamps
  art::Handle< std::vector<raw::RDTimeStamp> > RawdigitTSListHandle;
  std::vector<art::Ptr<raw::RDTimeStamp> > digitTSlist;
  

   // Get raw digits
  art::Handle< std::vector<raw::RawDigit> > RawdigitListHandle;
  std::vector<art::Ptr<raw::RawDigit> > digitlist;

  if (evt.getByLabel(fRawDigitTag, RawdigitListHandle)){
  	
    art::fill_ptr_vector(digitlist, RawdigitListHandle);  

    for(auto const & dptr : digitlist) {
    const raw::RawDigit& digit = *dptr;
    // Get the channel number for this digit
    uint32_t chan = digit.Channel();
    
    int iapa = channelMap->APAFromOfflineChannel(chan);
    if(iapa != apa) continue;
    // Get the channel FEMB and WIB
    int WIB = channelMap->WIBFromOfflineChannel(chan); // 0-4
    int FEMB = channelMap->FEMBFromOfflineChannel(chan); // 1-4
    //int FEMBchan = channelMap->FEMBChannelFromOfflineChannel(chan);
    int iFEMB = ((WIB*4)+(FEMB-1)); //index of the FEMB 0-19

    apaset.insert(iFEMB);
  }
}

else{ // if raw digits have been dropped use RDTimeStamps instead
	evt.getByLabel(fRawDigitTimeStampTag, RawdigitTSListHandle);
	
	art::fill_ptr_vector(digitTSlist, RawdigitTSListHandle);  
	
  for(auto const & dptr : digitTSlist) {
  	
    const raw::RDTimeStamp & digit = *dptr;
    
    // Get the channel number for this digit
    uint16_t chan = digit.GetFlags();
    
    int iapa = channelMap->APAFromOfflineChannel(chan);
    if(iapa != apa) continue;
    // Get the channel FEMB and WIB
    int WIB = channelMap->WIBFromOfflineChannel(chan); // 0-4
    int FEMB = channelMap->FEMBFromOfflineChannel(chan); // 1-4
    //int FEMBchan = channelMap->FEMBChannelFromOfflineChannel(chan);
    int iFEMB = ((WIB*4)+(FEMB-1)); //index of the FEMB 0-19

    apaset.insert(iFEMB);
  }
}

  return (apaset.size());

}




