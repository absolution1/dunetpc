#ifndef VDColdboxDataInterface_H
#define VDColdboxDataInterface_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "artdaq-core/Data/Fragment.hh"
#include "dune/DuneObj/PDSPTPCDataInterfaceParent.h"
#include "daqdataformats/Fragment.hpp"
#include <hdf5.h>

typedef dunedaq::daqdataformats::Fragment duneFragment;
typedef std::vector<duneFragment> duneFragments; 

class VDColdboxDataInterface : public PDSPTPCDataInterfaceParent {

 public:

  VDColdboxDataInterface(fhicl::ParameterSet const& ps);
  int retrieveData(art::Event &evt, std::string inputlabel,
                   std::vector<raw::RawDigit> &raw_digits,
                   std::vector<raw::RDTimeStamp> &rd_timestamps,
                   std::vector<raw::RDStatus> &rdstatuses );

  void readFragmentsForEvent (art::Event &evt);

  int retrieveDataAPAListWithLabels(
      art::Event &evt, std::string inputlabel,
      std::vector<raw::RawDigit> &raw_digits,
      std::vector<raw::RDTimeStamp> &rd_timestamps,
      std::vector<raw::RDStatus> &rdstatuses, 
      std::vector<int> &apalist);

  int retrieveDataForSpecifiedAPAs(
      art::Event &evt, std::vector<raw::RawDigit> &raw_digits,
      std::vector<raw::RDTimeStamp> &rd_timestamps,
      std::vector<raw::RDStatus> &rdstatuses,  
      std::vector<int> &apalist);



 private:

  std::map<int,std::vector<std::string>> _input_labels_by_apa;
  void _collectRDStatus(std::vector<raw::RDStatus> &rdstatuses){};

  //For nicer log syntax
  std::string logname = "VDColdboxDataInterface";
  bool fForceOpen;
  hid_t fPrevStoredHandle = -1;
  hid_t fHDFFile = -1;
};

#endif
