#ifndef VDColdboxDataInterface_H
#define VDColdboxDataInterface_H
#include "dune/DuneObj/PDSPTPCDataInterfaceParent.h"
#include "fhiclcpp/ParameterSet.h"

class VDColdboxDataInterface : public PDSPTPCDataInterfaceParent {

 public:

  VDColdboxDataInterface(fhicl::ParameterSet const& ps);
  int retrieveData(art::Event &evt, std::string inputlabel,
                   std::vector<raw::RawDigit> &raw_digits,
                   std::vector<raw::RDTimeStamp> &rd_timestamps,
                   std::vector<raw::RDStatus> &rdstatuses );

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



};

#endif
