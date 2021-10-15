#include "VDColdboxDataInterface.h"
VDColdboxDataInterface::VDColdboxDataInterface(fhicl::ParameterSet const& p)
{}

int VDColdboxDataInterface::retrieveData(
    art::Event &evt, 
    std::string inputLabel, 
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::vector<raw::RDStatus> &rdstatuses)
{
  return 0;
}

// get data for specified APAs.  Loop over labels specified in the fcl configuration looking for the data so the caller doesn't have to
// keep track of all the branch labels an APA's data might be on.

int VDColdboxDataInterface::retrieveDataForSpecifiedAPAs(
    art::Event &evt, 
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::vector<raw::RDStatus> &rdstatuses, 
    std::vector<int> &apalist)
{
  return 0;
}

int VDColdboxDataInterface::retrieveDataAPAListWithLabels(
    art::Event &evt, 
    std::string inputLabel, 
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::vector<raw::RDStatus> &rdstatuses, 
    std::vector<int> &apalist)
{
  return 0;
}
