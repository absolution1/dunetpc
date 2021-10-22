#include "VDColdboxDataInterface.h"

#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <set>

#include "art/Framework/Services/Registry/ServiceHandle.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/RceFragment.hh"
#include "dune-raw-data/Overlays/FelixFragment.hh"
#include "dune-raw-data/Overlays/RceFragment.hh"
#include "dune-raw-data/Overlays/FelixFragment.hh"
#include "dune-raw-data/Overlays/FragmentType.hh"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "dune/DuneObj/DUNEHDF5FileInfo.h"
#include "dune/VDColdbox/RawDecoding/VDColdboxHDF5Utils.h"



// Keep this here for now. This needs scrutiny regarding whether the input labels are fetching right values.

VDColdboxDataInterface::VDColdboxDataInterface(fhicl::ParameterSet const& p) {
  //_input_labels_by_apa[1] = p.get< std::vector<std::string> >("APA1InputLabels");
  //  _input_labels_by_apa[2] = p.get< std::vector<std::string> >("APA2InputLabels");
  //  _input_labels_by_apa[3] = p.get< std::vector<std::string> >("APA3InputLabels");
  //  _input_labels_by_apa[4] = p.get< std::vector<std::string> >("APA4InputLabels");
  //  _input_labels_by_apa[5] = p.get< std::vector<std::string> >("APA5InputLabels");
  // _input_labels_by_apa[6] = p.get< std::vector<std::string> >("APA6InputLabels");
  // _input_labels_by_apa[7] = p.get< std::vector<std::string> >("APA7InputLabels");
  // _input_labels_by_apa[8] = p.get< std::vector<std::string> >("APA8InputLabels");
  // _input_labels_by_apa[-1] = p.get< std::vector<std::string> >("MISCAPAInputLabels");
}


int VDColdboxDataInterface::retrieveData(
    art::Event &evt, 
    std::string inputLabel, 
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::vector<raw::RDStatus> &rdstatuses) {
  return 0;
}


// get data for specified APAs.  Loop over labels specified in the fcl configuration looking for the data so 
// the caller doesn't have to keep track of all the branch labels an APA's data might be on.

//Fhicl Config. must have to be taken care of in future

int VDColdboxDataInterface::retrieveDataForSpecifiedAPAs(
    art::Event &evt, 
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::vector<raw::RDStatus> &rdstatuses, 
    std::vector<int> &apalist) {

  auto infoHandle = evt.getHandle<raw::DUNEHDF5FileInfo>("daq");
  std::cout << "Got infos? " << infoHandle << std::endl;
  //Add check for infoHandle?

  //NOTE: this bit of code is just for testing/demonstration
  const std::string & event_group = infoHandle->GetEventGroupName();
  const std::string & file_name = infoHandle->GetFileName();
  //const hid_t & hdf_file = infoHandle->GetHDF5FileHandle();
  std::cout << "\t" << file_name << std::endl;
  std::cout << "\t" << infoHandle->GetFormatVersion() << std::endl;
  std::cout << "\t" << event_group << std::endl;

  hid_t hdf_file = H5Fopen(file_name.data(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t the_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
      hdf_file, event_group);


  std::list<std::string> det_types
      = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(the_group);
  std::cout << "\tDet types: " << det_types.size() << std::endl;
  /////////////////////////////////////////////


  //Here: need to form Fragments -- reimplement getFragmentsForEvent from HDFFileReader
  //Then: process these as in PDSPTPCDataInterface

  int totretcode = 0;
 
  /*
  for (size_t i=0; i<apalist.size(); ++i) 
  { 
    auto lli = _input_labels_by_apa.find(apalist.at(i)); 
    if (lli == _input_labels_by_apa.end())
      {
	MF_LOG_WARNING("PDSPTPCDataInterface:") << " No list of input labels known for APA " << apalist.at(i) << " Returning no data.";
      }
    for (size_t j=0; j<lli->second.size(); ++j)
      {
	int retcode = retrieveDataAPAListWithLabels(evt, lli->second.at(j), raw_digits, rd_timestamps, rdstatuses, apalist );
	if (retcode > totretcode) totretcode = retcode; // take most severe retcode of everything
      }
  }
  _collectRDStatus(rdstatuses);*/
 
  return totretcode;
}


// get data for a specific label, but only return those raw digits that correspond to APA's on the list
// 
int VDColdboxDataInterface::retrieveDataAPAListWithLabels(
    art::Event &evt, 
    std::string inputLabel, 
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::vector<raw::RDStatus> &rdstatuses, 
    std::vector<int> &apalist) {
  return 0;
}

DEFINE_ART_CLASS_TOOL(VDColdboxDataInterface)
