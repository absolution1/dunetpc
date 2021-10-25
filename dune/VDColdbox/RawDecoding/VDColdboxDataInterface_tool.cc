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

  //Turn "daq" --> fcl parameter defined within constructor
  auto infoHandle = evt.getHandle<raw::DUNEHDF5FileInfo>("daq");
  //Add check for infoHandle?

  //NOTE: this bit of code is just for testing/demonstration
  const std::string & event_group = infoHandle->GetEventGroupName();
  const std::string & file_name = infoHandle->GetFileName();
  std::cout << "\t" << file_name << std::endl;
  std::cout << "\t" << infoHandle->GetFormatVersion() << std::endl;
  std::cout << "\t" << event_group << std::endl;

  hid_t hdf_file = H5Fopen(file_name.data(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t the_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
      hdf_file, event_group);


  std::list<std::string> det_types
      = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(the_group);
  std::cout << "\tDet types: " << det_types.size() << std::endl << std::endl;
  std::string tpc_path = event_group + "/TPC";
  std::cout << "Attempting to open " << tpc_path << std::endl;
  hid_t tpc_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
      the_group, "TPC");
  std::list<std::string> subdet_types
      = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(tpc_group);
  std::cout << "\tSubdet types: " << subdet_types.size() << std::endl <<
               std::endl;
  for (const auto & t : subdet_types) {
    std::cout << "\t" << t << std::endl;
  }

  hid_t apa_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
      tpc_group, "APA000");
  std::list<std::string> link_names
      = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(apa_group);
  std::cout << "\tLink types: " << link_names.size() << std::endl <<
               std::endl;
  for (const auto & t : link_names) {
    std::cout << "\t" << t << std::endl;
  }
  /////////////////////////////////////////////

  /*auto frags
      = */dune::VDColdboxHDF5Utils::getFragmentsForEvent(hdf_file,
                                                         event_group,
                                                         raw_digits,
                                                         rd_timestamps);
  //std::map<std::string, std::unique_ptr<duneFragments>> frags;
  //std::map<std::string, duneFragments> frags;
  //getFragmentsForEvent(hdf_file, event_group, frags);

  //for (auto it = frags.begin(); it != frags.end(); ++it) {
    //std::cout << "TPC" << " has " << frags["TPC"]->size() << std::endl;
    //for (size_t i = 0; i < 8; ++i) {
    //  std::cout << "\t" << (frags["TPC"])[i].get_size() << std::endl;
    //}
    //processFragments(raw_digits, rd_timestamps, frags["TPC"]);
//  }

  //Here: need to form Fragments -- reimplement getFragmentsForEvent from HDFFileReader
  //$DUNE_RAW_DATA_DIR/source/dune-raw-data/HDFUtils/HDFFileReader.cc
  //There's a block for TPC data starting at 172, that should work for our needs
  //
  //
  //
  //Then: process these as in PDSPTPCDataInterface
  //dunetpc/dune/Protodune/singlephase/RawDecoding/PDSPTPCDataInterface_tool.cc
  //within retrieveDataForSpecifiedAPAs(...):
  //   retrieveDataAPAListWithLabels(...):
  //       _processFELIX(...):
  //           _felixProcContNCFrags(...):
  //               _process_FELIX_AUX(...) <-- This is where the fragments are actually converted 
  //                                           into RawDigits

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
 

  //Currently putting in dummy values for the RD Statuses
  rdstatuses.clear();
  rdstatuses.emplace_back(false, false, 0);

  return totretcode;
}

void VDColdboxDataInterface::processFragments(
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::unique_ptr<duneFragments> & frags) {
  for (const duneFragment & frag : *frags) {
    std::cout << frag.get_size() << std::endl;  
  }
}

// get data for a specific label, but only return those raw digits that correspond to APA's on the list
int VDColdboxDataInterface::retrieveDataAPAListWithLabels(
    art::Event &evt, 
    std::string inputLabel, 
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::vector<raw::RDStatus> &rdstatuses, 
    std::vector<int> &apalist) {
  return 0;
}

void VDColdboxDataInterface::getFragmentsForEvent(
    hid_t hdf_file, const std::string & group_name,
    std::map<std::string, duneFragments> & results) {
  std::cout << "Frag" << std::endl;
  hid_t the_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
      hdf_file, group_name);

  std::list<std::string> det_types
      = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(the_group);
  for (const auto & det : det_types) {
    if (det != "TPC") continue;
    results["TPC"] = duneFragments();
    std::cout << "\t" << det << std::endl;
    hid_t det_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
        the_group, det);
    std::list<std::string> subdet_types
        = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(det_group);
    for (const auto & subdet: subdet_types) {
      std::cout << "\t\t" << subdet << std::endl;
      hid_t subdet_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
          det_group, subdet);
      std::list<std::string> link_names
          = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(subdet_group);
      for (const auto & t : link_names) {
        std::cout << "\t\t\t" << t << std::endl;
        
        hid_t dataset = H5Dopen(subdet_group, t.data(), H5P_DEFAULT);
        hsize_t ds_size = H5Dget_storage_size(dataset);
        if (ds_size <= sizeof(dunedaq::daqdataformats::FragmentHeader)) continue; //Too small

        std::vector<char> ds_data(ds_size);   
        H5Dread(dataset, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                ds_data.data());
        H5Dclose(dataset);

        std::unique_ptr<duneFragment> frag
            = std::make_unique<duneFragment>(
                &ds_data[0], duneFragment::BufferAdoptionMode::kReadOnlyMode);
        std::cout << "Got fragment of size: " << frag->get_size() << std::endl;
        results["TPC"].emplace_back(std::move(*frag.release()));
        std::cout << results["TPC"].back().get_size() << std::endl;

        //Make 256 --> n_channels configurable?
        //size_t n_frames = (ds_size - sizeof(FragmentHeader))/256;
        //std::cout << "N frames: " << n_frames << std::endl;
        //for (size_t iFrame = 0; iFrame < n_frames; ++iFrame) {
        //  size_t start = iFrame*256;
        //  WIBFrame this_frame;
        //  for (size_t iChan = 0; iChan < 256; ++iChan) {
        //    this_frame.set_channel(iChan, ds_data[80 + start + iChan]);
        //  }
        //}
      }
      H5Gclose(subdet_group);
    }
    H5Gclose(det_group);
  }
  H5Gclose(the_group);
}

DEFINE_ART_CLASS_TOOL(VDColdboxDataInterface)
