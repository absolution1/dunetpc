#include "VDColdboxHDF5Utils.h"
#include "VDColdboxDataInterface.h"

#include <hdf5.h>
#include <iostream>
#include <list>
#include <set>
#include <sstream>
#include <cstring>
#include <string>
#include "TMath.h"
#include "TString.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

// artdaq and dune-raw-data includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "dune/DuneObj/DUNEHDF5FileInfo.h"
#include "dune/VDColdbox/RawDecoding/VDColdboxHDF5Utils.h"


//The file handle is from the raw::DUNEHDF5FileInfo data product that the source puts into the event. Art's getHandle<type> is usedto retrieve a data product from the event.  
// The idea is hand this function an art event and it will return you APA/CRU info FOR THE VDColdbox.

void readFragmentsForEvent (art::Event &evt)
{
  auto infoHandle = evt.getHandle <raw::DUNEHDF5FileInfo> ("daq");
  std::cout << "Got infos? " << infoHandle << std::endl;
  
  const std::string & toplevel_groupname = infoHandle->GetEventGroupName();
  const std::string & file_name = infoHandle->GetFileName();
  hid_t file_id = infoHandle->GetHDF5FileHandle();

  std::cout << "Top-Level Group Name: " << toplevel_groupname << std::endl;
  std::cout << "HDF5 FileName: " << file_name << std::endl;
  
  // now look inside those "Top-Level Group Name" for "Detector type".
  hid_t requestedGroup = dune::VDColdboxHDF5Utils::getGroupFromPath(file_id, toplevel_groupname);

  std::list<std::string> detectorTypeNames = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(requestedGroup);
  
  for (auto& detectorTypeName : detectorTypeNames)
    {
      if (detectorTypeName == "TPC" && detectorTypeName != "TriggerRecordHeader")
	{
	  std::cout << "  Detector type: " << detectorTypeName << std::endl;
	  std::string geoPath = toplevel_groupname + "/" + detectorTypeName;
	  hid_t geoGroup = dune::VDColdboxHDF5Utils::getGroupFromPath(file_id,geoPath);
	  std::list<std::string> apaNames = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(geoGroup);
	  
	  // loop over APAs
	  for (auto& apaName : apaNames)
	    {
	      std::string apaGroupPath = geoPath + "/" + apaName;
	      std::cout << "     Geo path: " << apaGroupPath << std::endl;
	      hid_t linkGroup = dune::VDColdboxHDF5Utils::getGroupFromPath(file_id,apaGroupPath);
	      std::list<std::string> linkNames = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(linkGroup);
	      
	      // loop over Links
	      for (auto& linkName : linkNames)
		{
		  std::string dataSetPath = apaGroupPath + "/" + linkName;
		  std::cout << "      Data Set Path: " << dataSetPath << std::endl;
		  hid_t datasetid = H5Dopen(linkGroup,linkName.data(),H5P_DEFAULT);
		  hsize_t ds_size = H5Dget_storage_size(datasetid);
		  std::cout << "      Data Set Size (bytes): " << ds_size << std::endl;
		  
		  if (ds_size < 80) continue;
		  
		  size_t narray = ds_size / sizeof(char);
		  size_t rdr = ds_size % sizeof(char);
		  if (rdr > 0 || narray == 0) narray++;
		  char *ds_data = new char[narray];
		  herr_t ecode = H5Dread(datasetid, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ds_data);
		  int firstbyte = ds_data[0];
		  firstbyte &= 0xFF;
		  int lastbyte = ds_data[narray-1];
		  lastbyte &= 0xFF;
		  
		  std::cout << std::hex << "      Retrieved data: ecode: " << ecode << "  first byte: " << firstbyte
			    << " last byte: " << lastbyte  << std::dec << std::endl;
		  

		  int magic_word = 0;
		  memcpy(&magic_word,&ds_data[0],4);
		  std::cout << "   Magic word: 0x" << std::hex << magic_word << std::dec << std::endl;
		  
		  int version = 0;
		  memcpy(&version, &ds_data[4],4);
		  std::cout << "   Version: " << std::dec << version << std::dec << std::endl;
		  
		  uint64_t fragsize=0;
		  memcpy(&fragsize, &ds_data[8],8);
		  std::cout << "   Frag Size: " << std::dec << fragsize << std::dec << std::endl;
		  
		  uint64_t trignum=0;
		  memcpy(&trignum, &ds_data[16],8);
		  std::cout << "   Trig Num: " << std::dec << trignum << std::dec << std::endl;
		  
		  uint64_t trig_timestamp=0;
		  memcpy(&trig_timestamp, &ds_data[24],8);
		  std::cout << "   Trig Timestamp: " << std::dec << trig_timestamp << std::dec << std::endl;
		  
		  uint64_t windowbeg=0;
		  memcpy(&windowbeg, &ds_data[32],8);
		  std::cout << "   Window Begin:   " << std::dec << windowbeg << std::dec << std::endl;
		  
		  uint64_t windowend=0;
		  memcpy(&windowend, &ds_data[40],8);
		  std::cout << "   Window End:     " << std::dec << windowend << std::dec << std::endl;
		  
		  int runno=0;
		  memcpy(&runno, &ds_data[48], 4);
		  std::cout << "   Run Number: " << std::dec << runno << std::endl;
		  
		  int errbits=0;
		  memcpy(&errbits, &ds_data[52], 4);
		  std::cout << "   Error bits: " << std::dec << errbits << std::endl;
		  
		  int fragtype=0;
		  memcpy(&fragtype, &ds_data[56], 4);
		  std::cout << "   Fragment type: " << std::dec << fragtype << std::endl;
		  
		  int fragpadding=0;
		  memcpy(&fragtype, &ds_data[60], 4);
		  std::cout << "   Fragment padding: " << std::dec << fragpadding << std::endl;
		  
		  int geoidversion=0;
		  memcpy(&geoidversion, &ds_data[64], 4);
		  std::cout << "   GeoID version: " << std::dec << geoidversion << std::endl;
		  
		  unsigned short geoidtype;
		  memcpy(&geoidtype, &ds_data[70], 1);
		  std::cout << "   GeoID type: " << geoidtype << std::endl;
		  
		  unsigned short geoidregion=0;
		  memcpy(&geoidregion, &ds_data[71], 1);
		  std::cout << "   GeoID region: " << std::dec << geoidregion << std::endl;
		  
		  int geoidelement=0;
		  memcpy(&geoidelement, &ds_data[72], 4);
		  std::cout << "   GeoID element: " << std::dec << geoidelement << std::endl;
		  
		  int geoidpadding=0;
		  memcpy(&geoidpadding, &ds_data[76], 4);
		  std::cout << "   GeoID padding: " << std::dec << geoidpadding << std::endl;
		  
		  delete[] ds_data;  // free up memory

		} 
	    }
	}
    }

}

// Keep this here for now. This needs scrutiny regarding whether the input labels are fetching right values.

VDColdboxDataInterface::VDColdboxDataInterface(fhicl::ParameterSet const& p)
  : fForceOpen(p.get<bool>("ForceOpen", false)) {
   
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


// TO DO :May be we should CONSIDER using "retrieveTPCData" function from above TO DEFINE "retrieveDataForSpecifiedAPAs" once the hdf5 info are placed in an art root file?? UNDER DISCUSSION.

// get data for specified APAs.  Loop over labels specified in the fcl configuration looking for the data so 
// the caller doesn't have to keep track of all the branch labels an APA's data might be on.

//Fhicl Config. must have to be taken care of in future

int VDColdboxDataInterface::retrieveDataForSpecifiedAPAs(
    art::Event &evt, 
    std::vector<raw::RawDigit> &raw_digits, 
    std::vector<raw::RDTimeStamp> &rd_timestamps,
    std::vector<raw::RDStatus> &rdstatuses, 
    std::vector<int> &apalist) {

  //if (apalist.size() > 0) {
  //  if (apalist[0] > 1) return 0;
  //}

  std::cout << "Retrieving Data for " << apalist.size() << " APAs: ";
  for (const int & i : apalist)
    std::cout << i << " ";
  std::cout << std::endl;

  //Turn "daq" --> fcl parameter defined within constructor
  auto infoHandle = evt.getHandle<raw::DUNEHDF5FileInfo>("daq");
  //Add check for infoHandle?

  //NOTE: this bit of code is just for testing/demonstration
  const std::string & event_group = infoHandle->GetEventGroupName();
  const std::string & file_name = infoHandle->GetFileName();
  std::cout << "\t" << file_name << std::endl;
  std::cout << "\t" << infoHandle->GetFormatVersion() << std::endl;
  std::cout << "\t" << event_group << std::endl;

  //If the fcl file said to force open the file
  //(i.e. because one is just running DataPrep), then open
  //but only if we are on a new file -- identified by if the handle
  //stored in the event is different
  hid_t stored_handle = infoHandle->GetHDF5FileHandle();
  if (fForceOpen && (stored_handle != fPrevStoredHandle)) {
    std::cout << "Opening" << std::endl;
    fHDFFile = H5Fopen(file_name.data(), H5F_ACC_RDONLY, H5P_DEFAULT);
  }//If the handle is the same, fHDFFile won't change 
  else if (!fForceOpen) {
    fHDFFile = stored_handle;
  }
  fPrevStoredHandle = stored_handle;

  hid_t the_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
      fHDFFile, event_group);


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

  dune::VDColdboxHDF5Utils::getFragmentsForEvent(fHDFFile, event_group,
                                                 raw_digits, rd_timestamps);
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
*/
 

  //Currently putting in dummy values for the RD Statuses
  rdstatuses.clear();
  rdstatuses.emplace_back(false, false, 0);

  return totretcode;
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

DEFINE_ART_CLASS_TOOL(VDColdboxDataInterface)
