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



// Under Development
// TO DO : 1. Convert the event TPC info into an art root file?? Depends on our future need. Still Under Development.
// TO DO : 2. Remove the hardcoded dataset and place a loop to go over hdf5 files.
// The idea is hand this function an hdf5 file and it will return you APA info.
int VDColdboxDataInterface::retrieveTPCData (hid_t fd)
{
  fd = H5Fopen("/dune/data/users/barnali/download",H5F_ACC_RDONLY,H5P_DEFAULT);
  hid_t grp = H5Gopen(fd,"/", H5P_DEFAULT);
  hid_t ga = H5Aopen_name(grp,"data_format_version");
  int dataformatversion=0;
  herr_t ecode = H5Aread(ga,H5Aget_type(ga),&dataformatversion);
  std::cout << "Data Format verison: " << dataformatversion << " Error code: " << ecode << std::endl;
  H5Aclose(ga);
  H5Gclose(grp);
  
  std::list<std::string> theList = dune::VDColdboxHDF5Utils::findTopLevelGroupNames(fd);
  for (auto i : theList)
    {
      std::cout << "Top-Level Group Name: " << i << std::endl;
      std::string topLevelGroupName = i;
      
      // now look inside those "Top-Level Group Name" for "Detector type".
      hid_t requestedGroup = dune::VDColdboxHDF5Utils::getGroupFromPath(fd,i);
      std::list<std::string> detectorTypeNames = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(requestedGroup);
      
      for (auto& detectorTypeName : detectorTypeNames)
	{
	  if (detectorTypeName == "TPC" && detectorTypeName != "TriggerRecordHeader")
	    {
	      std::cout << "  Detector type: " << detectorTypeName << std::endl;
	      std::string subdetGroupPath = i + "/" + detectorTypeName;
	      hid_t subdetGroup = dune::VDColdboxHDF5Utils::getGroupFromPath(fd,subdetGroupPath);
	      std::list<std::string> subdetGeoNames = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(subdetGroup);
	      
	      for (auto& subdetGeoName : subdetGeoNames) // loop over APAs
		{
		  std::string geoGroupPath = subdetGroupPath + "/" + subdetGeoName;
		  std::cout << "     Geo path: " << geoGroupPath << std::endl;
		  hid_t geoGroup = dune::VDColdboxHDF5Utils::getGroupFromPath(fd,geoGroupPath);
		  std::list<std::string> dataSetNames = dune::VDColdboxHDF5Utils::getMidLevelGroupNames(geoGroup);
		  
		  for (auto& dataSetName : dataSetNames) //// loop over Links
		    {
		      std::string dataSetPath = geoGroupPath + "/" + dataSetName;
		      std::cout << "      Data Set Path: " << dataSetPath << std::endl;
		      hid_t datasetid = H5Dopen(geoGroup,dataSetName.data(),H5P_DEFAULT);
		      hsize_t ds_size = H5Dget_storage_size(datasetid);
		      std::cout << "      Data Set Size (bytes): " << ds_size << std::endl;
		      
		      if (ds_size == 0) continue;
		      if (ds_size < 80)
			{
			  std::cout << "TPC datset too small for the fragment header" << std::endl;
			}
		      
		      size_t narray = ds_size / sizeof(char);
		      size_t rdr = ds_size % sizeof(char);
		      if (rdr > 0 || narray == 0) narray++;
		      char *ds_data = new char[narray];
		      ecode = H5Dread(datasetid, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ds_data);
		      int firstbyte = ds_data[0];
		      firstbyte &= 0xFF;
		      int lastbyte = ds_data[narray-1];
		      lastbyte &= 0xFF;
		      
		      std::cout << std::hex << "      Retrieved data: ecode: " << ecode << "  first byte: " << firstbyte
				<< " last byte: " << lastbyte  << std::dec << std::endl;
		      
		      H5Dclose(datasetid); 
		      
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
		  
		  H5Gclose(geoGroup);
                }
              
	      H5Gclose(subdetGroup);
            }
	}
      H5Gclose(requestedGroup);
    }
  H5Fclose(fd);

  return 0;  
}



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
