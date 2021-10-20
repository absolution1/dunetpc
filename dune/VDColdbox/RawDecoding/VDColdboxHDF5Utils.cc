#include "VDColdboxHDF5Utils.h"

namespace dune {
namespace VDColdboxHDF5Utils {

HDFFileInfoPtr openFile(const std::string& fileName) {
  HDFFileInfoPtr hdfFileInfoPtr(new HDFFileInfo());
  hdfFileInfoPtr->filePtr = H5Fopen(fileName.data(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hdfFileInfoPtr->bytesWritten = 0;
  hdfFileInfoPtr->fileName = fileName;
  hdfFileInfoPtr->runNumber = 0;   // don't know it yet.  Fill it when we get fragments.
  hdfFileInfoPtr->dataFormatVersion = 0;

  hid_t grp = H5Gopen(hdfFileInfoPtr->filePtr,"/", H5P_DEFAULT);
  if (attrExists(grp, "data_format_version"))
    {
      hid_t ga = H5Aopen_name(grp, "data_format_version");
      H5Aread(ga, H5Aget_type(ga), &hdfFileInfoPtr->dataFormatVersion);
      H5Aclose(ga);
    }
  H5Gclose(grp);
  return hdfFileInfoPtr;
}

void closeFile(HDFFileInfoPtr hdfFileInfoPtr) {
  H5Fclose(hdfFileInfoPtr->filePtr);
  hdfFileInfoPtr->filePtr = 0;
}

std::list<std::string> getTopLevelGroupNames(HDFFileInfoPtr& hdfFileInfoPtr) {
  hid_t grp = H5Gopen(hdfFileInfoPtr->filePtr,"/", H5P_DEFAULT);
  std::list<std::string> theList = getMidLevelGroupNames(grp);
  H5Gclose(grp);
  return theList;
}

std::list<std::string> getMidLevelGroupNames(hid_t grp) {
  std::list<std::string> theList;
  hsize_t nobj = 0;
  H5Gget_num_objs(grp, &nobj);
  for (hsize_t idx = 0; idx < nobj; ++idx) {
    hsize_t len = H5Gget_objname_by_idx(grp, idx, NULL, 0 );
    char *memb_name = new char(len+1);
    H5Gget_objname_by_idx(grp, idx, memb_name, len+1 );
    theList.emplace_back(memb_name);
    delete[] memb_name;
  }
  return theList;
}

bool attrExists(hid_t object, const std::string &attrname) {
  // Save old error handler 
  H5E_auto_t old_func;
  void *old_client_data;
  H5Eget_auto(H5E_DEFAULT,&old_func, &old_client_data);

  // Turn off error handling */
  H5Eset_auto(H5E_DEFAULT,NULL, NULL);

  // Probe. On failure, retval is supposed to be negative

  hid_t retval = H5Aopen_name(object, attrname.data());

  // Restore previous error handler 
  H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);

  bool result = (retval >= 0);
  return result;
}
    
}
}
