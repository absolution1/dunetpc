#ifndef VDColdboxHDF5Utils_h
#define VDColdboxHDF5Utils_h

#include <hdf5.h>
#include <string>
#include <memory>
#include <list>

namespace dune {
namespace VDColdboxHDF5Utils {

struct HDFFileInfo {
  hid_t filePtr;
  size_t bytesWritten;
  std::string fileName;
  int runNumber;
  int dataFormatVersion;
};

typedef std::unique_ptr<HDFFileInfo> HDFFileInfoPtr;
HDFFileInfoPtr openFile(const std::string& fileName);
void closeFile(HDFFileInfoPtr hdfFileInfoPtr);
std::list<std::string> getTopLevelGroupNames(HDFFileInfoPtr& hdfFileInfoPtr);
std::list<std::string> getMidLevelGroupNames(hid_t gid);
bool attrExists(hid_t object, const std::string& attrname);

}
}
#endif
