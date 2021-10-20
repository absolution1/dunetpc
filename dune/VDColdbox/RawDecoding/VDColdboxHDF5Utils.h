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

struct HeaderInfo {
  int magicWord = 0;
  int version = 0;
  uint64_t trigNum = 0;
  uint64_t trigTimestamp = 0;
  uint64_t nReq = 0;
  int runNum = 0;
  int errBits = 0;
  short triggerType = 0;
};

typedef std::unique_ptr<HDFFileInfo> HDFFileInfoPtr;
HDFFileInfoPtr openFile(const std::string& fileName);
void closeFile(HDFFileInfoPtr hdfFileInfoPtr);
std::list<std::string> getTopLevelGroupNames(HDFFileInfoPtr& hdfFileInfoPtr);
std::list<std::string> getMidLevelGroupNames(hid_t gid);
bool attrExists(hid_t object, const std::string& attrname);
hid_t getGroupFromPath(hid_t fd, const std::string &path);

void getHeaderInfo(hid_t the_group, const std::string & det_type,
                   HeaderInfo & info);


}
}
#endif
