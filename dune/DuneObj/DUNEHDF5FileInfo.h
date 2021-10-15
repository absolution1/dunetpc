#ifndef DUNEHDF5FileInfo_h
#define DUNEHDF5FileInfo_h
#include <hdf5.h>
namespace raw {
class DUNEHDF5FileInfo {
 private:
  std::string fFileName;
  hid_t fHDF5FileHandle;
  int fFormatVersion;
  std::string fEventGroupName;
};
}
#endif
