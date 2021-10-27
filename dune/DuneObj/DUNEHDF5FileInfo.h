#ifndef DUNEHDF5FileInfo_h
#define DUNEHDF5FileInfo_h
#include <hdf5.h>
namespace raw {
class DUNEHDF5FileInfo {
 public:
  DUNEHDF5FileInfo(std::string filename, hid_t file_handle, int version,
                   std::string group)
   : fFileName(filename), fHDF5FileHandle(file_handle),
     fFormatVersion(version), fEventGroupName(group) {};
  DUNEHDF5FileInfo(std::initializer_list<DUNEHDF5FileInfo>) {};

  const std::string & GetFileName() const {return fFileName;};
  const int & GetFormatVersion() const {return fFormatVersion;};
  const std::string & GetEventGroupName() const {return fEventGroupName;};
  const hid_t & GetHDF5FileHandle() const {return fHDF5FileHandle;};

 private:
  std::string fFileName;
  hid_t fHDF5FileHandle;
  int fFormatVersion;
  std::string fEventGroupName;
};
}
#endif
