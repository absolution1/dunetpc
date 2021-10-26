#ifndef VDColdboxHDF5Utils_h
#define VDColdboxHDF5Utils_h

#include "artdaq-core/Data/Fragment.hh"

#include <hdf5.h>
#include <list>
#include <map>
#include <memory>
#include <string>



#include "daqdataformats/Fragment.hpp"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"

namespace dune {
namespace VDColdboxHDF5Utils {

using dunedaq::daqdataformats::Fragment;
using dunedaq::daqdataformats::FragmentHeader;
typedef std::vector<raw::RawDigit> RawDigits;
typedef std::vector<raw::RDTimeStamp> RDTimeStamps;

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
std::list<std::string> findTopLevelGroupNames(hid_t fd);
std::list<std::string> getTopLevelGroupNames(HDFFileInfoPtr& hdfFileInfoPtr);
std::list<std::string> getMidLevelGroupNames(hid_t gid);
bool attrExists(hid_t object, const std::string& attrname);
hid_t getGroupFromPath(hid_t fd, const std::string &path);

void getHeaderInfo(hid_t the_group, const std::string & det_type,
                   HeaderInfo & info);

typedef std::vector<Fragment> Fragments;
typedef std::map<std::string, std::unique_ptr<Fragments>> FragmentListsByType;
void getFragmentsForEvent(hid_t hdf_file, const std::string & group_name,
                          RawDigits& raw_digits, RDTimeStamps &timestamps);
void getMedianSigma(const raw::RawDigit::ADCvector_t &v_adc, float &median,
                    float &sigma);
}
}
#endif
