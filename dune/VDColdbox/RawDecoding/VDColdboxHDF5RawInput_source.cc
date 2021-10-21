#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "dune/VDColdbox/RawDecoding/VDColdboxHDF5RawInput.h"

#include "dune/DuneObj/DUNEHDF5FileInfo.h"
#include "lardataobj/RawData/RDTimeStamp.h"

// Use the standard service interfaces (CatalogInterface and
// FileTransfer) to obtain files.
//namespace art {
//template <raw::VDColdboxHDF5RawInputDetail>
//struct Source_wantFileServices {
//  static constexpr bool value = false
//};
//}

raw::VDColdboxHDF5RawInputDetail::VDColdboxHDF5RawInputDetail(
    fhicl::ParameterSet const & ps,
    art::ProductRegistryHelper & rh,
    art::SourceHelper const & sh) 
  : pretend_module_name(ps.get<std::string>("raw_data_label", "daq")),
    pmaker(sh) {
  rh.reconstitutes<raw::DUNEHDF5FileInfo, art::InEvent>(pretend_module_name); 
  rh.reconstitutes<raw::RDTimeStamp, art::InEvent>(pretend_module_name, "trigger");
}

void raw::VDColdboxHDF5RawInputDetail::readFile(
    std::string const & filename, art::FileBlock*& fb) {
  hdf_file_ = dune::VDColdboxHDF5Utils::openFile(filename);
  unprocessedEventList_
      = dune::VDColdboxHDF5Utils::getTopLevelGroupNames(hdf_file_);
  MF_LOG_INFO("VDColdboxHDF5")
      << "VDColdboxHDF5 opened HDF file with run number " <<
         hdf_file_->runNumber  << " and " <<
         unprocessedEventList_.size() << " events";
  for (const auto & e : unprocessedEventList_)
    MF_LOG_INFO("VDColdboxHDF5") << e;

  fb = new art::FileBlock(art::FileFormatVersion(1, "RawEvent2011"),
                          filename); 
}

bool raw::VDColdboxHDF5RawInputDetail::readNext(
    art::RunPrincipal const* const inR,
    art::SubRunPrincipal const* const inSR,
    art::RunPrincipal*& outR,
    art::SubRunPrincipal*& outSR,
    art::EventPrincipal*& outE) {
  // Establish default 'results'
  outR = 0;
  outSR = 0;
  outE = 0;

  if (unprocessedEventList_.empty()) {return false;}
  std::string nextEventGroupName = unprocessedEventList_.front();
  unprocessedEventList_.pop_front();

  art::Timestamp currentTime = 0;
  timespec hi_res_time;
  int retcode = clock_gettime(CLOCK_REALTIME, &hi_res_time);
  MF_LOG_INFO("VDColdboxHDF5") << "hi_res_time tv_sec = " << hi_res_time.tv_sec
                             << " tv_nsec = " << hi_res_time.tv_nsec << " (retcode = " << retcode << ")";
  if (retcode == 0) {
    currentTime = ((hi_res_time.tv_sec & 0xffffffff) << 32) | (hi_res_time.tv_nsec & 0xffffffff);
  }
  else {
    MF_LOG_ERROR("VDColdboxHDF5")
      << "Unable to fetch a high-resolution time with clock_gettime for art::Event Timestamp. "
      << "The art::Event Timestamp will be zero for event ";
  }


  size_t run_id = -1; //runNumber can be 0,
                      //but this seems like an issue
                      //with art

  //Accessing run number
  hid_t the_group = dune::VDColdboxHDF5Utils::getGroupFromPath(
      hdf_file_->filePtr, nextEventGroupName);
  std::list<std::string> detector_types =
      dune::VDColdboxHDF5Utils::getMidLevelGroupNames(the_group);
  dune::VDColdboxHDF5Utils::HeaderInfo header_info;
  std::string det_type = "TriggerRecordHeader";
  dune::VDColdboxHDF5Utils::getHeaderInfo(the_group, det_type, header_info);
  std::cout << "   Magic word: 0x" << std::hex << header_info.magicWord <<
               std::dec << std::endl;
  std::cout << "   Version: " << std::dec << header_info.version <<
               std::dec << std::endl;
  std::cout << "   Trig Num: " << std::dec << header_info.trigNum <<
               std::dec << std::endl;
  std::cout << "   Trig Timestamp: " << std::dec <<
               header_info.trigTimestamp << std::dec << std::endl;
  std::cout << "   No. of requested components:   " << std::dec <<
               header_info.nReq << std::dec << std::endl;
  std::cout << "   Run Number: " << std::dec << header_info.runNum <<
               std::endl;
  std::cout << "   Error bits: " << std::dec << header_info.errBits <<
               std::endl;
  std::cout << "   Trigger type: " << std::dec << header_info.triggerType <<
               std::endl;

  run_id = header_info.runNum;
  std::unique_ptr<raw::RDTimeStamp> rd_timestamp(
    new raw::RDTimeStamp(header_info.trigTimestamp));

  // make new run if inR is 0 or if the run has changed
  if (inR == 0 || inR->run() != run_id) {
    outR = pmaker.makeRunPrincipal(run_id, currentTime);
  }

  // make new subrun if inSR is 0 or if the subrun has changed
  art::SubRunID subrun_check(run_id, 1);
  if (inSR == 0 || subrun_check != inSR->subRunID()) {
    outSR = pmaker.makeSubRunPrincipal(run_id, 1, currentTime);
  }

  //Where to get event number?
  //For now -- parse from input record
  std::string event_str = nextEventGroupName;
  std::string trig = "TriggerRecord";
  auto pos = event_str.begin() + event_str.find(trig);
  event_str.erase(pos, pos + trig.size());
  int event = std::stoi(event_str);
  outE = pmaker.makeEventPrincipal(run_id, 1, event, currentTime);

  
  std::unique_ptr<DUNEHDF5FileInfo> the_info(
      new DUNEHDF5FileInfo(hdf_file_->fileName, hdf_file_->filePtr,
                           0, nextEventGroupName));

  put_product_in_principal(std::move(the_info), *outE, pretend_module_name,
                           "");
  put_product_in_principal(std::move(rd_timestamp), *outE, pretend_module_name,
                           "trigger");

  return true;
}

//typedef for shorthand
namespace raw {
using VDColdboxHDF5RawInputSource = art::Source<VDColdboxHDF5RawInputDetail>;
}


DEFINE_ART_INPUT_SOURCE(raw::VDColdboxHDF5RawInputSource)
