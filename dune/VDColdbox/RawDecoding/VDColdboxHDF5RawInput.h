#ifndef VDColdboxHDF5RawInput_h
#define VDColdboxHDF5RawInput_h
#include "art/Framework/Core/InputSourceMacros.h" 
#include "art/Framework/IO/Sources/Source.h" 
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/Core/Frameworkfwd.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Provenance/FileFormatVersion.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune-raw-data/HDFUtils/HDFFileReader.hh"

#include "dune/DuneObj/DUNEHDF5FileInfo.h"


//Is raw a good namespace?
namespace raw {
//Forward declare the class
class VDColdboxHDF5RawInputDetail;
}

class raw::VDColdboxHDF5RawInputDetail {
 public:
  VDColdboxHDF5RawInputDetail(fhicl::ParameterSet const & ps,
                              art::ProductRegistryHelper & rh,
                              art::SourceHelper const & sh)
    : pretend_module_name(ps.get<std::string>("raw_data_label", "daq")),
      pmaker(sh) {
    rh.reconstitutes<raw::DUNEHDF5FileInfo, art::InEvent>(pretend_module_name); 
  };

  void readFile(std::string const & filename, art::FileBlock*& fb) {
    hdf_file_ = dune::HDFFileReader::openFile(filename);
    unprocessedEventList_
        = dune::HDFFileReader::getTopLevelGroupNames(hdf_file_);
    TLOG_INFO("VDColdboxHDF5")
        << "VDColdboxHDF5 opened HDF file with run number " <<
           hdf_file_->runNumber  << " and " <<
           unprocessedEventList_.size() << " events";
    for (const auto & e : unprocessedEventList_)
      TLOG_INFO("VDColdboxHDF5") << e;

    fb = new art::FileBlock(art::FileFormatVersion(1, "RawEvent2011"),
                            filename); 
  };

  bool readNext(art::RunPrincipal const* const inR,
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
    TLOG_INFO("VDColdboxHDF5") << "hi_res_time tv_sec = " << hi_res_time.tv_sec
                               << " tv_nsec = " << hi_res_time.tv_nsec << " (retcode = " << retcode << ")";
    if (retcode == 0) {
      currentTime = ((hi_res_time.tv_sec & 0xffffffff) << 32) | (hi_res_time.tv_nsec & 0xffffffff);
    }
    else {
      TLOG_ERROR("VDColdboxHDF5")
        << "Unable to fetch a high-resolution time with clock_gettime for art::Event Timestamp. "
        << "The art::Event Timestamp will be zero for event ";
    }


    size_t run_id = /*hdf_file_->runNumber*/1; //runNumber can be 0,
                                               //but this seems like an issue
                                               //with art
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

    //Replace with dune hdf 5 info
    put_product_in_principal(std::move(the_info), *outE, pretend_module_name, "");

    return true;
  };

  void closeCurrentFile() {};

 private:
  std::unique_ptr<dune::HDFFileReader::HDFFileInfo> hdf_file_;
  std::list<std::string> unprocessedEventList_;
  std::string pretend_module_name;
  art::SourceHelper const& pmaker;

 };
#endif
