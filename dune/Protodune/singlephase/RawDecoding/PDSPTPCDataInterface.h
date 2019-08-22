////////////////////////////////////////////////////////////////////////
// PDSPTPCDataInterface.h
//
// Tool to unpack RCE and FELIX fragments.  A restructuring from the PDSPTPCRawDecoder module
//
//
/////////////////////////////////////////////////////////////////////////
#ifndef PDSPTPCDataInterface_H
#define PDSPTPCDataInterface_H

#include <vector>

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
//#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/Protodune/singlephase/RawDecoding/data/RDStatus.h"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "PDSPTPCDataInterfaceParent.h"

class PDSPTPCDataInterface : public PDSPTPCDataInterfaceParent {

 public:

  PDSPTPCDataInterface(fhicl::ParameterSet const& ps);

  int retrieveData(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
		   art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, std::vector<raw::RDStatus> &rdstatuses );

  // method to get raw digits, RDTimeStamps, RDStatuses and associations from all input fragments specified by an input label (like "daq:ContainerTPC") but ony for
  // APA's (== crate numbers) on a list.  If the list contains a -1 in it, it returns all APA data found in the input label.

  int retrieveDataAPAListWithLabels(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
				    art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, std::vector<raw::RDStatus> &rdstatuses, 
				    std::vector<int> &apalist);

  // method to get raw digits, RDTimeStamps, RDStatuses and associations for a specified list of APA's.  The list of possible labels on which to find
  // APA data is proved by fcl configuration.

  int retrieveDataForSpecifiedAPAs(art::Event &evt, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
				   art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, std::vector<raw::RDStatus> &rdstatuses,  
				   std::vector<int> &apalist);

  // inputLabel examples:  "daq:TPC" or "daq::ContainerTPC" for RCE, "daq:FELIX" or "daq::ContainerFELIX" for FELIX
  // returns:  0:  success, or   1: discarded corrupted data, or 2: kept some corrupted data

  // same as above, but do not make the associations between raw digits and timestamps

  int retrieveDataForSpecifiedAPAs_NoAssoc(art::Event &evt, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
				           std::vector<raw::RDStatus> &rdstatuses,  
				           std::vector<int> &apalist);

 private:

  bool          _make_tsassocs;
  std::map<int,std::vector<std::string>> _input_labels_by_apa;
  bool          _enforce_same_tick_count;
  bool          _enforce_full_tick_count;
  unsigned int  _full_tick_count;
  bool          _enforce_error_free;
  bool          _enforce_no_duplicate_channels;
  bool          _drop_small_rce_frags;
  size_t        _rce_frag_small_size;
  bool          _rce_drop_frags_with_badcsf;
  bool          _rce_hex_dump;
  bool          _rce_save_frags_to_files;
  bool          _rce_check_buffer_size;
  size_t        _rce_buffer_size_checklimit;

  // flags for attempting to fix FEMB 110's misaligned data

  bool          _rce_fix110;
  unsigned int  _rce_fix110_nticks;

  bool          _felix_hex_dump;
  bool          _felix_drop_frags_with_badcsf;
  bool          _drop_small_felix_frags;
  size_t        _felix_frag_small_size;
  bool          _felix_check_buffer_size;
  size_t        _felix_buffer_size_checklimit;

  unsigned int  _tick_count_this_event; // for use in comparing tick counts for all channels
  bool          _initialized_tick_count_this_event;
  bool          _discard_data;
  bool          _DiscardedCorruptData;
  bool          _KeptCorruptData;

  // some convenience typedefs for porting old code

  typedef std::vector<raw::RawDigit> RawDigits;
  typedef std::vector<raw::RDTimeStamp> RDTimeStamps;
  typedef art::Assns<raw::RawDigit,raw::RDTimeStamp> RDTsAssocs;
  typedef art::PtrMaker<raw::RawDigit> RDPmkr;
  typedef art::PtrMaker<raw::RDTimeStamp> TSPmkr;
  typedef std::vector<raw::RDStatus> RDStatuses;

  // private methods
  bool _processRCE(art::Event &evt, 
		   std::string inputLabel, 
		   RawDigits& raw_digits, 
		   RDTimeStamps &timestamps, 
		   RDTsAssocs &tsassocs, 
		   std::vector<int> &apalist);

  bool _rceProcContNCFrags(art::Handle<artdaq::Fragments> frags, 
			   size_t &n_rce_frags, 
			   bool is_container, 
			   art::Event &evt, 
			   RawDigits& raw_digits, 
			   RDTimeStamps &timestamps, 
			   RDTsAssocs &tsassocs, 
			   std::vector<int> &apalist);

  bool _process_RCE_AUX(art::Event &evt,
			const artdaq::Fragment& frag, 
			RawDigits& raw_digits, 
			RDTimeStamps &timestamps, 
			RDTsAssocs &tsassocs, 
			std::vector<int> &apalist);

  bool _processFELIX(art::Event &evt, 
		     std::string inputLabel, 
		     RawDigits& raw_digits, 
		     RDTimeStamps &timestamps, 
		     RDTsAssocs &tsassocs, 
		     std::vector<int> &apalist);

  bool _felixProcContNCFrags(art::Handle<artdaq::Fragments> frags, 
			     size_t &n_felix_frags, 
			     bool is_container, 
			     art::Event &evt, 
			     RawDigits& raw_digits,
			     RDTimeStamps &timestamps, 
			     RDTsAssocs &tsassocs, 
			     std::vector<int> &apalist);

  bool _process_FELIX_AUX(art::Event &evt,
			  const artdaq::Fragment& frag, 
			  RawDigits& raw_digits, 
			  RDTimeStamps &timestamps, 
			  RDTsAssocs &tsassocs, 
			  std::vector<int> &apalist);

  void computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, 
			  float &median, 
			  float &sigma);

};

DEFINE_ART_CLASS_TOOL(PDSPTPCDataInterface)

#endif
