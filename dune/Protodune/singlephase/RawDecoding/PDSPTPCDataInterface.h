////////////////////////////////////////////////////////////////////////
// PDSPTPCDataInterface.h
//
// Tool to unpack RCE and FELIX fragments.  A restructuring from the PDSPTPCRawDecoder module
//
// These methods take references to vectors of raw::RawDigit, raw::RDTimeStamp, and raw::RDStatus data products as arguments.
// These vectors are not cleared on input, and so when data are retrieved, they are appended to any existing data already
// in the vectors.  The RDStatus vector is an exception, where just one RDStatus instance will be in the vector.  Previously
// accumulated RDStatus values from previous calls will be logically ORed into the single RDStatus instances contents.
//
//  Methods are provided to retrieve all data from fragments on an input label, or by specified APA list.  In cases where
//  data from a specified APA are requested but no labels are provided by the caller, labels are input via FCL parameters.
//  This is true because data from an APA may appear with different labels during the course of the ProtoDUNE-SP run.
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
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "dune/DuneObj/PDSPTPCDataInterfaceParent.h"

class PDSPTPCDataInterface : public PDSPTPCDataInterfaceParent {

 public:

  PDSPTPCDataInterface(fhicl::ParameterSet const& ps);

  int retrieveData(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
		   std::vector<raw::RDStatus> &rdstatuses );

  // method to get raw digits, RDTimeStamps, RDStatuses from all input fragments specified by an input label (like "daq:ContainerTPC") but ony for
  // APA's (== crate numbers) on a list.  If the list contains a -1 in it, it returns all APA data found in the input label.

  int retrieveDataAPAListWithLabels(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
				    std::vector<raw::RDStatus> &rdstatuses, 
				    std::vector<int> &apalist);

  // method to get raw digits, RDTimeStamps, RDStatuses for a specified list of APA's.  The list of possible labels on which to find
  // APA data is proved by fcl configuration.

  int retrieveDataForSpecifiedAPAs(art::Event &evt, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
				   std::vector<raw::RDStatus> &rdstatuses,  
				   std::vector<int> &apalist);

  // inputLabel examples:  "daq:TPC" or "daq::ContainerTPC" for RCE, "daq:FELIX" or "daq::ContainerFELIX" for FELIX
  // returns:  0:  success, or   1: discarded corrupted data, or 2: kept some corrupted data

 private:

  std::map<int,std::vector<std::string>> _input_labels_by_apa;

  // what to do with unexpected crate numbers

  unsigned int  _default_crate_if_unexpected;

  long int _min_offline_channel;  // min offline channel to decode.  <0: no limit
  long int _max_offline_channel;  // max offline channel to decode.  <0: no limit.  max<min: no limit

  bool          _enforce_same_tick_count;
  bool          _enforce_median_tick_count;
  bool          _enforce_full_tick_count;
  unsigned int  _full_tick_count;
  bool          _enforce_error_free;
  bool          _enforce_no_duplicate_channels;
  bool          _drop_small_rce_frags;
  size_t        _rce_frag_small_size;
  bool          _rce_drop_frags_with_badsf;
  bool          _rce_drop_frags_with_badc;
  bool          _rce_hex_dump;
  bool          _rce_save_frags_to_files;
  bool          _rce_check_buffer_size;
  size_t        _rce_buffer_size_checklimit;

  // flags for attempting to fix FEMB 110's misaligned data

  bool          _rce_fix110;
  unsigned int  _rce_fix110_nticks;

  bool          _felix_hex_dump;
  bool          _felix_drop_frags_with_badsf;
  bool          _felix_drop_frags_with_badc;
  bool          _drop_small_felix_frags;
  size_t        _felix_frag_small_size;
  bool          _felix_check_buffer_size;
  size_t        _felix_buffer_size_checklimit;

  unsigned int  _tick_count_this_event; // for use in comparing tick counts for all channels
  bool          _initialized_tick_count_this_event;
  bool          _DiscardedCorruptData;
  bool          _KeptCorruptData;

  // some convenience typedefs for porting old code

  typedef std::vector<raw::RawDigit> RawDigits;
  typedef std::vector<raw::RDTimeStamp> RDTimeStamps;
  typedef std::vector<raw::RDStatus> RDStatuses;

  // private methods

  void _collectRDStatus(std::vector<raw::RDStatus> &rdstatuses);

  bool _processRCE(art::Event &evt, 
		   std::string inputLabel, 
		   RawDigits& raw_digits, 
		   RDTimeStamps &timestamps, 
		   std::vector<int> &apalist);

  bool _rceProcContNCFrags(art::Handle<artdaq::Fragments> frags, 
			   size_t &n_rce_frags, 
			   bool is_container, 
			   art::Event &evt, 
			   RawDigits& raw_digits, 
			   RDTimeStamps &timestamps, 
			   std::vector<int> &apalist);

  bool _process_RCE_AUX(art::Event &evt,
			const artdaq::Fragment& frag, 
			RawDigits& raw_digits, 
			RDTimeStamps &timestamps, 
			std::vector<int> &apalist);

  bool _processFELIX(art::Event &evt, 
		     std::string inputLabel, 
		     RawDigits& raw_digits, 
		     RDTimeStamps &timestamps, 
		     std::vector<int> &apalist);

  bool _felixProcContNCFrags(art::Handle<artdaq::Fragments> frags, 
			     size_t &n_felix_frags, 
			     bool is_container, 
			     art::Event &evt, 
			     RawDigits& raw_digits,
			     RDTimeStamps &timestamps, 
			     std::vector<int> &apalist);

  bool _process_FELIX_AUX(art::Event &evt,
			  const artdaq::Fragment& frag, 
			  RawDigits& raw_digits, 
			  RDTimeStamps &timestamps, 
			  std::vector<int> &apalist);

  void computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, 
			  float &median, 
			  float &sigma);

};


#endif
