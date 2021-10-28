////////////////////////////////////////////////////////////////////////
// Class:       PDSPTPCRawDecoder
// Plugin Type: producer (art v2_10_03)
// File:        PDSPTPCRawDecoder_module.cc
//
// Generated at Fri Mar  2 15:36:20 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_02_00.
// Original code from Jingbo Wang in separate RCE and FELIX raw decoders
// *********************************************************************
// July 2018, Maggie Greenwood, Added histograms for error checking.
// March 2019 restructure to read in lists of module label/instances of inputs, including
//  an optional GetManyByType if we don't know in advance what labels we're going to see
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art_root_io/TFileService.h"

#include <memory>
#include <cmath>

// ROOT includes
#include "TH1.h"
#include "TStyle.h"
#include "TMath.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/RceFragment.hh"
#include "dune-raw-data/Overlays/FelixFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "dune-raw-data/Overlays/FragmentType.hh"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "dam/HeaderFragmentUnpack.hh"
#include "dam/DataFragmentUnpack.hh"
#include "dam/TpcFragmentUnpack.hh"
#include "dam/TpcStreamUnpack.hh"
#include "dam/access/WibFrame.hh"
#include "dam/access/Headers.hh"
#include "dam/access/TpcStream.hh"
#include "dam/access/TpcRanges.hh"
#include "dam/access/TpcToc.hh"
#include "dam/access/TpcPacket.hh"
#include "dam/RceFragmentUnpack.hh"

// larsoft includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/raw.h"

// DUNE includes
#include "dune/DuneObj/RDStatus.h"

class PDSPTPCRawDecoder;

class PDSPTPCRawDecoder : public art::EDProducer {

public:
  explicit PDSPTPCRawDecoder(fhicl::ParameterSet const & p);
  PDSPTPCRawDecoder(PDSPTPCRawDecoder const &) = delete;
  PDSPTPCRawDecoder(PDSPTPCRawDecoder &&) = delete;
  PDSPTPCRawDecoder & operator = (PDSPTPCRawDecoder const &) = delete;
  PDSPTPCRawDecoder & operator = (PDSPTPCRawDecoder &&) = delete;
  void produce(art::Event & e) override;

private:
  typedef std::vector<raw::RawDigit> RawDigits;
  typedef std::vector<raw::RDTimeStamp> RDTimeStamps;
  typedef art::Assns<raw::RawDigit,raw::RDTimeStamp> RDTsAssocs;
  typedef art::PtrMaker<raw::RawDigit> RDPmkr;
  typedef art::PtrMaker<raw::RDTimeStamp> TSPmkr;
  typedef std::vector<raw::RDStatus> RDStatuses;

  // configuration parameters

  std::vector<int> _apas_to_decode;

  long int _min_offline_channel;  // min offline channel to decode.  <0: no limit
  long int _max_offline_channel;  // max offline channel to decode.  <0: no limit.  max<min: no limit

  bool _rce_useInputLabels;
  bool _felix_useInputLabels;
  std::vector<std::string>   _rce_input_labels;   // input labels also include instances.  Example: "daq:TPC" or "daq::ContainerTPC"
  std::vector<std::string>   _felix_input_labels; // input labels also include instances.  Example: "daq:FELIX" or "daq::ContainerFELIX"
  bool          _rce_enforce_fragment_type_match;
  bool          _felix_enforce_fragment_type_match;
  int           _rce_fragment_type;
  int           _felix_fragment_type;
  std::string   _output_label;
  bool          _enforce_full_channel_count;
  unsigned int  _full_channel_count;
  bool          _enforce_same_tick_count;
  bool          _enforce_full_tick_count;
  unsigned int  _full_tick_count;
  bool          _enforce_error_free;
  bool          _enforce_no_duplicate_channels;
  bool          _drop_events_with_small_rce_frags;
  bool          _drop_small_rce_frags;
  size_t        _rce_frag_small_size;
  bool          _rce_drop_frags_with_badsf;
  bool          _rce_drop_frags_with_badc;
  bool          _rce_hex_dump;
  bool          _rce_save_frags_to_files;
  bool          _rce_check_buffer_size;
  size_t        _rce_buffer_size_checklimit;

  // what to do with unexpected crate numbers

  unsigned int  _default_crate_if_unexpected;

  // flags for attempting to fix FEMB 110's misaligned data

  bool          _rce_fix110;
  unsigned int  _rce_fix110_nticks;

  bool          _felix_hex_dump;
  bool          _felix_drop_frags_with_badsf;
  bool          _felix_drop_frags_with_badc;
  bool          _drop_events_with_small_felix_frags;
  bool          _drop_small_felix_frags;
  size_t        _felix_frag_small_size;
  bool          _felix_check_buffer_size;
  size_t        _felix_buffer_size_checklimit;

  bool          _compress_Huffman;
  bool          _print_coldata_convert_count;

  //declare histogram data memebers
  bool	_make_histograms;
  unsigned int 	duplicate_channels;
  unsigned int 	error_counter;
  unsigned int incorrect_ticks;
  unsigned int rcechans;
  unsigned int felixchans;
  TH1D * fIncorrectTickNumbers;
  //TH1I * fIncorrectTickNumbersZoomed;
  TH1I * fParticipRCE;
  TH1I * fParticipFELIX;
  TH1I * fDuplicatesNumber;
  TH1D * fErrorsNumber;
  TH1I * fFragSizeRCE;
  TH1I * fFragSizeFELIX;

  // flags and state needed for the data integrity enforcement mechanisms

  unsigned int  _tick_count_this_event; // for use in comparing tick counts for all channels
  bool          _initialized_tick_count_this_event;
  bool          _discard_data;
  const unsigned int  _duplicate_channel_checklist_size=15360;
  bool          _duplicate_channel_checklist[15360]; 
  bool          _DiscardedCorruptData;
  bool          _KeptCorruptData;

  // internal methods

  bool _processRCE(art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm);
  bool _rceProcContNCFrags(art::Handle<artdaq::Fragments> frags, size_t &n_rce_frags, bool is_container, 
			   art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm);
  bool _process_RCE_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm, size_t ntickscheck);
  void _process_RCE_nticksvf(const artdaq::Fragment& frag, std::vector<size_t> &nticksvec);

  bool _processFELIX(art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm);
  bool _felixProcContNCFrags(art::Handle<artdaq::Fragments> frags, size_t &n_felix_frags, bool is_container, art::Event &evt, RawDigits& raw_digits,
			     RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm);
  bool _process_FELIX_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm);

  void computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma);

  std::vector<int16_t> _buffer;
};


PDSPTPCRawDecoder::PDSPTPCRawDecoder(fhicl::ParameterSet const & p) : EDProducer{p}
{
  std::vector<int> emptyivec;
  _apas_to_decode = p.get<std::vector<int> >("APAsToDecode",emptyivec);
  _default_crate_if_unexpected = p.get<unsigned int>("DefaultCrateIfUnexpected",3);

  _rce_input_labels.resize(0);
  _rce_useInputLabels =  p.get_if_present<std::vector<std::string> >("RCERawDataLabels",_rce_input_labels);
  if (!_rce_useInputLabels)
    {
      _rce_input_labels.resize(0);
    }
  if (_rce_input_labels.size() == 0)  // handle case where array is in the fcl document, just is empty
    {
      _rce_useInputLabels = false;
    }
  _rce_enforce_fragment_type_match = p.get<bool>("RCEEnforceFragmentTypeMatch",false);
  _rce_fragment_type = p.get<int>("RCEFragmentType",2);
  _drop_events_with_small_rce_frags = p.get<bool>("RCEDropEventsWithSmallFrags",false);
  _drop_small_rce_frags = p.get<bool>("RCEDropSmallFrags",true);
  _rce_frag_small_size = p.get<unsigned int>("RCESmallFragSize",10000);
  _rce_drop_frags_with_badsf = p.get<bool>("RCEDropFragsWithBadSF",true);
  _rce_drop_frags_with_badc = p.get<bool>("RCEDropFragsWithBadC",true);
  _rce_hex_dump = p.get<bool>("RCEHexDump",false);  
  _rce_save_frags_to_files = p.get<bool>("RCESaveFragsToFiles",false);  
  _rce_check_buffer_size = p.get<bool>("RCECheckBufferSize",true);
  _rce_buffer_size_checklimit = p.get<unsigned int>("RCEBufferSizeCheckLimit",10000000);

  // parameters to steer the FEMB 110 band-aid

  _rce_fix110 = p.get<bool>("RCEFIX110",true);
  _rce_fix110_nticks = p.get<unsigned int>("RCEFIX110NTICKS",18);

  _felix_input_labels.resize(0);
  _felix_useInputLabels =  p.get_if_present<std::vector<std::string> >("FELIXRawDataLabels",_felix_input_labels);
  if (!_felix_useInputLabels)
    {
      _felix_input_labels.resize(0);
    }
  if (_felix_input_labels.size() == 0)  // handle case where array is in the fcl document, just is empty
    {
      _felix_useInputLabels = false;
    }

  _felix_enforce_fragment_type_match = p.get<bool>("FELIXEnforceFragmentTypeMatch",false);
  _felix_fragment_type = p.get<int>("FELIXFragmentType",8);
  _felix_drop_frags_with_badsf = p.get<bool>("FELIXDropFragsWithBadSF",true);
  _felix_drop_frags_with_badc = p.get<bool>("FELIXDropFragsWithBadC",true);
  _felix_hex_dump = p.get<bool>("FELIXHexDump",false);  
  _drop_events_with_small_felix_frags = p.get<bool>("FELIXDropEventsWithSmallFrags",false);
  _drop_small_felix_frags = p.get<bool>("FELIXDropSmallFrags",true);
  _felix_frag_small_size = p.get<unsigned int>("FELIXSmallFragSize",10000);
  _felix_check_buffer_size = p.get<bool>("FELIXCheckBufferSize",true);
  _felix_buffer_size_checklimit = p.get<unsigned int>("FELIXBufferSizeCheckLimit",10000000);

  _output_label = p.get<std::string>("OutputDataLabel");

  _enforce_full_channel_count = p.get<bool>("EnforceFullChannelCount", false);
  _full_channel_count = p.get<unsigned int>("FullChannelCount", 15360);
  _enforce_same_tick_count = p.get<bool>("EnforceSameTickCount",false);
  _enforce_full_tick_count = p.get<bool>("EnforceFullTickCount",false);
  _full_tick_count = p.get<unsigned int>("FullTickCount",6000);
  _enforce_error_free = p.get<bool>("EnforceErrorFree",false);
  _enforce_no_duplicate_channels = p.get<bool>("EnforceNoDuplicateChannels", true);

  _compress_Huffman = p.get<bool>("CompressHuffman",false);
  _print_coldata_convert_count = p.get<bool>("PrintColdataConvertCount",false);

  _min_offline_channel = p.get<long int>("MinOfflineChannel",-1);
  _max_offline_channel = p.get<long int>("MaxOfflineChannel",-1);

  produces<RawDigits>( _output_label ); //the strings in <> are the typedefs defined above
  produces<RDTimeStamps>( _output_label );
  produces<RDTsAssocs>( _output_label );
  produces<RDStatuses>( _output_label );

  //Initialize Histograms if the tag is present
  //art::ServiceHandle<art::TFileService> fs;
  _make_histograms = p.get<bool>("MakeHistograms",false);

  if (_make_histograms)
    {
      art::ServiceHandle<art::TFileService> tFileService;

      //Number of channels with wrong number of tics plotted to have an adjusted log2 scale on x axis
      fIncorrectTickNumbers = tFileService->make<TH1D>("fIncorrectTickNumbers","Channels with Unexpected Number of Ticks",  45, -0.5, 14.5);
      fIncorrectTickNumbers->GetXaxis()->SetTitle("Channels with an Unexpected Number of Ticks");
      fIncorrectTickNumbers->GetXaxis()->SetBinLabel(2,"1");
      fIncorrectTickNumbers->GetXaxis()->SetBinLabel(8,"4");
      fIncorrectTickNumbers->GetXaxis()->SetBinLabel(14,"16");
      fIncorrectTickNumbers->GetXaxis()->SetBinLabel(23,"128");
      fIncorrectTickNumbers->GetXaxis()->SetBinLabel(32,"1024");
      fIncorrectTickNumbers->GetXaxis()->SetBinLabel(38,"4096");
      fIncorrectTickNumbers->GetXaxis()->SetBinLabel(44, "16384");
	
      //same as fIncorrectTickNumbers but with a zoomed domain
      //fIncorrectTickNumbersZoomed = tFileService->make<TH1I>("fIncorrectTickNumbersZoomed","Channels with Unexpected Number of Ticks", , 0.5, 100.5);
      //fIncorrectTickNumbersZoomed->GetXaxis()->SetTitle("Channels with an Unexpected Number of Tics (log2)");

      //number of participating RCE channels per event
      fParticipRCE = tFileService->make<TH1I>("fParticipRCE","Participating RCE channels", 130, 0.5, 15000.5); //expected value 128000
      fParticipRCE->GetXaxis()->SetTitle("RCE channels");

      //number of participating FELIX channels per event
      fParticipFELIX = tFileService->make<TH1I>("fParticipFELIX","Participating FELIX channels", 100, 0.5, 3000.5); //expected value 2560
      fParticipFELIX->GetXaxis()->SetTitle("FELIX channels");

      //number of duplicated channels per event
      fDuplicatesNumber = tFileService->make<TH1I>("fDuplicatesNumber", "Number of Duplucated Channels", 200, 0.5, 200.5);
      fDuplicatesNumber->GetXaxis()->SetTitle("Number of Duplicates");
      //gStyle->SetOptStat("nemro");

      //number of channels with error returns
      fErrorsNumber = tFileService->make<TH1D>("fErrorsNumber", "Channels with Errors", 45, -0.5, 14.5);
      fErrorsNumber->GetXaxis()->SetTitle("Number of channels with errors");
      fErrorsNumber->GetXaxis()->SetBinLabel(2,"1");
      fErrorsNumber->GetXaxis()->SetBinLabel(8,"4");
      fErrorsNumber->GetXaxis()->SetBinLabel(14,"16");
      fErrorsNumber->GetXaxis()->SetBinLabel(23,"128");
      fErrorsNumber->GetXaxis()->SetBinLabel(32,"1024");
      fErrorsNumber->GetXaxis()->SetBinLabel(38,"4096");
      fErrorsNumber->GetXaxis()->SetBinLabel(44, "16384");

      //total fragment sizes
      fFragSizeRCE = tFileService->make<TH1I>("fFragSizeRCE", "RCE Fragment Size", 100, 0.5, 288000000.5);
      fFragSizeRCE->GetXaxis()->SetTitle("Size of RCE Fragments (bytes)");

      fFragSizeFELIX = tFileService->make<TH1I>("fFragSizeFELIX", "FELIX Fragment Size", 100, 0.5, 57600000.5);
      fFragSizeFELIX->GetXaxis()->SetTitle("Size of FELIX Fragments (bytes)");
    }
}

void PDSPTPCRawDecoder::produce(art::Event &e)
{
  RawDigits raw_digits;
  RDTimeStamps rd_timestamps;
  RDTsAssocs rd_ts_assocs;

  RDPmkr rdpm(e,_output_label);
  TSPmkr tspm(e,_output_label);

  error_counter = 0; //reset the errors to zero for each run
  incorrect_ticks = 0;
  duplicate_channels = 0;
  rcechans = 0;
  felixchans = 0;

  _initialized_tick_count_this_event = false;
  _discard_data = false;   // true if we're going to drop the whole event's worth of data
  _DiscardedCorruptData = false;   // can be set to true if we drop some of the event's data
  _KeptCorruptData = false;      // true if we identify a corruption candidate but are skipping the test to drop it

  for (size_t i=0; i<_duplicate_channel_checklist_size; ++i) _duplicate_channel_checklist[i] = false;
  
  _processRCE(e,raw_digits,rd_timestamps,rd_ts_assocs,rdpm,tspm);
  _processFELIX(e,raw_digits,rd_timestamps,rd_ts_assocs,rdpm,tspm);

  //Make the histograms for error checking. (other histograms are filled within the _process and _AUX functions)
  if(_make_histograms)
    {
      fErrorsNumber->Fill(log2(error_counter));
      fDuplicatesNumber->Fill(duplicate_channels);
      fIncorrectTickNumbers->Fill(log2(incorrect_ticks));
      fParticipFELIX->Fill(felixchans);
      fParticipRCE->Fill(rcechans);
      //fIncorrectTickNumbersZoomed->Fill(incorrect_ticks);
    }

  if (_enforce_full_channel_count && raw_digits.size() != _full_channel_count) 
    {
      MF_LOG_WARNING("PDSPTPCRawDecoder:") << "Wrong Total number of Channels " << raw_digits.size()  
					   << " which is not " << _full_channel_count << ". Discarding Data";
      _DiscardedCorruptData = true;
      _discard_data = true;
    }

  if (_discard_data)
    {
      RawDigits empty_raw_digits;
      RDTimeStamps empty_rd_timestamps;
      RDTsAssocs empty_rd_ts_assocs;
      RDStatuses statuses;
      statuses.emplace_back(true,false,1);
      e.put(std::make_unique<decltype(empty_raw_digits)>(std::move(empty_raw_digits)),_output_label);
      e.put(std::make_unique<decltype(empty_rd_timestamps)>(std::move(empty_rd_timestamps)),_output_label);
      e.put(std::make_unique<decltype(empty_rd_ts_assocs)>(std::move(empty_rd_ts_assocs)),_output_label);
      e.put(std::make_unique<decltype(statuses)>(std::move(statuses)),_output_label);
    }
  else
    {
      RDStatuses statuses;
      unsigned int statword=0;
      if (_DiscardedCorruptData) statword |= 1;
      if (_KeptCorruptData) statword |= 2;
      statuses.emplace_back(_DiscardedCorruptData,_KeptCorruptData,statword);
      e.put(std::make_unique<decltype(raw_digits)>(std::move(raw_digits)),_output_label);
      e.put(std::make_unique<decltype(rd_timestamps)>(std::move(rd_timestamps)),_output_label);
      e.put(std::make_unique<decltype(rd_ts_assocs)>(std::move(rd_ts_assocs)),_output_label);
      e.put(std::make_unique<decltype(statuses)>(std::move(statuses)),_output_label);
    }
}

bool PDSPTPCRawDecoder::_processRCE(art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm)
{
  size_t n_rce_frags = 0;
  bool have_data=false;
  bool have_data_nc=false;

  if (_rce_useInputLabels)
    {
      for (size_t ilabel = 0; ilabel < _rce_input_labels.size(); ++ ilabel)
	{
	  if (_rce_input_labels.at(ilabel).find("Container") != std::string::npos)
	    {
              auto cont_frags = evt.getHandle<artdaq::Fragments>(_rce_input_labels.at(ilabel));
	      if (cont_frags)
		{
		  have_data = true;
	          if (! _rceProcContNCFrags(cont_frags, n_rce_frags, true, evt, raw_digits, timestamps, tsassocs, rdpm, tspm))
		    {
		      return false;
		    }
		}
	    }
	  else
	    {
	      auto frags = evt.getHandle<artdaq::Fragments>(_rce_input_labels.at(ilabel));
	      if (frags)
		{
		  have_data_nc = true;
	          if (! _rceProcContNCFrags(frags, n_rce_frags, false, evt, raw_digits, timestamps, tsassocs, rdpm, tspm))
		    {
		      return false;
		    }
		}
	    }
	}
    }
  else  // get all the fragments in the event and look for the ones that say TPC in them
    {
      //std::vector<art::Handle<artdaq::Fragments> > fraghv;  // fragment handle vector
      //evt.getManyByType(fraghv);
      auto fraghv = evt.getMany<artdaq::Fragments>();

      for (size_t ihandle=0; ihandle<fraghv.size(); ++ihandle)
	{
	  if (fraghv.at(ihandle).provenance()->inputTag().instance().find("TPC") != std::string::npos)
	    {
	      if (fraghv.at(ihandle).isValid())
		{
	          if (fraghv.at(ihandle).provenance()->inputTag().instance().find("Container") != std::string::npos)
		    {
		      have_data = true;
		      if (! _rceProcContNCFrags(fraghv.at(ihandle), n_rce_frags, true, evt, raw_digits, timestamps, tsassocs, rdpm, tspm) )
			{
			  return false;
			}
		    }
		  else
		    {
		      have_data_nc = true;
		      if (! _rceProcContNCFrags(fraghv.at(ihandle), n_rce_frags, false, evt, raw_digits, timestamps, tsassocs, rdpm, tspm))
			{
			  return false;
			}
		    }
		}
	    }
	}
    }

  //MF_LOG_INFO("_processRCE")
  //<< " Processed " << n_rce_frags
  //<< " RCE Fragments, making "
  //<< raw_digits.size()
  //<< " RawDigits.";

  // returns true if we want to add to the number of fragments processed.  Separate flag used
  // for data error conditions (_discard_data).

  return have_data || have_data_nc;
}

bool PDSPTPCRawDecoder::_rceProcContNCFrags(art::Handle<artdaq::Fragments> frags, size_t &n_rce_frags, bool is_container, 
					    art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm)
{
  //size of RCE fragments into histogram
  if(_make_histograms)
    {
      size_t rcebytes = 0;
      for (auto const& frag : *frags)
	{
	  rcebytes = rcebytes + (frag.sizeBytes());
	}
      fFragSizeRCE->Fill(rcebytes);
    }
    

  // figure out what the median number of ticks is

  std::vector<size_t> nticksvec;
  for (auto const& frag : *frags)
    {
      //std::cout << "RCE fragment size bytes: " << frag.sizeBytes() << std::endl; 

      // skip small fragments even here

      if (frag.sizeBytes() >= _rce_frag_small_size || (!_drop_small_rce_frags && !_drop_events_with_small_rce_frags))
	{
	  if (is_container)
	    {
	      artdaq::ContainerFragment cont_frag(frag);
	      for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
		{
		  _process_RCE_nticksvf(*cont_frag[ii], nticksvec);
		}
	    }
	  else
	    {
	      _process_RCE_nticksvf(frag, nticksvec);
	    }
	}
    }
  if (nticksvec.size() == 0)
    {
      MF_LOG_WARNING("_process_RCE:") << " No valid nticks to check.  Discarding Event.";
      _discard_data = true; 
      _DiscardedCorruptData = true;
      evt.removeCachedProduct(frags);
      return false;
    }
  size_t nticksmedian = TMath::Median(nticksvec.size(),nticksvec.data()) + 0.01;  // returns a double -- want to make sure it gets truncated to the right integer


  // actually process the 
  for (auto const& frag : *frags)
    {
      //std::cout << "RCE fragment size bytes: " << frag.sizeBytes() << std::endl; 

      bool process_flag = true;
      if (frag.sizeBytes() < _rce_frag_small_size)
	{
	  if ( _drop_events_with_small_rce_frags )
	    { 
	      MF_LOG_WARNING("_process_RCE:") << " Small RCE fragment size: " << frag.sizeBytes() << " Discarding Event on request.";
	      _discard_data = true; 
	      _DiscardedCorruptData = true;
	      evt.removeCachedProduct(frags);
	      return false;
	    }
	  if ( _drop_small_rce_frags )
	    { 
	      MF_LOG_WARNING("_process_RCE:") << " Small RCE fragment size: " << frag.sizeBytes() << " Discarding just this fragment on request.";
	      _DiscardedCorruptData = true;
	      process_flag = false;
	    }
	  _KeptCorruptData = true;
	}
      if (process_flag)
	{
	  if (is_container)
	    {
	      artdaq::ContainerFragment cont_frag(frag);
	      for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
		{
		  if (_process_RCE_AUX(*cont_frag[ii], raw_digits, timestamps, tsassocs, rdpm, tspm, nticksmedian)) ++n_rce_frags;
		}
	    }
	  else
	    {
	      if (_process_RCE_AUX(frag, raw_digits, timestamps,tsassocs, rdpm, tspm, nticksmedian)) ++n_rce_frags;
	    }
	}
    }
  evt.removeCachedProduct(frags);
  return true;
}

void PDSPTPCRawDecoder::_process_RCE_nticksvf(
					      const artdaq::Fragment& frag, 
					      std::vector<size_t> &nticksvec
					      )
{

  if (_rce_hex_dump)
    {
      std::ios oldState(nullptr);
      oldState.copyfmt(std::cout);

      std::cout << "RCE Fragment: all numbers in hex "  << std::hex
		<< "   SequenceID = " << frag.sequenceID()
		<< "   fragmentID = " << frag.fragmentID()
		<< "   fragmentType = " << (unsigned)frag.type()
		<< "   Timestamp =  " << frag.timestamp() << std::endl;
      std::cout << "Offset      Data";
      artdaq::Fragment fragloc(frag);
      unsigned char *dbegin = reinterpret_cast<unsigned char *>(fragloc.dataAddress());
      size_t dsize = fragloc.dataSizeBytes();
      size_t offcounter=0;
      for (size_t bcounter=0; bcounter<dsize;++bcounter)
	{
	  if ( (offcounter % 8) == 0 )
	    {
	      std::cout << std::endl << std::hex << std::setfill('0') << std::setw(8) << offcounter << " ";
	    }
	  std::cout << std::hex << std::setfill('0') << std::setw(2) << (int) *dbegin << " ";
	  dbegin++;
	  offcounter++;
	}
      std::cout << std::endl;
      std::cout.copyfmt(oldState);
    }

  dune::RceFragment rce(frag);

  if (_rce_save_frags_to_files)
    {
      TString outfilename="rce";
      outfilename += frag.sequenceID();
      outfilename += "_";
      outfilename += frag.fragmentID();
      outfilename+=".fragment";
      rce.save(outfilename.Data());
      std::cout << "Saved an RCE fragment with " << rce.size() << " streams: " << outfilename << std::endl;
    }

  artdaq::Fragment cfragloc(frag);
  size_t cdsize = cfragloc.dataSizeBytes();
  const uint64_t* cdptr = (uint64_t const*) (cfragloc.dataBeginBytes() + 12);  // see dune-raw-data/Overlays/RceFragment.cc
  HeaderFragmentUnpack const cdheader(cdptr);
  //bool isOkay = RceFragmentUnpack::isOkay(cdptr,cdsize+sizeof(cdheader));
  if (cdsize>16) cdsize -= 16;
  bool isOkay = RceFragmentUnpack::isOkay(cdptr,cdsize);
  if (!isOkay) return;

  // // skip if damaged but not FEMB 302
  // if (cdheader.isData())
  //   {
  //     DataFragmentUnpack df(cdptr);
  //     if (df.isTpcDamaged())
  // 	{
  // 	  bool found302 = false;
  // 	   for (int i = 0; i < rce.size(); ++i)
  //            {
  //              auto const * rce_stream = rce.get_stream(i);
  //              auto const identifier = rce_stream->getIdentifier();
  //              uint32_t crateNumber = identifier.getCrate();
  //              uint32_t slotNumber = identifier.getSlot();
  //              uint32_t fiberNumber = identifier.getFiber();
  // 	       if (crateNumber == 3 && slotNumber == 3 && fiberNumber == 2) found302 = true;
  //            } 
  // 	   if (!found302) return;
  // 	}
  //   }

  for (int i = 0; i < rce.size(); ++i)
    {
      auto const * rce_stream = rce.get_stream(i);
      //size_t n_ch = rce_stream->getNChannels();
      size_t n_ticks = rce_stream->getNTicks();
      nticksvec.push_back(n_ticks);
    }

}

bool PDSPTPCRawDecoder::_process_RCE_AUX(
					 const artdaq::Fragment& frag, 
					 RawDigits& raw_digits,
					 RDTimeStamps &timestamps,
					 RDTsAssocs &tsassocs,
					 RDPmkr &rdpm, TSPmkr &tspm,
					 size_t ntickscheck
					 )
{

  if (_rce_enforce_fragment_type_match && (frag.type() != _rce_fragment_type)) 
    {
      MF_LOG_WARNING("_process_RCE_AUX:") << " RCE fragment type " << (int) frag.type() << " doesn't match expected value: " << _rce_fragment_type << " Discarding RCE fragment";
      _DiscardedCorruptData = true;
      return false;
    }
  //MF_LOG_INFO("_Process_RCE_AUX")
  //<< "   SequenceID = " << frag.sequenceID()
  //<< "   fragmentID = " << frag.fragmentID()
  //<< "   fragmentType = " << (unsigned)frag.type()
  //<< "   Timestamp =  " << frag.timestamp();
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;

  dune::RceFragment rce(frag);
  artdaq::Fragment cfragloc(frag);
  size_t cdsize = cfragloc.dataSizeBytes();
  const uint64_t* cdptr = (uint64_t const*) (cfragloc.dataBeginBytes() + 12);  // see dune-raw-data/Overlays/RceFragment.cc
  HeaderFragmentUnpack const cdheader(cdptr);
  //bool isOkay = RceFragmentUnpack::isOkay(cdptr,cdsize+sizeof(cdheader));
  if (cdsize>16) cdsize -= 16;
  bool isOkay = RceFragmentUnpack::isOkay(cdptr,cdsize);
  if (!isOkay)
    {
      MF_LOG_WARNING("_process_RCE_AUX:") << "RCE Fragment isOkay failed: " << cdsize << " Discarding this fragment"; 
      error_counter++;
      _DiscardedCorruptData = true;
      return false; 
    }

  //DataFragmentUnpack df(cdptr);
  //std::cout << "isTPpcNormal: " << df.isTpcNormal() << " isTpcDamaged: " << df.isTpcDamaged() << " isTpcEmpty: " << df.isTpcEmpty() << std::endl;


  uint32_t ch_counter = 0;
  for (int i = 0; i < rce.size(); ++i)
    {
      auto const * rce_stream = rce.get_stream(i);
      size_t n_ch = rce_stream->getNChannels();
      size_t n_ticks = rce_stream->getNTicks();
      if (n_ticks == 0) continue;  // on David Adams's request
      auto const identifier = rce_stream->getIdentifier();
      uint32_t crateNumber = identifier.getCrate();
      uint32_t slotNumber = identifier.getSlot();
      uint32_t fiberNumber = identifier.getFiber();

      //std::cout << "Processing an RCE Stream: " << crateNumber << " " << slotNumber << " " << fiberNumber << " " << n_ticks << " " << n_ch << std::endl;

      size_t adsiz = _apas_to_decode.size();
      bool apafound = true;  // default if no vector to test against
      if (adsiz)
	{
	  apafound = false;
	  for (unsigned int j=0; j<adsiz; ++j)
	    {
	      if (  ((int) crateNumber == _apas_to_decode[j])         || 
		    (_apas_to_decode[j] < 0)                          ||
		    ( (crateNumber == 0 || crateNumber > 6) && _apas_to_decode[j] == 7) ) 
		{
		  apafound = true;
		  break;
		}
	    }
	}
      if (!apafound) 
	{
	  return false;
	}

      // check for bad crate numbers -- default empty list of APAs to check, default check is on, and if crate number isn't
      // one of the DAQ ones.  If we asked for a crate and it's there, even if it's unusual, let it pass.
 

      if ( (crateNumber == 0 || crateNumber > 6) && _rce_drop_frags_with_badc && adsiz == 0)
	{
	  MF_LOG_WARNING("_process_RCE:") << "Bad crate number, discarding fragment on request: " 
					  << (int) crateNumber;
	  return false;
	}

      if (slotNumber > 4 || fiberNumber == 0 || fiberNumber > 4)
	{
	  if (_rce_drop_frags_with_badsf)
	    {
	      MF_LOG_WARNING("_process_RCE:") << "Bad  slot, fiber number, discarding fragment on request: " 
					      << " " << slotNumber << " " << fiberNumber;
              _DiscardedCorruptData = true;
	      return false;
	    }
	  _KeptCorruptData = true;
	}

      if (_print_coldata_convert_count)
	{
	  // from JJ's PdReaderTest.cc
	  using namespace pdd;
	  using namespace pdd::access;
	  bool printed=false;
	  TpcStream const        &stream = rce_stream->getStream ();
	  TpcToc           toc    (stream.getToc    ());
	  TpcPacket        pktRec (stream.getPacket ());
	  TpcPacketBody    pktBdy (pktRec.getRecord ());
	  int   npkts = toc.getNPacketDscs ();
	  for (int ipkt = 0; ipkt < npkts; ++ipkt)
	    {
	      TpcTocPacketDsc pktDsc (toc.getPacketDsc (ipkt));
	      unsigned int      o64 = pktDsc.getOffset64 ();
	      unsigned int  pktType = pktDsc.getType ();
	      unsigned nWibFrames = pktDsc.getNWibFrames ();
	      WibFrame const *wf = pktBdy.getWibFrames (pktType, o64);
	      for (unsigned iwf = 0; iwf < nWibFrames; ++iwf)
		{
		  auto const &colddata = wf->getColdData ();
		  auto cvt0 = colddata[0].getConvertCount ();
		  //auto cvt1 = colddata[1].getConvertCount ();
		  std::cout << "RCE coldata convert count: " << cvt0 << std::endl;
		  printed = true;
		  ++wf;  // in case we were looping over WIB frames, but let's stop at the first
		  break;
		}
	      if (printed) break;
	    }
	}


      if(_make_histograms)
	{
	  //log the participating RCE channels
	  rcechans=rcechans+n_ch;
	}

      // check the number of ticks and allow FEMB302 to have 10% fewer
      size_t ntc10 = ( 0.9 * (float) ntickscheck );
      if (!(
	    n_ticks == ntickscheck ||
	    (
	     (crateNumber == 3 && slotNumber == 3 && fiberNumber == 2) &&
	     n_ticks < ntickscheck && n_ticks > ntc10)
	    )
	  )
	{
	  MF_LOG_WARNING("_process_RCE_AUX:") << "Nticks differs from median or FEMB302 nticks not expected: " << n_ticks << " " 
					      << ntickscheck << " Discarding this fragment";
	  _DiscardedCorruptData = true;
	  return false;
	} 

      if (n_ticks != _full_tick_count)
	{
	  if (_enforce_full_tick_count)
	    {
	      MF_LOG_WARNING("_process_RCE_AUX:") << "Nticks not the required value: " << n_ticks << " " 
						  << _full_tick_count << " Discarding Data";
	      error_counter++;
	      incorrect_ticks++;
	      _discard_data = true;
              _DiscardedCorruptData = true;
	      return false; 
	    }
	  _KeptCorruptData = true;
	}

      if (!_initialized_tick_count_this_event)
	{
	  _initialized_tick_count_this_event = true;
	  _tick_count_this_event = n_ticks;
	}
      else
	{
	  if (n_ticks != _tick_count_this_event)
	    {
	      if (_enforce_same_tick_count)
		{
		  MF_LOG_WARNING("_process_RCE_AUX:") << "Nticks different for two channel streams: " << n_ticks 
						      << " vs " << _tick_count_this_event << " Discarding Data";
		  error_counter++;
		  _discard_data = true;
		  _DiscardedCorruptData = true;
		  return false;
		}
	    }
	  _KeptCorruptData = true;
	}


      //MF_LOG_INFO("_Process_RCE_AUX")
      //<< "RceFragment timestamp: " << rce_stream->getTimeStamp()
      //<< ", NChannels: " << n_ch
      //<< ", NTicks: " << n_ticks;

      // TODO -- speed this up!!  Remove one buffer copy

      size_t buffer_size = n_ch * n_ticks;

      if (buffer_size > _rce_buffer_size_checklimit)
	{
	  if (_rce_check_buffer_size)
	    {
	      MF_LOG_WARNING("_process_RCE_AUX:") << "n_ch*nticks too large: " << n_ch << " * " << n_ticks << " = " << 
		buffer_size << " larger than: " <<  _rce_buffer_size_checklimit << ".  Discarding this fragment";
	      _DiscardedCorruptData = true;
	      return false;
	    }
	  else
	    {
	      _KeptCorruptData = true;
	    }
	}

      if (_buffer.capacity() < buffer_size)
	{
	  //  MF_LOG_INFO("_process_RCE_AUX")
	  //<< "Increase buffer size from " << _buffer.capacity()
	  //<< " to " << buffer_size;

	  _buffer.reserve(buffer_size);
	}

      int16_t* adcs = _buffer.data();
      bool sgmcdretcode = rce_stream->getMultiChannelData(adcs);
      if (!sgmcdretcode)
	{
	  if (_enforce_error_free)
	    {
	      MF_LOG_WARNING("_process_RCE_AUX:") << "getMutliChannelData returns error flag: " 
						  << " c:s:f:ich: " << crateNumber << " " << slotNumber << " " << fiberNumber << " Discarding Data";
	      error_counter++;
              _DiscardedCorruptData = true;
	      return false;
	    }
	  _KeptCorruptData = true;
	}

      //std::cout << "RCE raw decoder trj: " << crateNumber << " " << slotNumber << " " << fiberNumber << std::endl;

      // David Adams's request for channels to start at zero for coldbox test data
      unsigned int crateloc = crateNumber;
      if (crateNumber == 0 || crateNumber > 6) crateloc = _default_crate_if_unexpected;

      raw::RawDigit::ADCvector_t v_adc;
      for (size_t i_ch = 0; i_ch < n_ch; i_ch++)
	{
	  unsigned int offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crateloc, slotNumber, fiberNumber, i_ch, dune::PdspChannelMapService::kRCE);

	  // skip this channel if we are asked to.

	  if (_max_offline_channel >= 0 && _min_offline_channel >= 0 && _max_offline_channel >= _min_offline_channel && 
	      (offlineChannel < (size_t) _min_offline_channel || offlineChannel > (size_t) _max_offline_channel) ) continue;
 
	  v_adc.clear();

	  if (_rce_fix110 && crateNumber == 1 && slotNumber == 0 && fiberNumber == 1 && channelMap->ChipFromOfflineChannel(offlineChannel) == 4 && n_ticks > _rce_fix110_nticks)
	    {
	      for (size_t i_tick = 0; i_tick < n_ticks-_rce_fix110_nticks; i_tick++)
		{
		  v_adc.push_back(adcs[i_tick+_rce_fix110_nticks]);
		}
	      for (size_t i_tick=0; i_tick<_rce_fix110_nticks; ++i_tick)
		{
		  v_adc.push_back(v_adc.back());
		}
	      
	    }
	  else
	    {
	      for (size_t i_tick = 0; i_tick < n_ticks; i_tick++)
		{
		  v_adc.push_back(adcs[i_tick]);
		}
	    }
	  adcs += n_ticks;

	  ch_counter++;

	  if (offlineChannel < _duplicate_channel_checklist_size)
	    {
	      if (_duplicate_channel_checklist[offlineChannel])
		{
		  if(_make_histograms)
		    {
		      duplicate_channels++;
		    }

		  if (_enforce_no_duplicate_channels)
		    {
		      MF_LOG_WARNING("_process_RCE_AUX:") << "Duplicate Channel: " << offlineChannel
							  << " c:s:f:ich: " << crateNumber << " " << slotNumber << " " << fiberNumber << " " << i_ch << " Discarding Data";
		      error_counter++;
		      _discard_data = true;
		      _DiscardedCorruptData = true;
		      return false;
		    }
		  _KeptCorruptData = true;
		}
	      _duplicate_channel_checklist[offlineChannel] = true;
	    }
	  
	  float median=0;
	  float sigma=0;
	  computeMedianSigma(v_adc,median,sigma);

	  /// FEMB 302 IS crate 3, slot 3, fiber 2

	  auto uncompressed_nticks = v_adc.size();  // can be different from n_ticks due to padding of FEMB 302

	  raw::Compress_t cflag=raw::kNone;
	  if (_compress_Huffman)
	    {
	      cflag = raw::kHuffman;
	      raw::Compress(v_adc,cflag);
	    }
	  // here n_ticks is the uncompressed size as required by the constructor
	  raw::RawDigit raw_digit(offlineChannel, uncompressed_nticks, v_adc, cflag);
	  raw_digit.SetPedestal(median,sigma);
	  raw_digits.push_back(raw_digit);  

	  raw::RDTimeStamp rdtimestamp(rce_stream->getTimeStamp(),offlineChannel);
	  timestamps.push_back(rdtimestamp);

	  //associate the raw digit and the timestamp data products
	  auto const rawdigitptr = rdpm(raw_digits.size()-1);
	  auto const rdtimestampptr = tspm(timestamps.size()-1);
	  tsassocs.addSingle(rawdigitptr,rdtimestampptr);            
	}
    }

  return true;
}



bool PDSPTPCRawDecoder::_processFELIX(art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm)
{
  size_t n_felix_frags = 0;
  bool have_data=false;
  bool have_data_nc=false;

  if (_felix_useInputLabels)
    {
      for (size_t ilabel = 0; ilabel < _felix_input_labels.size(); ++ ilabel)
	{
	  if (_felix_input_labels.at(ilabel).find("Container") != std::string::npos)
	    {
              auto cont_frags = evt.getHandle<artdaq::Fragments>(_felix_input_labels.at(ilabel));
	      if (cont_frags)
		{
		  have_data = true;
	          if (! _felixProcContNCFrags(cont_frags, n_felix_frags, true, evt, raw_digits, timestamps, tsassocs, rdpm, tspm))
		    {
		      return false;
		    }
		}
	    }
	  else
	    {
	      auto frags = evt.getHandle<artdaq::Fragments>(_felix_input_labels.at(ilabel));
	      if (frags)
		{
		  have_data_nc = true;
	          if (! _felixProcContNCFrags(frags, n_felix_frags, false, evt, raw_digits, timestamps, tsassocs, rdpm, tspm))
		    {
		      return false;
		    }
		}
	    }
	}
    }
  else  // get all the fragments in the event and look for the ones that say FELIX in them
    {
      //std::vector<art::Handle<artdaq::Fragments> > fraghv;  // fragment handle vector
      //evt.getManyByType(fraghv);
      auto fraghv = evt.getMany<artdaq::Fragments>();

      for (size_t ihandle=0; ihandle<fraghv.size(); ++ihandle)
	{
	  if (fraghv.at(ihandle).provenance()->inputTag().instance().find("FELIX") != std::string::npos)
	    {
	      if (fraghv.at(ihandle).isValid())
		{
	          if (fraghv.at(ihandle).provenance()->inputTag().instance().find("Container") != std::string::npos)
		    {
		      have_data = true;
		      if (! _felixProcContNCFrags(fraghv.at(ihandle), n_felix_frags,true, evt, raw_digits, timestamps, tsassocs, rdpm, tspm) )
			{
			  return false;
			}
		    }
		  else
		    {
		      have_data_nc = true;
		      if (! _felixProcContNCFrags(fraghv.at(ihandle), n_felix_frags, false, evt, raw_digits, timestamps, tsassocs, rdpm, tspm))
			{
			  return false;
			}
		    }
		}

	      // we had swept in all the TPC fragments, possibly for a second time, so remove them
	      // be even more aggressive and remove all cached fragments -- anyone who needs them
	      // should read them in, and getManyByType had swept them all into memory.  Awaiting
	      // a getManyLabelsByType if we go this route

	      else // if (fraghv.at(ihandle).provenance()->inputTag().instance().find("TPC") != std::string::npos) 
		{
		  evt.removeCachedProduct(fraghv.at(ihandle));
		}
	    }
	}
    }

  //MF_LOG_INFO("_processFELIX")
  //<< " Processed " << n_felix_frags
  //<< " FELIX Fragments, making "
  //<< raw_digits.size()
  //<< " RawDigits.";

  // returns true if we want to add to the number of fragments processed.  Separate flag used
  // for data error conditions (_discard_data).

  return have_data || have_data_nc;
}

bool PDSPTPCRawDecoder::_felixProcContNCFrags(art::Handle<artdaq::Fragments> frags, size_t &n_felix_frags, bool is_container, 
					      art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm)
{
  //size of FELIX fragments into histogram
  if(_make_histograms)
    {
      size_t felixbytes = 0;
      for (auto const& frag : *frags)
	{
	  felixbytes = felixbytes + (frag.sizeBytes());
	}
      fFragSizeFELIX->Fill(felixbytes);
    }
    
  for (auto const& frag : *frags)
    {
      //std::cout << "FELIX fragment size bytes: " << frag.sizeBytes() << std::endl; 

      bool process_flag = true;
      if (frag.sizeBytes() < _felix_frag_small_size)
	{
	  if ( _drop_events_with_small_felix_frags )
	    { 
	      MF_LOG_WARNING("_process_FELIX:") << " Small FELIX fragment size: " << frag.sizeBytes() << " Discarding Event on request.";
	      _discard_data = true; 
	      _DiscardedCorruptData = true;
	      evt.removeCachedProduct(frags);
	      return false;
	    }
	  if ( _drop_small_felix_frags )
	    { 
	      MF_LOG_WARNING("_process_FELIX:") << " Small FELIX fragment size: " << frag.sizeBytes() << " Discarding just this fragment on request.";
	      _DiscardedCorruptData = true;
	      process_flag = false;
	    }
	  _KeptCorruptData = true;
	}
      if (process_flag)
	{
	  if (is_container)
	    {
	      artdaq::ContainerFragment cont_frag(frag);
	      for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
		{
		  if (_process_FELIX_AUX(*cont_frag[ii], raw_digits, timestamps, tsassocs, rdpm, tspm)) ++n_felix_frags;
		}
	    }
	  else
	    {
	      if (_process_FELIX_AUX(frag, raw_digits, timestamps,tsassocs, rdpm, tspm)) ++n_felix_frags;
	    }
	}
    }
  evt.removeCachedProduct(frags);
  return true;
}


bool PDSPTPCRawDecoder::_process_FELIX_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits,
					   RDTimeStamps &timestamps,
					   RDTsAssocs &tsassocs,
					   RDPmkr &rdpm, TSPmkr &tspm)
{

  //std::cout 
  //<< "   SequenceID = " << frag.sequenceID()
  //<< "   fragmentID = " << frag.fragmentID()
  //<< "   fragmentType = " << (unsigned)frag.type()
  //<< "   Timestamp =  " << frag.timestamp() << std::endl;

  if (_felix_hex_dump)
    {
      std::ios oldState(nullptr);
      oldState.copyfmt(std::cout);

      std::cout << "FELIX Fragment: all numbers in hex "  << std::hex
		<< "   SequenceID = " << frag.sequenceID()
		<< "   fragmentID = " << frag.fragmentID()
		<< "   fragmentType = " << (unsigned)frag.type()
		<< "   Timestamp =  " << frag.timestamp() << std::endl;
      std::cout << "Offset      Data";
      artdaq::Fragment fragloc(frag);
      unsigned char *dbegin = reinterpret_cast<unsigned char *>(fragloc.dataAddress());
      size_t dsize = fragloc.dataSizeBytes();
      size_t offcounter=0;
      for (size_t bcounter=0; bcounter<dsize;++bcounter)
	{
	  if ( (offcounter % 8) == 0 )
	    {
	      std::cout << std::endl << std::hex << std::setfill('0') << std::setw(8) << offcounter << " ";
	    }
	  std::cout << std::hex << std::setfill('0') << std::setw(2) << (int) *dbegin << " ";
	  dbegin++;
	  offcounter++;
	}
      std::cout << std::endl;
      std::cout.copyfmt(oldState);
    }

  // check against _felix_fragment_type
  if ( _felix_enforce_fragment_type_match && (frag.type() != _felix_fragment_type) )
    {
      _DiscardedCorruptData = true;
      MF_LOG_WARNING("_process_FELIX_AUX:") << " FELIX fragment type " << (int) frag.type() << " doesn't match expected value: " << _felix_fragment_type << " Discarding FELIX fragment";
      return false;
    }

  art::ServiceHandle<dune::PdspChannelMapService> channelMap;

  //Load overlay class.
  dune::FelixFragment felix(frag);

  //Get detector element numbers from the fragment

  uint8_t crate = felix.crate_no(0);
  uint8_t slot = felix.slot_no(0);
  uint8_t fiber = felix.fiber_no(0); // decode this one later 

  size_t adsiz = _apas_to_decode.size();
  bool apafound = true;  // default if no vector to test against
  if (adsiz)
    {
      apafound = false;
      for (unsigned int j=0; j<adsiz; ++j)
	{
	  if (
	      ((int) crate == _apas_to_decode[j])     || 
	      (_apas_to_decode[j] < 0)                ||
	      ( (crate == 0 || crate > 6) && _apas_to_decode[j] == 7) )  // on David Adams's request for coldbox test data 
	    {
	      apafound = true;
	      break;
	    }
	}
    }
  if (!apafound) 
    {
      return false;
    }

  // check for bad crate numbers -- default empty list of APAs to check, default check is on, and if crate number isn't
  // one of the DAQ ones.  If we asked for a crate and it's there, even if it's unusual, let it pass.
 
  if ( (crate == 0 || crate > 6) && _felix_drop_frags_with_badc && adsiz == 0)
    {
      MF_LOG_WARNING("_process_FELIX:") << "Bad crate number, discarding fragment on request: " 
					<< (int) crate;
      return false;
    }

  if ( slot > 4) 
    {
      if (_felix_drop_frags_with_badsf)  // we'll check the fiber later
	{
	  _DiscardedCorruptData = true;
	  MF_LOG_WARNING("_process_FELIX_AUX:") << "Invalid slot:  s=" << (int) slot << " discarding FELIX data.";
	  return false;
	}
      _KeptCorruptData = true;
    }

  if (_print_coldata_convert_count)
    {
      uint16_t first_coldata_convert_count = felix.coldata_convert_count(0,0);
      std::cout << "FELIX Coldata convert count: " << (int) first_coldata_convert_count << std::endl;
    }

  //std::cout << "FELIX raw decoder trj: " << (int) crate << " " << (int) slot << " " << (int) fiber << std::endl;

  const unsigned n_frames = felix.total_frames(); // One frame contains 25 felix (20 ns-long) ticks.  A "frame" is an offline tick
  //std::cout<<" Nframes = "<<n_frames<<std::endl;
  //_h_nframes->Fill(n_frames);
  if (n_frames == 0) return true;
  const unsigned n_channels = dune::FelixFrame::num_ch_per_frame;// should be 256


  if (n_frames*n_channels > _felix_buffer_size_checklimit)
    {
      if (_felix_check_buffer_size)
	{
	  MF_LOG_WARNING("_process_FELIX_AUX:") << "n_channels*n_frames too large: " << n_channels << " * " << n_frames << " = " << 
	    n_frames*n_channels << " larger than: " <<  _felix_buffer_size_checklimit << ".  Discarding this fragment";
	  _DiscardedCorruptData = true;
	  return false;
	}
      else
	{
	  _KeptCorruptData = true;
	}
    }

  if(_make_histograms)
    {
      felixchans=felixchans+n_channels;
    }

  for (unsigned int iframe=0; iframe<n_frames; ++iframe)
    {
      if ( felix.wib_errors(iframe) != 0)
	{
	  if (_enforce_error_free )
	    {
	      _DiscardedCorruptData = true;
	      MF_LOG_WARNING("_process_FELIX_AUX:") << "WIB Errors on frame: " << iframe << " : " << felix.wib_errors(iframe)
						    << " Discarding Data";
	      error_counter++;
	      // drop just this fragment
	      //_discard_data = true;
	      return true;
	    }
	  _KeptCorruptData = true;
	}
    }

  // check optimization of this -- size not reserved

  raw::RawDigit::ADCvector_t v_adc;
  //v_adc.reserve(n_frames*n_channels);
  // Fill the adc vector.  

  for(unsigned ch = 0; ch < n_channels; ++ch) {

    // handle 256 channels on two fibers -- use the channel map that assumes 128 chans per fiber (=FEMB)
    
    unsigned int fiberloc = 0;
    if (fiber == 1) 
      {
	fiberloc = 1;
      }
    else if (fiber == 2)
      {
	fiberloc = 3;
      }
    else
      {
	MF_LOG_WARNING("_process_FELIX_AUX:") << " Fiber number " << (int) fiber << " is expected to be 1 or 2 -- revisit logic";
	fiberloc = 1;
	error_counter++;
	if (_felix_drop_frags_with_badsf) 
	  {
	    MF_LOG_WARNING("_process_FELIX_AUX:") << " Dropping FELIX Data";
	    return false;
	  }
      }

    unsigned int chloc = ch;
    if (chloc > 127)
      {
	chloc -= 128;
	fiberloc++;
      }
    unsigned int crateloc = crate;  

    // David Adams's request for channels to start at zero for coldbox test data
    if (crateloc == 0 || crateloc > 6) crateloc = _default_crate_if_unexpected;  

    unsigned int offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crateloc, slot, fiberloc, chloc, dune::PdspChannelMapService::kFELIX); 

    // skip this channel if we are asked to.

    if (_max_offline_channel >= 0 && _min_offline_channel >= 0 && _max_offline_channel >= _min_offline_channel && 
	(offlineChannel < (size_t) _min_offline_channel || offlineChannel > (size_t) _max_offline_channel) ) continue;

    v_adc.clear();
    //std::cout<<"crate:slot:fiber = "<<crate<<", "<<slot<<", "<<fiber<<std::endl;
    std::vector<dune::adc_t> waveform( felix.get_ADCs_by_channel(ch) );
    for(unsigned int nframe=0;nframe<waveform.size();nframe++){
      // if(ch==0 && nframe<100) {
      //  if(nframe==0) std::cout<<"Print the first 100 ADCs of Channel#1"<<std::endl;  
      //  std::cout<<waveform.at(nframe)<<"  ";
      //  if(nframe==99) std::cout<<std::endl;
      // }
      v_adc.push_back(waveform.at(nframe));  
    }

    if ( v_adc.size() != _full_tick_count)
      {
	if (_enforce_full_tick_count)
	  {
	    MF_LOG_WARNING("_process_FELIX_AUX:") << "Nticks not the required value: " << v_adc.size() << " " 
						  << _full_tick_count << " Discarding Data";
	    error_counter++;
	    incorrect_ticks++;
	    _discard_data = true;
	    _DiscardedCorruptData = true;
	    return true; 
	  }
	_KeptCorruptData = true;
      }

    if (!_initialized_tick_count_this_event)
      {
	_initialized_tick_count_this_event = true;
	_tick_count_this_event = v_adc.size();
      }
    else
      {
	if (_enforce_same_tick_count)
	  {
	    if (v_adc.size() != _tick_count_this_event)
	      {
		MF_LOG_WARNING("_process_FELIX_AUX:") << "Nticks different for two channel streams: " << v_adc.size() 
						      << " vs " << _tick_count_this_event << " Discarding Data";
		error_counter++;
		_discard_data = true;
		_DiscardedCorruptData = true;
		return true;
	      }
	    _KeptCorruptData = true;
	  }
      }

    if (offlineChannel < _duplicate_channel_checklist_size)
      {
	if (_duplicate_channel_checklist[offlineChannel])
	  {
	    if(_make_histograms)
	      {
		duplicate_channels++;
	      }
	    if (_enforce_no_duplicate_channels)
	      {
		MF_LOG_WARNING("_process_FELIX_AUX:") << "Duplicate Channel: " << offlineChannel
						      << " c:s:f:ich: " << (int) crate << " " << (int) slot << " " << (int) fiber << " " << (int) ch << " Discarding Data";
		error_counter++;
		_discard_data = true;
		_DiscardedCorruptData = true;
		return true;
	      }
	    _KeptCorruptData = true;	    
	  }
	_duplicate_channel_checklist[offlineChannel] = true;
      }

    float median=0;
    float sigma=0;
    computeMedianSigma(v_adc,median,sigma);

    auto n_ticks = v_adc.size();
    raw::Compress_t cflag=raw::kNone;
    if (_compress_Huffman)
      {
	cflag = raw::kHuffman;
	raw::Compress(v_adc,cflag);
      }
    // here n_ticks is the uncompressed size as required by the constructor
    raw::RawDigit raw_digit(offlineChannel, n_ticks, v_adc, cflag);
    raw_digit.SetPedestal(median,sigma);
    raw_digits.push_back(raw_digit);

    raw::RDTimeStamp rdtimestamp(felix.timestamp(),offlineChannel);
    timestamps.push_back(rdtimestamp);

    //associate the raw digit and the timestamp data products
    auto const rawdigitptr = rdpm(raw_digits.size()-1);
    auto const rdtimestampptr = tspm(timestamps.size()-1);
    tsassocs.addSingle(rawdigitptr,rdtimestampptr);
  }
  return true;
}


// compute median and sigma.  Sigma is half the distance between the upper and lower bounds of the
// 68% region where 34% is above the median and 34% is below ("centered" on the median).

void PDSPTPCRawDecoder::computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma)
{
  size_t asiz = v_adc.size();
  int imed=0;
  if (asiz == 0)
    {
      median = 0;
      sigma = 0;
    }
  else
    {
      // this is actually faster than the code below by about one second per event.
      // the RMS includes tails from bad samples and signals and may not be the best RMS calc.

      imed = TMath::Median(asiz,v_adc.data()) + 0.01;  // add an offset to make sure the floor gets the right integer
      median = imed;
      sigma = TMath::RMS(asiz,v_adc.data());

      // add in a correction suggested by David Adams, May 6, 2019

      size_t s1 = 0;
      size_t sm = 0;
      for (size_t i=0; i<asiz; ++i)
	{
	  if (v_adc[i] < imed) s1++;
	  if (v_adc[i] == imed) sm++;
	}
      if (sm > 0)
	{
	  float mcorr = (-0.5 + (0.5*(float) asiz - (float) s1)/ ((float) sm) );
	  //if (std::abs(mcorr)>1.0) std::cout << "mcorr: " << mcorr << std::endl;
	  median += mcorr;
	}
    }



  // never do this, but keep the code around in case we want it later

  // if (asiz > 100000000)
  //   {
  //     size_t mednum = asiz/2;
  //     size_t m1snum = mednum - ( (float) asiz )*0.34;
  //     size_t p1snum = mednum + ( (float) asiz )*0.34;

  //     std::map<size_t,size_t> adcmap;
  //     for (auto const adc : v_adc)
  // 	{
  // 	  auto mapiter = adcmap.find(adc);
  // 	  if (mapiter != adcmap.end())
  // 	    {
  // 	      mapiter->second ++;
  // 	    }
  // 	  else
  // 	    {
  // 	      adcmap[adc] = 1;
  // 	    }
  // 	}

  //     // find quantiles, -1 sigma, median, plus 1 sigma
  //     size_t sum = 0;
  //     size_t m1s = 0;
  //     size_t p1s = 0;
  //     size_t m = 0;
  //     for (auto const &mv : adcmap)
  // 	{
  // 	  sum += mv.second;
  // 	  if (m1s == 0 && sum >= m1snum) 
  // 	    {
  // 	      m1s = mv.first;
  // 	    }
  // 	  if (m == 0 && sum >= mednum)
  // 	    {
  // 	      m = mv.first;
  // 	    }
  // 	  if (p1s == 0 && sum >= p1snum)
  // 	    {
  // 	      p1s = mv.first;
  // 	      break;
  // 	    }
  // 	}
  //     median = (float) m;
  //     sigma = ((float) (p1s - m1s))/2.0;
  //   }
}

DEFINE_ART_MODULE(PDSPTPCRawDecoder)
