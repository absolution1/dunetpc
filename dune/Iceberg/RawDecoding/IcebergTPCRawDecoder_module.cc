////////////////////////////////////////////////////////////////////////
// Class:       IcebergTPCRawDecoder
// Plugin Type: producer (art v2_10_03)
// File:        IcebergTPCRawDecoder_module.cc
//
// Generated at Fri Mar  2 15:36:20 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_02_00.  Original code from Jingbo Wang for ProtoDUNE-SP
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include <memory>
#include <cmath>

// ROOT includes
#include "TH1.h"
#include "TStyle.h"
#include "TMath.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/RceFragment.hh"
#include "dune-raw-data/Overlays/FelixFragment.hh"
#include "dune-raw-data/Overlays/Frame14Fragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "dune-raw-data/Overlays/FragmentType.hh"
#include "dune-raw-data/Services/ChannelMap/IcebergChannelMapService.h"
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

class IcebergTPCRawDecoder : public art::EDProducer {

public:
  explicit IcebergTPCRawDecoder(fhicl::ParameterSet const & p);
  IcebergTPCRawDecoder(IcebergTPCRawDecoder const &) = delete;
  IcebergTPCRawDecoder(IcebergTPCRawDecoder &&) = delete;
  IcebergTPCRawDecoder & operator = (IcebergTPCRawDecoder const &) = delete;
  IcebergTPCRawDecoder & operator = (IcebergTPCRawDecoder &&) = delete;
  void produce(art::Event & e) override;

private:
  typedef std::vector<raw::RawDigit> RawDigits;
  typedef std::vector<raw::RDTimeStamp> RDTimeStamps;
  typedef art::Assns<raw::RawDigit,raw::RDTimeStamp> RDTsAssocs;
  typedef art::PtrMaker<raw::RawDigit> RDPmkr;
  typedef art::PtrMaker<raw::RDTimeStamp> TSPmkr;
  typedef std::vector<raw::RDStatus> RDStatuses;

  // configuration parameters

  std::string   _rce_input_label; 
  std::string   _rce_input_container_instance;
  std::string   _rce_input_noncontainer_instance;
  std::string   _felix_input_label; 
  std::string   _felix_input_container_instance;
  std::string   _felix_input_noncontainer_instance;
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
  bool          _rce_drop_frags_with_badcsf;
  bool          _rce_hex_dump;
  bool          _rce_save_frags_to_files;
  bool          _rce_check_buffer_size;
  size_t        _rce_buffer_size_checklimit;

  bool          _felix_hex_dump;
  bool          _felix_drop_frags_with_badcsf;
  bool          _felix_enforce_exact_crate_number;
  int           _felix_crate_number_to_check;
  bool          _drop_events_with_small_felix_frags;
  bool          _drop_small_felix_frags;
  size_t        _felix_frag_small_size;
  bool          _felix_check_buffer_size;
  size_t        _felix_buffer_size_checklimit;

  bool          _compress_Huffman;
  bool          _print_coldata_convert_count;

  //declare histogram data memebers
  bool  _make_histograms;
  unsigned int  duplicate_channels;
  unsigned int  error_counter;
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
  TH1I * fDeltaTimestamp;

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
  bool _processFELIX(art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm);
  bool _process_RCE_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm, uint32_t runNumber);
  bool _process_FELIX_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm, uint32_t runNumber);

  void computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma);

  std::vector<int16_t> _buffer;
};


IcebergTPCRawDecoder::IcebergTPCRawDecoder(fhicl::ParameterSet const & p)
  : EDProducer(p)
{
  std::vector<int> emptyivec;
  _rce_input_label = p.get<std::string>("RCERawDataLabel","daq");
  _rce_input_container_instance = p.get<std::string>("RCERawDataContainerInstance","ContainerTPC");
  _rce_input_noncontainer_instance = p.get<std::string>("RCERawDataNonContainerInstance","TPC");
  _rce_fragment_type = p.get<int>("RCEFragmentType",2);
  _drop_events_with_small_rce_frags = p.get<bool>("RCEDropEventsWithSmallFrags",false);
  _drop_small_rce_frags = p.get<bool>("RCEDropSmallFrags",true);
  _rce_frag_small_size = p.get<unsigned int>("RCESmallFragSize",10000);
  _rce_drop_frags_with_badcsf = p.get<bool>("RCEDropFragsWithBadCSF",true);
  _rce_hex_dump = p.get<bool>("RCEHexDump",false);  
  _rce_save_frags_to_files = p.get<bool>("RCESaveFragsToFiles",false);  
  _rce_check_buffer_size = p.get<bool>("RCECheckBufferSize",true);
  _rce_buffer_size_checklimit = p.get<unsigned int>("RCEBufferSizeCheckLimit",10000000);

  _felix_input_label = p.get<std::string>("FELIXRawDataLabel");
  _felix_input_container_instance = p.get<std::string>("FELIXRawDataContainerInstance","ContainerFELIX");
  _felix_input_noncontainer_instance = p.get<std::string>("FELIXRawDataNonContainerInstance","FELIX");
  _felix_fragment_type = p.get<int>("FELIXFragmentType",8);
  _felix_drop_frags_with_badcsf = p.get<bool>("FELIXDropFragsWithBadCSF",true);
  _felix_hex_dump = p.get<bool>("FELIXHexDump",false);  
  _felix_enforce_exact_crate_number = p.get<bool>("FELIXEnforceExactCrateNumber",false);  
  _felix_crate_number_to_check = p.get<int>("FELIXCrateNumberToCheck",6);  
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

      fDeltaTimestamp = tFileService->make<TH1I>("fDeltaTimestamp", "Delta Timestamp Between Frames", 101, -0.5, 100.5);

    }

}

void IcebergTPCRawDecoder::produce(art::Event &e)
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
      MF_LOG_WARNING("IcebergTPCRawDecoder:") << "Wrong Total number of Channels " << raw_digits.size()  
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

bool IcebergTPCRawDecoder::_processRCE(art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm)
{
  size_t n_rce_frags = 0;
  art::InputTag itag1(_rce_input_label, _rce_input_container_instance);
  auto cont_frags = evt.getHandle<artdaq::Fragments>(itag1);

  bool have_data=false;
  bool have_data_nc=false;

  uint32_t runNumber = evt.run();

  if (cont_frags)
    {
      have_data = true;

      //size of RCE fragments into histogram
      if(_make_histograms)
        {
          size_t rcebytes = 0;
          for (auto const& cont : *cont_frags)
            {
              rcebytes = rcebytes + (cont.sizeBytes());
            }
          fFragSizeRCE->Fill(rcebytes);
        }
    
      for (auto const& cont : *cont_frags)
        {
          //std::cout << "RCE container fragment size bytes: " << cont.sizeBytes() << std::endl; 
          bool process_flag = true;
          if (cont.sizeBytes() < _rce_frag_small_size)
            {
              if ( _drop_events_with_small_rce_frags )
                { 
                  MF_LOG_WARNING("_process_RCE:") << " Small RCE fragment size: " << cont.sizeBytes() << " Discarding Event on request.";
                  _discard_data = true; 
                  _DiscardedCorruptData = true;
                  evt.removeCachedProduct(cont_frags);
                  return false;
                }
              if ( _drop_small_rce_frags )
                { 
                  MF_LOG_WARNING("_process_RCE:") << " Small RCE fragment size: " << cont.sizeBytes() << " Discarding just this fragment on request.";
                  _DiscardedCorruptData = true;
                  process_flag = false;
                }
              _KeptCorruptData = true;
            }
          if (process_flag)
            {
              artdaq::ContainerFragment cont_frag(cont);
              for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
                {
                  if (_process_RCE_AUX(*cont_frag[ii], raw_digits, timestamps, tsassocs, rdpm, tspm, runNumber)) ++n_rce_frags;
                }
            }
        }
      evt.removeCachedProduct(cont_frags);
    }

  //noncontainer frags

  art::InputTag itag2(_rce_input_label, _rce_input_noncontainer_instance);
  auto frags = evt.getHandle<artdaq::Fragments>(itag2);

  if (frags)
    {
      have_data_nc = true;

      //size of RCE fragments into histogram
      if(_make_histograms)
        {
          size_t rcebytes = 0;
          for (auto const& frag: *frags)
            {
              rcebytes = rcebytes + (frag.sizeBytes());
            }
          fFragSizeRCE->Fill(rcebytes);
        }

      for(auto const& frag: *frags)
        {
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
              if (_process_RCE_AUX(frag, raw_digits, timestamps,tsassocs, rdpm, tspm, runNumber)) ++n_rce_frags;
            }
        }
      evt.removeCachedProduct(frags);
    }

  //MF_LOG_INFO("_processRCE")
  //<< " Processed " << n_rce_frags
  //<< " RCE Fragments, making "
  //<< raw_digits.size()
  //<< " RawDigits.";
  return have_data || have_data_nc;
}

// returns true if we want to add to the number of fragments processed.  Separate flag used
// for data error conditions (_discard_data).

bool IcebergTPCRawDecoder::_process_RCE_AUX(
                                            const artdaq::Fragment& frag, 
                                            RawDigits& raw_digits,
                                            RDTimeStamps &timestamps,
                                            RDTsAssocs &tsassocs,
                                            RDPmkr &rdpm, 
                                            TSPmkr &tspm,
                                            uint32_t runNumber
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

  if(frag.type() != _rce_fragment_type) 
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
  art::ServiceHandle<dune::IcebergChannelMapService> channelMap;
  
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
  DataFragmentUnpack df(cdptr);
  //std::cout << "isTPpcNormal: " << df.isTpcNormal() << " isTpcDamaged: " << df.isTpcDamaged() << " isTpcEmpty: " << df.isTpcEmpty() << std::endl;

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

  uint32_t ch_counter = 0;
  for (int i = 0; i < rce.size(); ++i)
    {
      auto const * rce_stream = rce.get_stream(i);
      size_t n_ch = rce_stream->getNChannels();
      size_t n_ticks = rce_stream->getNTicks();
      auto const identifier = rce_stream->getIdentifier();
      uint32_t crateNumber = identifier.getCrate();
      uint32_t slotNumber = identifier.getSlot();
      uint32_t fiberNumber = identifier.getFiber();

      //std::cout << "Processing an RCE Stream: " << crateNumber << " " << slotNumber << " " << fiberNumber << " " << n_ticks << " " << n_ch << std::endl;

      if (crateNumber == 0 || crateNumber > 6 || slotNumber > 4 || fiberNumber == 0 || fiberNumber > 4)
        {
          if (_rce_drop_frags_with_badcsf)
            {
              MF_LOG_WARNING("_process_RCE:") << "Bad crate, slot, fiber number, discarding fragment on request: " 
                                              << crateNumber << " " << slotNumber << " " << fiberNumber;
              _DiscardedCorruptData = true;
              return false;
            }
          _KeptCorruptData = true;
        }

      // inverted ordering on back side, Run 2c (=Run 3)
      // note Shekhar's FEMB number is fiber-1, and WIB is slot+1

      if (runNumber > 2572)
        {
          auto oldfiber = fiberNumber;
          auto oldslot = slotNumber;

          if (oldslot == 0 && oldfiber == 4)
            {
              slotNumber = 1;
              fiberNumber = 3;
            }
          if (oldslot == 1 && oldfiber == 4)
            {
              slotNumber = 0;
              fiberNumber = 3;
            }
          if (oldslot == 1 && oldfiber == 3)
            {
              slotNumber = 0;
              fiberNumber = 4;
            }
          if (oldslot == 0 && oldfiber == 3)
            {
              slotNumber = 1;
              fiberNumber = 4;
            }
        }

      // two cable swaps on June 20, 2019, and go back to the original on Jan 22, 2020

      if (runNumber > 1530 && runNumber < 2572)
        {
          auto oldfiber = fiberNumber;
          auto oldslot = slotNumber;

          // second swap, June 21, 2019 -- see Slack

          if (oldslot == 2 && oldfiber == 1)
            {
              slotNumber = 2;
              fiberNumber = 3;
            }
          if (oldslot == 1 && oldfiber == 1)
            {
              slotNumber = 1;
              fiberNumber = 3;
            }
          if (oldslot == 2 && oldfiber == 3)
            {
              slotNumber = 2;
              fiberNumber = 1;
            }
          if (oldslot == 1 && oldfiber == 3)
            {
              slotNumber = 1;
              fiberNumber = 1;
            }

          oldfiber = fiberNumber;
          oldslot = slotNumber;

          if (oldslot == 0 && oldfiber == 4)
            {
              slotNumber = 1;
              fiberNumber = 3;
            }
          if (oldslot == 1 && oldfiber == 4)
            {
              slotNumber = 0;
              fiberNumber = 3;
            }
          if (oldslot == 0 && oldfiber == 3)
            {
              slotNumber = 1;
              fiberNumber = 4;
            }
          if (oldslot == 1 && oldfiber == 3)
            {
              slotNumber = 0;
              fiberNumber = 4;
            }
        }

      // skip the fake TPC data

      if ( slotNumber == 1 && fiberNumber == 1 ) 
        {
          continue;
        }

      if ( slotNumber == 2 && fiberNumber == 1 )
        {
          continue;
        }

      //std::cout << "  After cable swap : WIB: " << slotNumber+1 << " FEMB: " << fiberNumber-1 << std::endl;

      if (_print_coldata_convert_count)
        {
          std::cout << "Printing coldata convert counts for slot: " << slotNumber << " fiber: " << fiberNumber << std::endl;
          // from JJ's PdReaderTest.cc
          using namespace pdd;
          using namespace pdd::access;
          //bool printed=false;
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
                  int diff = (int) cvt0 - (int) iwf; 
                  std::cout << "Packet: " << ipkt << " WIB frame: " << iwf << " RCE coldata convert count: " << cvt0 << " Difference: " << diff << std::endl;
                  //printed = true;
                  ++wf;  // in case we were looping over WIB frames, but let's stop at the first
                  //break;
                }
              //if (printed) break;
            }
        }


      if(_make_histograms)
        {
          //log the participating RCE channels
          rcechans=rcechans+n_ch;
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

      //std::cout << "RCE raw decoder trj -- adjusted slot and fibers after run 1332: " << crateNumber << " " << slotNumber << " " << fiberNumber << std::endl;

      
      raw::RawDigit::ADCvector_t v_adc;
      for (size_t i_ch = 0; i_ch < n_ch; i_ch++)
        {
          // hardcode crate number 1 so we don't get warning messages
          unsigned int offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(1, slotNumber, fiberNumber, i_ch, dune::IcebergChannelMapService::kRCE);

          v_adc.clear();

          for (size_t i_tick = 0; i_tick < n_ticks; i_tick++)
            {
              v_adc.push_back(adcs[i_tick]);
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


bool IcebergTPCRawDecoder::_processFELIX(art::Event &evt, RawDigits& raw_digits, RDTimeStamps &timestamps, RDTsAssocs &tsassocs, RDPmkr &rdpm, TSPmkr &tspm)
{

  // TODO Use MF_LOG_DEBUG
  //MF_LOG_INFO("_processFELIX") << "-------------------- FELIX RawDecoder -------------------";

  unsigned int n_felix_frags = 0;  

  art::InputTag itag3(_felix_input_label, _felix_input_container_instance);
  auto cont_frags = evt.getHandle<artdaq::Fragments>(itag3);

  bool have_data = false;
  bool have_data_nc = false;

  uint32_t runNumber = evt.run();

  if (cont_frags)
    {
      have_data = true;

      //size of felix fragments into histogram
      if(_make_histograms)
        {
          size_t felixbytes = 0;
          for (auto const& cont : *cont_frags)
            {
              felixbytes = felixbytes + (cont.sizeBytes());
            }
          fFragSizeFELIX->Fill(felixbytes);
        }
    
      for (auto const& cont : *cont_frags)
        {
          bool process_flag = true;
          if (cont.sizeBytes() < _felix_frag_small_size)
            {
              if ( _drop_events_with_small_felix_frags )
                { 
                  MF_LOG_WARNING("_process_FELIX:") << " Small FELIX fragment size: " << cont.sizeBytes() << " Discarding Event on request.";
                  _discard_data = true; 
                  _DiscardedCorruptData = true;
                  evt.removeCachedProduct(cont_frags);
                  return false;
                }
              if ( _drop_small_felix_frags )
                { 
                  MF_LOG_WARNING("_process_FELIX:") << " Small FELIX fragment size: " << cont.sizeBytes() << " Discarding just this fragment on request.";
                  _DiscardedCorruptData = true;
                  process_flag = false;
                }
              _KeptCorruptData = true;
            }
          if (process_flag)
            {
              artdaq::ContainerFragment cont_frag(cont);
              for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
                {
                  if (_process_FELIX_AUX(*cont_frag[ii], raw_digits, timestamps, tsassocs, rdpm, tspm, runNumber)) ++n_felix_frags;
                }
            }
        }
      evt.removeCachedProduct(cont_frags);
    }

  // noncontainer frags

  art::InputTag itag4(_felix_input_label, _felix_input_noncontainer_instance);
  auto frags = evt.getHandle<artdaq::Fragments>(itag4);

  if (frags)
    {

      if(_make_histograms)
        {
          size_t felixbytes = 0;
          for (auto const& frag: *frags)
            {
              felixbytes = felixbytes + (frag.sizeBytes());
            }
          fFragSizeFELIX->Fill(felixbytes);
        }

      for(auto const& frag: *frags)
        {
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
              if (_process_FELIX_AUX(frag, raw_digits,timestamps, tsassocs, rdpm, tspm, runNumber)) ++n_felix_frags;
            }
        }
      evt.removeCachedProduct(frags);
    }


  //MF_LOG_INFO("_processFELIX")
  //<< " Processed " << n_felix_frags
  //<< " FELIX Fragments, total size of raw digits is now "
  //<< raw_digits.size()
  //<< " RawDigits.";

  return have_data || have_data_nc;
}

bool IcebergTPCRawDecoder::_process_FELIX_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits,
                                              RDTimeStamps &timestamps,
                                              RDTsAssocs &tsassocs,
                                              RDPmkr &rdpm, TSPmkr &tspm,
                                              uint32_t runNumber)
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
  if(frag.type() != _felix_fragment_type) 
    {
      _DiscardedCorruptData = true;
      MF_LOG_WARNING("_process_FELIX_AUX:") << " FELIX fragment type " << (int) frag.type() << " doesn't match expected value: " << _felix_fragment_type << " Discarding FELIX fragment";
      return false;
    }
  //if (frag.fragmentID() == 501) 
  //  {
  //    std::cout << "Temporary hack: discarding fragment ID 501" << std::endl;
  //    return false;
  //  }

  art::ServiceHandle<dune::IcebergChannelMapService> channelMap;

  //Load overlay class.

  bool is14 = _felix_input_noncontainer_instance == "FRAME14" || _felix_input_container_instance == "FRAME14";
  std::unique_ptr<dune::FelixFragment> felixptr;
  std::unique_ptr<dune::Frame14FragmentUnordered> frame14ptr;
  if (is14)
    {
      std::unique_ptr<dune::Frame14FragmentUnordered> ftmp(new dune::Frame14FragmentUnordered(frag));
      frame14ptr = std::move(ftmp);
    }
  else
    {
      std::unique_ptr<dune::FelixFragment> ftmp(new dune::FelixFragment(frag));
      felixptr = std::move(ftmp);
    }

  //Get detector element numbers from the fragment

  uint8_t crate = 0;
  uint8_t slot  = 0;
  uint8_t fiber = 0;

  if (is14)
    {
      crate = frame14ptr->crate_no(0);
      slot  = frame14ptr->slot_no(0);
      fiber = frame14ptr->fiber_no(0); // decode this one later
      int frame_version = frame14ptr->frame_version(0);
      //std::cout << "ICEBERG frame_version, crate, slot, fiber, fragID: " << frame_version << " "
      // << (int) crate << " " << (int) slot << " " << (int) fiber << " " << (int) frag.fragmentID() << std::endl;
      if (frame_version == 0) 
        {
          std::cout << "ICEBERG frame_version = 0; skipping this fragment" << std::endl;
          return false;
        }
      fiber ++;  // read in 0 to 1, go from 1 to 2

      // temporary hacks
      crate = 1;
      //int fragid = (int) frag.fragmentID();
      // if (fragid == 600)
      //        {
      //          slot = 0;
      //          fiber = 1;
      //        }
      // if (fragid == 601)
      //        {
      //          slot = 0;
      //          fiber = 2;
      //        }
      // if (fragid == 700)
      //        {
      //          slot = 2;
      //          fiber = 1;
      //        }
      // if (fragid == 701)
      //        {
      //          slot = 2;
      //          fiber = 2;
      //        }
      // if (fragid == 800)
      //        {
      //          slot = 1;
      //          fiber = 1;
      //        }
      // if (fragid == 801)
      //        {
      //          slot = 1;
      //          fiber = 2;
      //        }

    }
  else
    {
      crate = felixptr->crate_no(0);
      slot  = felixptr->slot_no(0);
      fiber = felixptr->fiber_no(0); // decode this one later 
    }

  if (crate == 0 || crate > 6 || slot > 4) 
    {
      if (_felix_drop_frags_with_badcsf)  // we'll check the fiber later
        {
          _DiscardedCorruptData = true;
          MF_LOG_WARNING("_process_FELIX_AUX:") << "Invalid crate or slot: c=" << (int) crate << " s=" << (int) slot << " discarding FELIX data.";
          return false;
        }
      _KeptCorruptData = true;
    }
  if ( _felix_crate_number_to_check > -1 && (int) crate != _felix_crate_number_to_check )
    {
      if (_felix_enforce_exact_crate_number)
        {
          _DiscardedCorruptData = true;
          MF_LOG_WARNING("_process_FELIX_AUX:") << "Crate c=" << (int) crate << " mismatches required crate: " << _felix_crate_number_to_check << " discarding FELIX data.";
          return false;  
        }
      _KeptCorruptData = true;
    }

  if (_print_coldata_convert_count)
    {
      size_t first_coldata_convert_count = 0;
      if (is14)
        {
          first_coldata_convert_count = felixptr->coldata_convert_count(0,0);
        }
      else
        {
          first_coldata_convert_count = frame14ptr->total_frames();
        }
      std::cout << "FELIX Coldata convert count: " << (int) first_coldata_convert_count << std::endl;
    }

  //std::cout << "FELIX raw decoder crate, slot, fiber: " << (int) crate << " " << (int) slot << " " << (int) fiber << std::endl;
  // One frame contains 25 felix (20 ns-long) ticks.  A "frame" is an offline tick
  unsigned n_frames = 0;
  if (is14)
    {
      n_frames = frame14ptr->total_frames();
    }
  else
    {
      n_frames = felixptr->total_frames(); 
    }

  //std::cout<<" Nframes = "<<n_frames<<std::endl;
  //_h_nframes->Fill(n_frames);
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

  // this test does not yet exist for Frame14

  if (!is14)
    {
      for (unsigned int iframe=0; iframe<n_frames; ++iframe)
        {
          if ( felixptr->wib_errors(iframe) != 0)
            {
              if (_enforce_error_free )
                {
                  _DiscardedCorruptData = true;
                  MF_LOG_WARNING("_process_FELIX_AUX:") << "WIB Errors on frame: " << iframe << " : " << felixptr->wib_errors(iframe)
                                                        << " Discarding Data";
                  error_counter++;
                  // drop just this fragment
                  //_discard_data = true;
                  return true;
                }
              _KeptCorruptData = true;
            }
        }
    }

  // check optimization of this -- size not reserved

  raw::RawDigit::ADCvector_t v_adc;
  //v_adc.reserve(n_frames*n_channels);
  // Fill the adc vector.  

  for(unsigned ch = 0; ch < n_channels; ++ch) {
    v_adc.clear();
    //std::cout<<"crate:slot:fiber = "<<crate<<", "<<slot<<", "<<fiber<<std::endl;
    std::vector<dune::adc_t> waveform( is14 ? 
                                       frame14ptr->get_ADCs_by_channel(ch) : 
                                       felixptr->get_ADCs_by_channel(ch) );
    for(unsigned int nframe=0;nframe<waveform.size();nframe++){
      // if(ch==0 && nframe<100) {
      //  if(nframe==0) std::cout<<"Print the first 100 ADCs of Channel#1"<<std::endl;  
      //  std::cout<<waveform.at(nframe)<<"  ";
      //  if(nframe==99) std::cout<<std::endl;
      // }
      v_adc.push_back(waveform.at(nframe));  
    }

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
        if (_felix_drop_frags_with_badcsf) 
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
    //unsigned int crateloc = crate;  

    // conversion of Run 7 slot and fibers to Run 5

    if (runNumber > 9745)
      {
        auto slotloc3 = slot;
        auto fiberloc3 = fiberloc;

        // conversion map
        // run 7 slot, run 7 fiber, run 5 slot, run 5 fiber
        unsigned int sfmap[10][4] = 
          {
            {4, 1,  0, 4},
	    {3, 4,  1, 4},
            {3, 3,  2, 4},
	    {3, 2,  0, 3},
	    {3, 1,  1, 3},
	    {0, 1,  0, 1},
	    {0, 2,  2, 2},
	    {0, 3,  1, 2},
	    {0, 4,  0, 2},
	    {1, 1,  2, 3}
          };
        bool found = false;
        for (size_t imap = 0; imap<10; ++imap)
          {
            if (slot == sfmap[imap][0] && fiberloc == sfmap[imap][1])
              {
                slotloc3 = sfmap[imap][2];
                fiberloc3 = sfmap[imap][3];
                found = true;
                break;
              }
            if (!found)
              {
                std::cout << "Slot, fiber not understood in mapping from Run 7 to Run 5: " << (int) slot << " " << (int) fiberloc << std::endl;
                slotloc3 = 0;
                fiberloc3 = 4;
              }
          }
        slot = slotloc3;
        fiberloc = fiberloc3;
      }

    // inverted ordering on back side, Run 2c (=Run 3)
    // note Shekhar's FEMB number is fiber-1, and WIB is slot+1
    // This goes up to run 5

    auto slotloc2 = slot;
    auto fiberloc2 = fiberloc;

    if (runNumber > 2572)
      {

        if (slot == 0 && fiberloc == 4)
          {
            slotloc2 = 1;
            fiberloc2 = 3;
          }
        if (slot == 1 && fiberloc == 4)
          {
            slotloc2 = 0;
            fiberloc2 = 3;
          }
        if (slot == 1 && fiberloc == 3)
          {
            slotloc2 = 0;
            fiberloc2 = 4;
          }
        if (slot == 0 && fiberloc == 3)
          {
            slotloc2 = 1;
            fiberloc2 = 4;
          }
      }

    // skip the fake TPC data

    if ( slotloc2 == 1 && fiberloc2 == 1 ) 
      {
        continue;
      }

    if ( slotloc2 == 2 && fiberloc2 == 1 )
      {
        continue;
      }

    // for iceberg, hardcode the crate number to suppress warnings
    unsigned int offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(1, slotloc2, fiberloc2, chloc, dune::IcebergChannelMapService::kFELIX); 

    //std::cout << "Calling channel map: " << (int) slotloc2 << " " << fiberloc2 << " " << chloc << " " << offlineChannel << std::endl;

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

    raw::RDTimeStamp rdtimestamp( is14 ? frame14ptr->timestamp() : felixptr->timestamp(),offlineChannel);
    timestamps.push_back(rdtimestamp);

    if (_make_histograms)
      {
        uint64_t last_timestamp = (is14 ? frame14ptr->timestamp(0) : felixptr->timestamp(0));
        for (size_t itick=1; itick<n_ticks; ++itick)
          {
            uint64_t timestamp = (is14 ? frame14ptr->timestamp(itick) : felixptr->timestamp(itick));
            uint64_t tdiff = 0;
            if (timestamp > last_timestamp)
              {
                tdiff = timestamp - last_timestamp;
              }
            else
              {
                tdiff = last_timestamp - timestamp;
              }
            fDeltaTimestamp->Fill(tdiff);
            last_timestamp = timestamp;
          }
      }

    //associate the raw digit and the timestamp data products
    auto const rawdigitptr = rdpm(raw_digits.size()-1);
    auto const rdtimestampptr = tspm(timestamps.size()-1);
    tsassocs.addSingle(rawdigitptr,rdtimestampptr);
  }


  return true;
}


// compute median and sigma.  Sigma is half the distance between the upper and lower bounds of the
// 68% region where 34% is above the median and 34% is below ("centered" on the median).

void IcebergTPCRawDecoder::computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma)
{
  size_t asiz = v_adc.size();
  if (asiz == 0)
    {
      median = 0;
      sigma = 0;
    }
  else
    {
      // this is actually faster than the code below by about one second per event.
      // the RMS includes tails from bad samples and signals and may not be the best RMS calc.

      median = TMath::Median(asiz,v_adc.data());
      sigma = TMath::RMS(asiz,v_adc.data());
    }

  //  std::cout << "sigma: " << sigma << std::endl;
}

DEFINE_ART_MODULE(IcebergTPCRawDecoder)
