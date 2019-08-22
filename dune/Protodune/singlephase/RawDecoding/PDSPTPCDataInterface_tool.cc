// PDSPTPCDataInterface_tool.cc

#include "PDSPTPCDataInterface.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <set>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"


// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/RceFragment.hh"
#include "dune-raw-data/Overlays/FelixFragment.hh"
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

PDSPTPCDataInterface::PDSPTPCDataInterface(fhicl::ParameterSet const& p)
{
  _input_labels_by_apa[1] = p.get< std::vector<std::string> >("APA1InputLabels");
  _input_labels_by_apa[2] = p.get< std::vector<std::string> >("APA2InputLabels");
  _input_labels_by_apa[3] = p.get< std::vector<std::string> >("APA3InputLabels");
  _input_labels_by_apa[4] = p.get< std::vector<std::string> >("APA4InputLabels");
  _input_labels_by_apa[5] = p.get< std::vector<std::string> >("APA5InputLabels");
  _input_labels_by_apa[6] = p.get< std::vector<std::string> >("APA6InputLabels");
  _input_labels_by_apa[-1] = p.get< std::vector<std::string> >("MISCAPAInputLabels");
  
  _drop_small_rce_frags = p.get<bool>("RCEDropSmallFrags",true);
  _rce_frag_small_size = p.get<unsigned int>("RCESmallFragSize",10000);
  _rce_drop_frags_with_badcsf = p.get<bool>("RCEDropFragsWithBadCSF",true);
  _rce_hex_dump = p.get<bool>("RCEHexDump",false);  
  _rce_save_frags_to_files = p.get<bool>("RCESaveFragsToFiles",false);  
  _rce_check_buffer_size = p.get<bool>("RCECheckBufferSize",true);
  _rce_buffer_size_checklimit = p.get<unsigned int>("RCEBufferSizeCheckLimit",10000000);

  // parameters to steer the FEMB 110 band-aid

  _rce_fix110 = p.get<bool>("RCEFIX110",true);
  _rce_fix110_nticks = p.get<unsigned int>("RCEFIX110NTICKS",18);

  _felix_drop_frags_with_badcsf = p.get<bool>("FELIXDropFragsWithBadCSF",true);
  _felix_hex_dump = p.get<bool>("FELIXHexDump",false);  
  _drop_small_felix_frags = p.get<bool>("FELIXDropSmallFrags",true);
  _felix_frag_small_size = p.get<unsigned int>("FELIXSmallFragSize",10000);
  _felix_check_buffer_size = p.get<bool>("FELIXCheckBufferSize",true);
  _felix_buffer_size_checklimit = p.get<unsigned int>("FELIXBufferSizeCheckLimit",10000000);

  _enforce_same_tick_count = p.get<bool>("EnforceSameTickCount",false);
  _enforce_full_tick_count = p.get<bool>("EnforceFullTickCount",false);
  _full_tick_count = p.get<unsigned int>("FullTickCount",6000);
  _enforce_error_free = p.get<bool>("EnforceErrorFree",false);
  _enforce_no_duplicate_channels = p.get<bool>("EnforceNoDuplicateChannels", true);
}

// wrapper for backward compatibility.  Return data for all APA's represented in the fragments on these labels

int PDSPTPCDataInterface::retrieveData(art::Event &evt, 
				       std::string inputLabel, 
				       std::vector<raw::RawDigit> &raw_digits, 
				       std::vector<raw::RDTimeStamp> &rd_timestamps,
		                       art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, 
				       std::vector<raw::RDStatus> &rdstatuses)
{
  std::vector<int> apalist;
  apalist.push_back(-1);
  int retcode = retrieveDataAPAListWithLabels(evt, inputLabel, raw_digits, rd_timestamps, rd_ts_assocs, rdstatuses, apalist );
  return retcode;
}

// get data for specified APAs.  Loop over labels specified in the fcl configuration looking for the data so the caller doesn't have to
// keep track of all the branch labels an APA's data might be on.

int PDSPTPCDataInterface::retrieveDataForSpecifiedAPAs(art::Event &evt, 
							std::vector<raw::RawDigit> &raw_digits, 
							std::vector<raw::RDTimeStamp> &rd_timestamps,
							art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, 
							std::vector<raw::RDStatus> &rdstatuses, 
							std::vector<int> &apalist)
{
  int totretcode = 0;

  for (size_t i=0; i<apalist.size(); ++i)
    { 
      auto lli = _input_labels_by_apa.find(apalist.at(i));
      if (lli == _input_labels_by_apa.end())
	{
	  MF_LOG_WARNING("PDSPTPCDataInterface:") << " No list of input labels known for APA " << apalist.at(i) << " Returning no data.";
	}
      for (size_t j=0; j<lli->second.size(); ++j)
	{
	  int retcode = retrieveDataAPAListWithLabels(evt, lli->second.at(j), raw_digits, rd_timestamps, rd_ts_assocs, rdstatuses, apalist );
	  if (retcode > totretcode) totretcode = retcode; // take most severe retcode of everything
	}
    }
  return totretcode;
}

// get data for a specific label, but only return those raw digits that correspond to APA's on the list

int PDSPTPCDataInterface::retrieveDataAPAListWithLabels(art::Event &evt, 
							std::string inputLabel, 
							std::vector<raw::RawDigit> &raw_digits, 
							std::vector<raw::RDTimeStamp> &rd_timestamps,
							art::Assns<raw::RawDigit,raw::RDTimeStamp> rd_ts_assocs, 
							std::vector<raw::RDStatus> &rdstatuses, 
							std::vector<int> &apalist)
{

  RDPmkr rdpm(evt);
  TSPmkr tspm(evt);

  _initialized_tick_count_this_event = false;
  _discard_data = false;   // true if we're going to drop the whole event's worth of data
  _DiscardedCorruptData = false;   // can be set to true if we drop some of the event's data
  _KeptCorruptData = false;      // true if we identify a corruption candidate but are skipping the test to drop it

  if (inputLabel.find("TPC") != std::string::npos)
    {
      _processRCE(evt, inputLabel, raw_digits, rd_timestamps, rd_ts_assocs, rdpm, tspm, apalist);
    }
  else if (inputLabel.find("FELIX") != std::string::npos)
    {
      _processFELIX(evt, inputLabel, raw_digits, rd_timestamps, rd_ts_assocs, rdpm, tspm, apalist);
    }
  else
    {
      throw cet::exception("PDSPTPCInterface_tool") << "ununderstood fragment branch label: \"" << inputLabel << "\"";
    }

  bool flagged_duplicate = false;

  if (_enforce_no_duplicate_channels)
    {
      std::set<unsigned int> channels_seen;

      for (const auto& rd : raw_digits)
	{
	  unsigned int ichan = rd.Channel();
	  if (channels_seen.find(ichan) == channels_seen.end())
	    {
	      channels_seen.insert(ichan);
	    }
	  else
	    {
	      MF_LOG_WARNING("PDSPTPCDataInterface:") << " Duplicate channel detected: " << ichan << " Discarding TPC data for this chunk: " << inputLabel;
	      _discard_data = true;
	      raw_digits.clear();
	      rd_timestamps.clear();
	      RDTsAssocs nullassocs;
	      rd_ts_assocs = nullassocs;
	      flagged_duplicate = true;
	      break;
	    }
	}
    }

  unsigned int statword=0;
  if (_DiscardedCorruptData) statword |= 1;
  if (_KeptCorruptData) statword |= 2;
  rdstatuses.emplace_back(_DiscardedCorruptData,_KeptCorruptData,statword);
  if (flagged_duplicate) statword = 4;  // a flag to the caller indicating that the entire event's worth of raw digits is to be dropped
  return statword;
}



bool PDSPTPCDataInterface::_processRCE(art::Event &evt, 
				       std::string inputLabel,  
				       RawDigits& raw_digits, 
				       RDTimeStamps &timestamps, 
				       RDTsAssocs &tsassocs, 
				       RDPmkr &rdpm, 
				       TSPmkr &tspm,
				       std::vector<int> &apalist)
{
  size_t n_rce_frags = 0;
  bool have_data=false;
  bool have_data_nc=false;

  if (inputLabel.find("Container") != std::string::npos)
    {
      art::Handle<artdaq::Fragments> cont_frags;
      evt.getByLabel(inputLabel,cont_frags);  
      if (cont_frags.isValid())
	{
	  have_data = true;
	  if (! _rceProcContNCFrags(cont_frags, n_rce_frags, true, evt, raw_digits, timestamps, tsassocs, rdpm, tspm, apalist))
	    {
	      return false;
	    }
	}
    }
  else
    {
      art::Handle<artdaq::Fragments> frags;
      evt.getByLabel(inputLabel, frags); 

      if (frags.isValid())
	{
	  have_data_nc = true;
	  if (! _rceProcContNCFrags(frags, n_rce_frags, false, evt, raw_digits, timestamps, tsassocs, rdpm, tspm, apalist))
	    {
	      return false;
	    }
	}
    }

  // returns true if we want to add to the number of fragments processed.  Separate flag used
  // for data error conditions (_discard_data).

  return have_data || have_data_nc;
}

bool PDSPTPCDataInterface::_rceProcContNCFrags(art::Handle<artdaq::Fragments> frags, 
					       size_t &n_rce_frags, 
					       bool is_container, 
					       art::Event &evt, 
					       RawDigits& raw_digits, 
					       RDTimeStamps &timestamps, 
					       RDTsAssocs &tsassocs, 
					       RDPmkr &rdpm, 
					       TSPmkr &tspm,
					       std::vector<int> &apalist)
{
    
  for (auto const& frag : *frags)
    {
      //std::cout << "RCE fragment size bytes: " << frag.sizeBytes() << std::endl; 

      bool process_flag = true;
      if (frag.sizeBytes() < _rce_frag_small_size)
	{
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
		  if (_process_RCE_AUX(*cont_frag[ii], raw_digits, timestamps, tsassocs, rdpm, tspm, apalist)) ++n_rce_frags;
		}
	    }
	  else
	    {
	      if (_process_RCE_AUX(frag, raw_digits, timestamps,tsassocs, rdpm, tspm, apalist)) ++n_rce_frags;
	    }
	}
    }
  evt.removeCachedProduct(frags);  // do this always, even if we need to re-read a TBranch
  return true;
}


bool PDSPTPCDataInterface::_process_RCE_AUX(
					    const artdaq::Fragment& frag, 
					    RawDigits& raw_digits,
					    RDTimeStamps &timestamps,
					    RDTsAssocs &tsassocs,
					    RDPmkr &rdpm, TSPmkr &tspm,
					    std::vector<int> &apalist
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

  //MF_LOG_INFO("_Process_RCE_AUX")
  //<< "   SequenceID = " << frag.sequenceID()
  //<< "   fragmentID = " << frag.fragmentID()
  //<< "   fragmentType = " << (unsigned)frag.type()
  //<< "   Timestamp =  " << frag.timestamp();
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;
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
  if (!isOkay)
    {
      MF_LOG_WARNING("_process_RCE_AUX:") << "RCE Fragment isOkay failed: " << cdsize << " Discarding this fragment";
      _DiscardedCorruptData = true;
      return false;
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

      // only take this rce stream if it has data from an APA we are interested in
      bool foundapainlist = false;
      for (size_t ialist=0; ialist < apalist.size(); ++ ialist)
	{
	  if (apalist[ialist] == -1 || apalist[ialist] == (int) crateNumber)
	    {
	      foundapainlist = true;
	      break;
	    }
	}
      if (!foundapainlist) continue;
      
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

      if (n_ticks != _full_tick_count)
	{
	  if (_enforce_full_tick_count)
	    {
	      MF_LOG_WARNING("_process_RCE_AUX:") << "Nticks not the required value: " << n_ticks << " " 
						  << _full_tick_count << " Discarding Data";
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

      std::vector<int16_t> _buffer(buffer_size);

      int16_t* adcs = _buffer.data();
      bool sgmcdretcode = rce_stream->getMultiChannelData(adcs);
      if (!sgmcdretcode)
	{
	  if (_enforce_error_free)
	    {
	      MF_LOG_WARNING("_process_RCE_AUX:") << "getMutliChannelData returns error flag: " 
						  << " c:s:f:ich: " << crateNumber << " " << slotNumber << " " << fiberNumber << " Discarding Data";
              _DiscardedCorruptData = true;
	      return false;
	    }
	  _KeptCorruptData = true;
	}

      //std::cout << "RCE raw decoder trj: " << crateNumber << " " << slotNumber << " " << fiberNumber << std::endl;

      raw::RawDigit::ADCvector_t v_adc;
      for (size_t i_ch = 0; i_ch < n_ch; i_ch++)
	{
	  unsigned int offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crateNumber, slotNumber, fiberNumber, i_ch, dune::PdspChannelMapService::kRCE);

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

	  float median=0;
	  float sigma=0;
	  computeMedianSigma(v_adc,median,sigma);

	  /// FEMB 302 IS crate 3, slot 3, fiber 2

	  auto uncompressed_nticks = v_adc.size();  // can be different from n_ticks due to padding of FEMB 302

	  raw::Compress_t cflag=raw::kNone;
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



bool PDSPTPCDataInterface::_processFELIX(art::Event &evt, 
					 std::string inputLabel, 
					 RawDigits& raw_digits, 
					 RDTimeStamps &timestamps, 
					 RDTsAssocs &tsassocs, 
					 RDPmkr &rdpm, 
					 TSPmkr &tspm,
					 std::vector<int> &apalist)
{
  size_t n_felix_frags = 0;
  bool have_data=false;
  bool have_data_nc=false;

  if (inputLabel.find("Container") != std::string::npos)
    {
      art::Handle<artdaq::Fragments> cont_frags;
      evt.getByLabel(inputLabel,cont_frags);  
      if (cont_frags.isValid())
	{
	  have_data = true;
	  if (! _felixProcContNCFrags(cont_frags, n_felix_frags, true, evt, raw_digits, timestamps, tsassocs, rdpm, tspm, apalist))
	    {
	      return false;
	    }
	}
    }
  else
    {
      art::Handle<artdaq::Fragments> frags;
      evt.getByLabel(inputLabel, frags); 

      if (frags.isValid())
	{
	  have_data_nc = true;
	  if (! _felixProcContNCFrags(frags, n_felix_frags, false, evt, raw_digits, timestamps, tsassocs, rdpm, tspm, apalist))
	    {
	      return false;
	    }
	}
    }

  // returns true if we want to add to the number of fragments processed.  Separate flag used
  // for data error conditions (_discard_data).

  return have_data || have_data_nc;
}

bool PDSPTPCDataInterface::_felixProcContNCFrags(art::Handle<artdaq::Fragments> frags, 
						 size_t &n_felix_frags, 
						 bool is_container, 
						 art::Event &evt, 
						 RawDigits& raw_digits, 
						 RDTimeStamps &timestamps, 
						 RDTsAssocs &tsassocs, 
						 RDPmkr &rdpm, 
						 TSPmkr &tspm,
						 std::vector<int> &apalist)
{
  for (auto const& frag : *frags)
    {
      //std::cout << "FELIX fragment size bytes: " << frag.sizeBytes() << std::endl; 

      bool process_flag = true;
      if (frag.sizeBytes() < _felix_frag_small_size)
	{
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
		  if (_process_FELIX_AUX(*cont_frag[ii], raw_digits, timestamps, tsassocs, rdpm, tspm, apalist)) ++n_felix_frags;
		}
	    }
	  else
	    {
	      if (_process_FELIX_AUX(frag, raw_digits, timestamps,tsassocs, rdpm, tspm, apalist)) ++n_felix_frags;
	    }
	}
    }
  evt.removeCachedProduct(frags);
  return true;
}


bool PDSPTPCDataInterface::_process_FELIX_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits,
					      RDTimeStamps &timestamps,
					      RDTsAssocs &tsassocs,
					      RDPmkr &rdpm, TSPmkr &tspm,
					      std::vector<int> &apalist)
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

  art::ServiceHandle<dune::PdspChannelMapService> channelMap;

  //Load overlay class.
  dune::FelixFragment felix(frag);

  //Get detector element numbers from the fragment

  uint8_t crate = felix.crate_no(0);
  uint8_t slot = felix.slot_no(0);
  uint8_t fiber = felix.fiber_no(0); // decode this one later 

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

  // only take this felix fragment if it has data from an APA we are interested in
  bool foundapainlist = false;
  for (size_t ialist=0; ialist < apalist.size(); ++ ialist)
    {
      if (apalist[ialist] == -1 || apalist[ialist] == (int) crate)
	{
	  foundapainlist = true;
	  break;
	}
    }
  if (!foundapainlist) return true;

  //std::cout << "FELIX raw decoder trj: " << (int) crate << " " << (int) slot << " " << (int) fiber << std::endl;

  const unsigned n_frames = felix.total_frames(); // One frame contains 25 felix (20 ns-long) ticks.  A "frame" is an offline tick
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

  for (unsigned int iframe=0; iframe<n_frames; ++iframe)
    {
      if ( felix.wib_errors(iframe) != 0)
	{
	  if (_enforce_error_free )
	    {
	      _DiscardedCorruptData = true;
	      MF_LOG_WARNING("_process_FELIX_AUX:") << "WIB Errors on frame: " << iframe << " : " << felix.wib_errors(iframe)
						    << " Discarding Data";
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
    v_adc.clear();
    std::vector<dune::adc_t> waveform( felix.get_ADCs_by_channel(ch) );
    for(unsigned int nframe=0;nframe<waveform.size();nframe++){
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
    unsigned int crateloc = crate;  

    unsigned int offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crateloc, slot, fiberloc, chloc, dune::PdspChannelMapService::kFELIX); 

    if ( v_adc.size() != _full_tick_count)
      {
	if (_enforce_full_tick_count)
	  {
	    MF_LOG_WARNING("_process_FELIX_AUX:") << "Nticks not the required value: " << v_adc.size() << " " 
						  << _full_tick_count << " Discarding Data";
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
		_discard_data = true;
		_DiscardedCorruptData = true;
		return true;
	      }
	    _KeptCorruptData = true;
	  }
      }

    float median=0;
    float sigma=0;
    computeMedianSigma(v_adc,median,sigma);

    auto n_ticks = v_adc.size();
    raw::Compress_t cflag=raw::kNone;
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

void PDSPTPCDataInterface::computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma)
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
}
