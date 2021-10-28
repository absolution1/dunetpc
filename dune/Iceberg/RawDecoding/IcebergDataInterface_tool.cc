// IcebergDataInterface_tool.cc

#include "IcebergDataInterface.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <set>

#include "art/Framework/Services/Registry/ServiceHandle.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/RceFragment.hh"
#include "dune-raw-data/Overlays/FelixFragment.hh"
#include "dune-raw-data/Overlays/Frame14Fragment.hh"
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

IcebergDataInterface::IcebergDataInterface(fhicl::ParameterSet const& p)
{
  _input_labels_by_apa[1] = p.get< std::vector<std::string> >("APA1InputLabels");
  _input_labels_by_apa[2] = p.get< std::vector<std::string> >("APA2InputLabels");
  _input_labels_by_apa[3] = p.get< std::vector<std::string> >("APA3InputLabels");
  _input_labels_by_apa[4] = p.get< std::vector<std::string> >("APA4InputLabels");
  _input_labels_by_apa[5] = p.get< std::vector<std::string> >("APA5InputLabels");
  _input_labels_by_apa[6] = p.get< std::vector<std::string> >("APA6InputLabels");
  _input_labels_by_apa[7] = p.get< std::vector<std::string> >("APA7InputLabels");
  _input_labels_by_apa[8] = p.get< std::vector<std::string> >("APA8InputLabels");
  _input_labels_by_apa[-1] = p.get< std::vector<std::string> >("MISCAPAInputLabels");
  
  _default_crate_if_unexpected = p.get<unsigned int>("DefaultCrateIfUnexpected",3);

  _min_offline_channel = p.get<long int>("MinOfflineChannel",-1);
  _max_offline_channel = p.get<long int>("MaxOfflineChannel",-1);

  _drop_small_rce_frags = p.get<bool>("RCEDropSmallFrags",true);
  _rce_frag_small_size = p.get<unsigned int>("RCESmallFragSize",10000);
  _rce_drop_frags_with_badsf = p.get<bool>("RCEDropFragsWithBadSF",true);
  _rce_drop_frags_with_badc = p.get<bool>("RCEDropFragsWithBadC",true);
  _rce_hex_dump = p.get<bool>("RCEHexDump",false);  
  _rce_save_frags_to_files = p.get<bool>("RCESaveFragsToFiles",false);  
  _rce_check_buffer_size = p.get<bool>("RCECheckBufferSize",true);
  _rce_buffer_size_checklimit = p.get<unsigned int>("RCEBufferSizeCheckLimit",10000000);

  _felix_drop_frags_with_badsf = p.get<bool>("FELIXDropFragsWithBadSF",true);
  _felix_drop_frags_with_badc = p.get<bool>("FELIXDropFragsWithBadC",true);
  _felix_hex_dump = p.get<bool>("FELIXHexDump",false);  
  _drop_small_felix_frags = p.get<bool>("FELIXDropSmallFrags",true);
  _felix_frag_small_size = p.get<unsigned int>("FELIXSmallFragSize",10000);
  _felix_check_buffer_size = p.get<bool>("FELIXCheckBufferSize",true);
  _felix_buffer_size_checklimit = p.get<unsigned int>("FELIXBufferSizeCheckLimit",10000000);

  _enforce_same_tick_count = p.get<bool>("EnforceSameTickCount",false);
  _enforce_median_tick_count = p.get<bool>("EnforceMedianTickCount",false);
  _enforce_full_tick_count = p.get<bool>("EnforceFullTickCount",false);
  _full_tick_count = p.get<unsigned int>("FullTickCount",6000);
  _enforce_error_free = p.get<bool>("EnforceErrorFree",false);
  _enforce_no_duplicate_channels = p.get<bool>("EnforceNoDuplicateChannels", true);
}

// wrapper for backward compatibility.  Return data for all APA's represented in the fragments on these labels

int IcebergDataInterface::retrieveData(art::Event &evt, 
                                       std::string inputLabel, 
                                       std::vector<raw::RawDigit> &raw_digits, 
                                       std::vector<raw::RDTimeStamp> &rd_timestamps,
                                       std::vector<raw::RDStatus> &rdstatuses)
{
  std::vector<int> apalist;
  apalist.push_back(-1);
  int retcode = retrieveDataAPAListWithLabels(evt, inputLabel, raw_digits, rd_timestamps, rdstatuses, apalist );
  _collectRDStatus(rdstatuses);
  return retcode;
}

// get data for specified APAs.  Loop over labels specified in the fcl configuration looking for the data so the caller doesn't have to
// keep track of all the branch labels an APA's data might be on.

int IcebergDataInterface::retrieveDataForSpecifiedAPAs(art::Event &evt, 
                                                       std::vector<raw::RawDigit> &raw_digits, 
                                                       std::vector<raw::RDTimeStamp> &rd_timestamps,
                                                       std::vector<raw::RDStatus> &rdstatuses, 
                                                       std::vector<int> &apalist)
{
  int totretcode = 0;

  for (size_t i=0; i<apalist.size(); ++i)
    { 
      auto lli = _input_labels_by_apa.find(apalist.at(i));
      if (lli == _input_labels_by_apa.end())
        {
          MF_LOG_WARNING("IcebergDataInterface:") << " No list of input labels known for APA " << apalist.at(i) << " Returning no data.";
        }
      for (size_t j=0; j<lli->second.size(); ++j)
        {
          int retcode = retrieveDataAPAListWithLabels(evt, lli->second.at(j), raw_digits, rd_timestamps, rdstatuses, apalist );
          if (retcode > totretcode) totretcode = retcode; // take most severe retcode of everything
        }
    }
  _collectRDStatus(rdstatuses);
  return totretcode;
}

// get data for a specific label, but only return those raw digits that correspond to APA's on the list

int IcebergDataInterface::retrieveDataAPAListWithLabels(art::Event &evt, 
                                                        std::string inputLabel, 
                                                        std::vector<raw::RawDigit> &raw_digits, 
                                                        std::vector<raw::RDTimeStamp> &rd_timestamps,
                                                        std::vector<raw::RDStatus> &rdstatuses, 
                                                        std::vector<int> &apalist)
{

  art::ServiceHandle<dune::IcebergChannelMapService> cmap;

  _initialized_tick_count_this_event = false;
  _DiscardedCorruptData = false;   // can be set to true if we drop some of the event's data
  _KeptCorruptData = false;      // true if we identify a corruption candidate but are skipping the test to drop it

  if (inputLabel.find("TPC") != std::string::npos)
    {
      _processRCE(evt, inputLabel, raw_digits, rd_timestamps, apalist);
    }
  else if ( (inputLabel.find("FELIX") != std::string::npos) ||
            (inputLabel.find("FRAME14") != std::string::npos) )
    {
      _processFELIX(evt, inputLabel, raw_digits, rd_timestamps, apalist);
    }
  else
    {
      throw cet::exception("IcebergInterface_tool") << "ununderstood fragment branch label: \"" << inputLabel << "\"";
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
              MF_LOG_WARNING("IcebergDataInterface:") << " Duplicate channel detected: " << ichan << " Discarding TPC data for this chunk: " << inputLabel;
              raw_digits.clear();
              rd_timestamps.clear();
              flagged_duplicate = true;
              break;
            }
        }
    }

  if (_enforce_median_tick_count)
    {

      // find the median tick count and if a channel has a different tick count from median, remove it from the list.  Special dispensation for FEMB302,
      // allowing it to have up to 10% missing ticks.

      std::vector<size_t> ticklist;
      for (size_t i=0; i<raw_digits.size(); ++i)
        {
          ticklist.push_back(raw_digits.at(i).Samples());
        }
      size_t tls = ticklist.size();
      size_t tickmed = 0;
      if (tls != 0)
        {
          tickmed = TMath::Median(tls,ticklist.data());
        }
      size_t tickexample=0;
      std::vector<size_t> dlist;
      for (size_t i=0; i<raw_digits.size(); ++i)
        {
          if (ticklist.at(i) != tickmed)
            {
              unsigned int channel = raw_digits.at(i).Channel();
              unsigned int crate = cmap->APAFromOfflineChannel(channel);
              unsigned int slot = cmap->WIBFromOfflineChannel(channel);
              unsigned int fiber = cmap->FEMBFromOfflineChannel(channel);
              //std::cout << "tick not at median: " << channel << " " << crate << " " << slot << " " << fiber << " " << ticklist.at(i) << " " << tickmed << std::endl;
              if ( (crate == 3) && (slot == 3) && (fiber == 2) && ( ticklist.at(i) > 0.9*tickmed && ticklist.at(i) < tickmed ) ) continue;  // FEMB 302
              dlist.push_back(i);
              tickexample = ticklist.at(i);
            }
        }

      if (dlist.size() != 0)
        {
          //std::cout << "IcebergDataInterface_tool: Discarding data with n_ticks not at the median. Example invalid ticks: " << tickexample << std::endl;
          MF_LOG_WARNING("IcebergDataInterface_tool:") << " Discarding data with n_ticks not at the median. Example invalid ticks: " << tickexample;
          for (size_t i=dlist.size(); i>0; --i)
            {
              raw_digits.erase(raw_digits.begin() + dlist.at(i-1));
              rd_timestamps.erase(rd_timestamps.begin() + dlist.at(i-1));
            } 
        }
    }

  unsigned int statword=0;
  if (_DiscardedCorruptData) statword |= 1;
  if (_KeptCorruptData) statword |= 2;
  rdstatuses.emplace_back(_DiscardedCorruptData,_KeptCorruptData,statword);
  if (flagged_duplicate) statword = 4;  // a flag to the caller indicating that the entire event's worth of raw digits is to be dropped
  _collectRDStatus(rdstatuses);
  return statword;
}


void IcebergDataInterface::_collectRDStatus(std::vector<raw::RDStatus> &rdstatuses)
{
  if (rdstatuses.size() < 2) return; 
  unsigned int statword=0;
  bool dcflag = false;
  bool kcflag = false;
  for (size_t i=0; i<rdstatuses.size(); ++i)
    {
      statword |= rdstatuses.at(i).GetStatWord();
      dcflag |= rdstatuses.at(i).GetCorruptDataDroppedFlag();
      kcflag |= rdstatuses.at(i).GetCorruptDataKeptFlag();
    }
  rdstatuses.clear();
  rdstatuses.emplace_back(dcflag,kcflag,statword);
}


bool IcebergDataInterface::_processRCE(art::Event &evt, 
                                       std::string inputLabel,  
                                       RawDigits& raw_digits, 
                                       RDTimeStamps &timestamps, 
                                       std::vector<int> &apalist)
{
  size_t n_rce_frags = 0;
  bool have_data=false;
  bool have_data_nc=false;

  if (inputLabel.find("Container") != std::string::npos)
    {
      auto cont_frags = evt.getHandle<artdaq::Fragments>(inputLabel);
      if (cont_frags)
        {
          have_data = true;
          if (! _rceProcContNCFrags(cont_frags, n_rce_frags, true, evt, raw_digits, timestamps, apalist))
            {
              return false;
            }
        }
    }
  else
    {
      auto frags = evt.getHandle<artdaq::Fragments>(inputLabel);
      if (frags)
        {
          have_data_nc = true;
          if (! _rceProcContNCFrags(frags, n_rce_frags, false, evt, raw_digits, timestamps, apalist))
            {
              return false;
            }
        }
    }

  // returns true if we want to add to the number of fragments processed. 

  return have_data || have_data_nc;
}

bool IcebergDataInterface::_rceProcContNCFrags(art::Handle<artdaq::Fragments> frags, 
                                               size_t &n_rce_frags, 
                                               bool is_container, 
                                               art::Event &evt, 
                                               RawDigits& raw_digits, 
                                               RDTimeStamps &timestamps, 
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
          else
            {
              _KeptCorruptData = true;
            }
        }
      if (process_flag)
        {
          if (is_container)
            {
              artdaq::ContainerFragment cont_frag(frag);
              for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
                {
                  if (_process_RCE_AUX(evt,*cont_frag[ii], raw_digits, timestamps, apalist)) ++n_rce_frags;
                }
            }
          else
            {
              if (_process_RCE_AUX(evt,frag, raw_digits, timestamps, apalist)) ++n_rce_frags;
            }
        }
    }
  evt.removeCachedProduct(frags);  // do this always, even if we need to re-read a TBranch
  return true;
}


bool IcebergDataInterface::_process_RCE_AUX(
                                            art::Event &evt,
                                            const artdaq::Fragment& frag, 
                                            RawDigits& raw_digits,
                                            RDTimeStamps &timestamps,
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
      _DiscardedCorruptData = true;
      return false;
    }

  //DataFragmentUnpack df(cdptr);
  //std::cout << "isTPpcNormal: " << df.isTpcNormal() << " isTpcDamaged: " << df.isTpcDamaged() << " isTpcEmpty: " << df.isTpcEmpty() << std::endl;

  dune::RceFragment rce(frag);
  if (_rce_save_frags_to_files)
    {
      TString outfilename="rce_";
      outfilename += evt.run();
      outfilename += "_";
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
      if (n_ticks == 0) continue;  // on David Adams's request.
      auto const identifier = rce_stream->getIdentifier();
      uint32_t crateNumber = identifier.getCrate();
      // always use crate #3 for iceberg to fool the ProtoDUNE data prep into thinking it should
      // have the low-numbered channels
      crateNumber = 3;
      uint32_t slotNumber = identifier.getSlot();
      uint32_t fiberNumber = identifier.getFiber();

      // only take this rce stream if it has data from an APA we are interested in
      bool foundapainlist = false;
      for (size_t ialist=0; ialist < apalist.size(); ++ ialist)
        {
          if (  ( (apalist[ialist] == -1) && (!_rce_drop_frags_with_badc || (crateNumber >0 && crateNumber < 7)) )  || 
                (apalist[ialist] == (int) crateNumber) ||
                (apalist[ialist] == 7 && (crateNumber == 0 || crateNumber > 6)) )
            {
              foundapainlist = true;
              break;
            }
        }
      if (!foundapainlist) continue;
      
      //std::cout << "Processing an RCE Stream: " << crateNumber << " " << slotNumber << " " << fiberNumber << " " << n_ticks << " " << n_ch << std::endl;

      //------------------------------------------------------------------------------------------
      // Handle run-dependent re-cabling of the detector
      //------------------------------------------------------------------------------------------

      uint32_t runNumber = evt.run();

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

      //------------------------------------------------------------------------------------------
      // End of handling of run-dependent re-cabling of the detector
      //------------------------------------------------------------------------------------------

      if (slotNumber > 4 || fiberNumber == 0 || fiberNumber > 4)
        {
          if (_rce_drop_frags_with_badsf)
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

      // David Adams's request for channels to start at zero for coldbox test data
      unsigned int crateloc = crateNumber;
      //if (crateNumber == 0 || crateNumber > 6) crateloc = _default_crate_if_unexpected;
      crateloc = 1; // always use crate 1 for Iceberg

      raw::RawDigit::ADCvector_t v_adc;
      for (size_t i_ch = 0; i_ch < n_ch; i_ch++)
        {
          unsigned int offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crateloc, slotNumber, fiberNumber, i_ch, dune::IcebergChannelMapService::kRCE);

          if (_max_offline_channel >= 0 && _min_offline_channel >= 0 && _max_offline_channel >= _min_offline_channel && 
              (offlineChannel < (size_t) _min_offline_channel || offlineChannel > (size_t) _max_offline_channel) ) continue;

          v_adc.clear();

          for (size_t i_tick = 0; i_tick < n_ticks; i_tick++)
            {
              v_adc.push_back(adcs[i_tick]);
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

        } // end loop over channels
    }  // end loop over RCE streams (1 per FEMB)
  // end processing one RCE fragment
  return true;
}



bool IcebergDataInterface::_processFELIX(art::Event &evt, 
                                         std::string inputLabel, 
                                         RawDigits& raw_digits, 
                                         RDTimeStamps &timestamps, 
                                         std::vector<int> &apalist)
{
  size_t n_felix_frags = 0;
  bool have_data=false;
  bool have_data_nc=false;

  if (inputLabel.find("Container") != std::string::npos)
    {
      auto cont_frags = evt.getHandle<artdaq::Fragments>(inputLabel);
      if (cont_frags)
        {
          have_data = true;
          if (! _felixProcContNCFrags(cont_frags, n_felix_frags, true, evt, raw_digits, timestamps, apalist, inputLabel))
            {
              return false;
            }
        }
    }
  else
    {
      auto frags = evt.getHandle<artdaq::Fragments>(inputLabel);

      if (frags)
        {
          have_data_nc = true;
          if (! _felixProcContNCFrags(frags, n_felix_frags, false, evt, raw_digits, timestamps, apalist, inputLabel))
            {
              return false;
            }
        }
    }

  // returns true if we want to add to the number of fragments processed. 

  return have_data || have_data_nc;
}

bool IcebergDataInterface::_felixProcContNCFrags(art::Handle<artdaq::Fragments> frags, 
                                                 size_t &n_felix_frags, 
                                                 bool is_container, 
                                                 art::Event &evt, 
                                                 RawDigits& raw_digits, 
                                                 RDTimeStamps &timestamps, 
                                                 std::vector<int> &apalist,
                                                 std::string inputLabel)
{
  uint32_t runNumber = evt.run();

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
          else
            {
              _KeptCorruptData = true;
            }
        }
      if (process_flag)
        {
          if (is_container)
            {
              artdaq::ContainerFragment cont_frag(frag);
              for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
                {
                  if (_process_FELIX_AUX(evt,*cont_frag[ii], raw_digits, timestamps, apalist, runNumber, inputLabel)) ++n_felix_frags;
                }
            }
          else
            {
              if (_process_FELIX_AUX(evt,frag, raw_digits, timestamps, apalist, runNumber, inputLabel)) ++n_felix_frags;
            }
        }
    }
  evt.removeCachedProduct(frags);
  return true;
}


bool IcebergDataInterface::_process_FELIX_AUX(art::Event &evt,
                                              const artdaq::Fragment& frag, RawDigits& raw_digits,
                                              RDTimeStamps &timestamps,
                                              std::vector<int> &apalist,
                                              uint32_t runNumber,
                                              std::string inputLabel)
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

  art::ServiceHandle<dune::IcebergChannelMapService> channelMap;

  // Load overlay class.   Either a felix or a frame14 overlay, depending on the
  // input instance name

  bool is14 = (inputLabel.find("FRAME14") != std::string::npos);
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
      //   << (int) crate << " " << (int) slot << " " << (int) fiber << " " << (int) frag.fragmentID() << std::endl;
      if (frame_version == 0) 
        {
          std::cout << "ICEBERG frame_version = 0; skipping this fragment" << std::endl;
          return false;
        }
      fiber ++;  // read in 0 to 1, go from 1 to 2
      crate = 1;  // ignored anyway for Iceberg -- just one crate.
    }
  else
    {
      crate = felixptr->crate_no(0);
      slot  = felixptr->slot_no(0);
      fiber = felixptr->fiber_no(0); // decode this one later 
    }

  if (slot > 4) 
    {
      if (_felix_drop_frags_with_badsf)  
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
      if ( ( (apalist[ialist] == -1)    && (!_rce_drop_frags_with_badc || (crate >0 && crate < 7)) )      || 
           (apalist[ialist] == (int) crate) || 
           (apalist[ialist] == 7 && (crate == 0 || crate > 6)) )
        {
          foundapainlist = true;
          break;
        }
    }
  if (!foundapainlist) return true;

  //std::cout << "FELIX raw decoder trj: " << (int) crate << " " << (int) slot << " " << (int) fiber << std::endl;

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
  if (n_frames ==0) return true;

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


    unsigned int offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crateloc, slotloc2, fiberloc2, chloc, dune::IcebergChannelMapService::kFELIX); 

    // skip this channel if we are asked to.

    if (_max_offline_channel >= 0 && _min_offline_channel >= 0 && _max_offline_channel >= _min_offline_channel && 
        (offlineChannel < (size_t) _min_offline_channel || offlineChannel > (size_t) _max_offline_channel) ) continue;

    v_adc.clear();
    std::vector<dune::adc_t> waveform( is14 ? 
                                       frame14ptr->get_ADCs_by_channel(ch) : 
                                       felixptr->get_ADCs_by_channel(ch) );

    for(unsigned int nframe=0;nframe<waveform.size();nframe++){
      v_adc.push_back(waveform.at(nframe));  
    }

    if ( v_adc.size() != _full_tick_count)
      {
        if (_enforce_full_tick_count)
          {
            MF_LOG_WARNING("_process_FELIX_AUX:") << "Nticks not the required value: " << v_adc.size() << " " 
                                                  << _full_tick_count << " Discarding Data";
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

    raw::RDTimeStamp rdtimestamp( is14 ? frame14ptr->timestamp() : felixptr->timestamp(),offlineChannel);
    timestamps.push_back(rdtimestamp);
  }

  return true;
}


// compute median and sigma.  

void IcebergDataInterface::computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma)
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

DEFINE_ART_CLASS_TOOL(IcebergDataInterface)
