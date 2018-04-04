////////////////////////////////////////////////////////////////////////
// Class:       PDSPTPCRawDecoder
// Plugin Type: producer (art v2_10_03)
// File:        PDSPTPCRawDecoder_module.cc
//
// Generated at Fri Mar  2 15:36:20 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_02_00.
// Original code from Jingbo Wang in separate RCE and FELIX raw decoders
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

#include <memory>

// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/RceFragment.hh"
#include "dune-raw-data/Overlays/FelixFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "dune-raw-data/Overlays/FragmentType.hh"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

// larsoft includes
#include "lardataobj/RawData/RawDigit.h"

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

  std::string _rce_input_label; 
  std::string _rce_input_container_instance;
  std::string _rce_input_noncontainer_instance;
  std::string _felix_input_label; 
  std::string _felix_input_container_instance;
  std::string _felix_input_noncontainer_instance;
  int _rce_fragment_type;
  //int _felix_fragment_type; // unused

  std::string _output_label;
  bool _expect_rce_container_fragments;  
  bool _expect_felix_container_fragments;
  bool _drop_RCE_damaged_data;  
  bool _drop_FELIX_damaged_data;

  bool _processRCE(art::Event &evt, RawDigits& raw_digits);
  bool _processFELIX(art::Event &evt, RawDigits& raw_digits);
  bool _process_RCE_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits);
  bool _process_FELIX_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits);

 std::vector<uint16_t> _buffer;
};


PDSPTPCRawDecoder::PDSPTPCRawDecoder(fhicl::ParameterSet const & p)
{
  _rce_input_label = p.get<std::string>("RCERawDataLabel","daq");
  _rce_input_container_instance = p.get<std::string>("RCERawDataContainerInstance","ContainerTPC");
  _rce_input_noncontainer_instance = p.get<std::string>("RCERawDataNonContainerInstance","TPC");
  _rce_fragment_type = p.get<int>("RCEFragmentType",2);

  _felix_input_label = p.get<std::string>("FELIXRawDataLabel");
  _felix_input_container_instance = p.get<std::string>("FELIXRawDataContainerInstance","ContainerFELIX");
  _felix_input_noncontainer_instance = p.get<std::string>("FELIXRawDataNonContainerInstance","FELIX");
  _rce_fragment_type = p.get<int>("FELIXFragmentType",2);

  _output_label = p.get<std::string>("OutputDataLabel");
  _expect_rce_container_fragments = p.get<bool>("ExpectRCEContainerFragments", true);
  _expect_felix_container_fragments = p.get<bool>("ExpectFELIXContainerFragments", false);
  _drop_RCE_damaged_data = p.get<bool>("DropRCEDamagedData", true);
  _drop_FELIX_damaged_data = p.get<bool>("DropFELIXDamagedData", true);

  produces<RawDigits>( _output_label ); 
}

void PDSPTPCRawDecoder::produce(art::Event &e)
{
  RawDigits raw_digits;
  _processRCE(e,raw_digits);
  _processFELIX(e,raw_digits);
  e.put(std::make_unique<decltype(raw_digits)>(std::move(raw_digits)),
          _output_label);
}

bool PDSPTPCRawDecoder::_processRCE(art::Event &evt, RawDigits& raw_digits)
{
  size_t n_rce_frags = 0;
  if (_expect_rce_container_fragments) {
    art::Handle<artdaq::Fragments> cont_frags;
    evt.getByLabel(_rce_input_label, _rce_input_container_instance, cont_frags);  // hardcoded label .. maybe fix
    try { cont_frags->size(); }
    catch(std::exception e) {
      LOG_DEBUG("_processRCE") << "Container TPC/RCE data not found " 
                                       << "Run: " << evt.run()
		                       << ", SubRun: " << evt.subRun()
				       << ", Event: " << evt.event();
      return false;
    }
    //Check that the data is valid
    if(!cont_frags.isValid()){
      LOG_ERROR("_processRCE") << "Container TPC/RCE fragments Not Valid " 
                                       << "Run: " << evt.run()
		                       << ", SubRun: " << evt.subRun()
				       << ", Event: " << evt.event();
      return false;
    }
    
    for (auto const& cont : *cont_frags)
      {
	artdaq::ContainerFragment cont_frag(cont);
	for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
	  {
	    if (_process_RCE_AUX(*cont_frag[ii], raw_digits)) ++n_rce_frags;
	  }
      }
  }
  else
    {
      art::Handle<artdaq::Fragments> frags;
      evt.getByLabel(_rce_input_label, _rce_input_noncontainer_instance, frags);  // hardcoded label?
      try { frags->size(); }
      catch(std::exception e) {
	LOG_DEBUG("_process_RCE_AUX") << "TPC/RCE fragment data not found " 
					 << "Run: " << evt.run()
					 << ", SubRun: " << evt.subRun()
					 << ", Event: " << evt.event();
	return false;
      }

      if(!frags.isValid()){
	LOG_ERROR("_process_RCE_AUX") << "TPC/RCE fragments Not Valid " 
					 << "Run: " << evt.run()
					 << ", SubRun: " << evt.subRun()
					 << ", Event: " << evt.event();
	return false;
      }

      for(auto const& frag: *frags)
	{
	  if (_process_RCE_AUX(frag, raw_digits)) ++n_rce_frags;
	}
    }

  LOG_INFO("_processRCE_AUX")
    << " Processed " << n_rce_frags
    << " RCE Fragments, making "
    << raw_digits.size()
    << " RawDigits.";
  return true;
}

bool PDSPTPCRawDecoder::_process_RCE_AUX(
					   const artdaq::Fragment& frag, 
					   RawDigits& raw_digits
					   )
{
  // FIXME: Remove hard-coded fragment type
  if((unsigned)frag.type() != 2) return false;

  LOG_INFO("_Process_RCE_AUX")
    << "   SequenceID = " << frag.sequenceID()
    << "   fragmentID = " << frag.fragmentID()
    << "   fragmentType = " << (unsigned)frag.type()
    << "   Timestamp =  " << frag.timestamp();
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;
  dune::RceFragment rce(frag);
  
  uint32_t ch_counter = 0;
  for (int i = 0; i < rce.size(); ++i)
    {
      auto const * rce_stream = rce.get_stream(i);
      int n_ch = rce_stream->getNChannels();
      int n_ticks = rce_stream->getNTicks();
      auto const identifier = rce_stream->getIdentifier();
      uint32_t crateNumber = identifier.getCrate();
      uint32_t slotNumber = identifier.getSlot();
      uint32_t fiberNumber = identifier.getFiber();

      LOG_INFO("_Process_RCE_AUX")
	<< "RceFragment timestamp: " << rce_stream->getTimeStamp()
	<< ", NChannels: " << n_ch
	<< ", NTicks: " << n_ticks;

      // TODO -- speed this up!!  Remove one buffer copy

      size_t buffer_size = n_ch * n_ticks;
      if (_buffer.capacity() < buffer_size)
	{
	  LOG_INFO("_process_RCE_AUX")
	    << "Increase buffer size from " << _buffer.capacity()
	    << " to " << buffer_size;

	  _buffer.reserve(buffer_size);
	}

      uint16_t* adcs = _buffer.data();
      rce_stream->getMultiChannelData(adcs);

      raw::RawDigit::ADCvector_t v_adc;
      for (int i_ch = 0; i_ch < n_ch; i_ch++)
	{
	  if(i==0 && i_ch ==0) std::cout<<" ADCs for the the 100 ticks in the 1st channel of the 1st RCE "<<std::endl;
	  v_adc.clear();
	  for (int i_tick = 0; i_tick < n_ticks; i_tick++)
	    {
	      v_adc.push_back(adcs[i_tick]);
	    }
	  adcs += n_ticks;

	  ch_counter++;
	  int offlineChannel = -1;
	  offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crateNumber, slotNumber, fiberNumber, i_ch);
	  raw::RawDigit raw_digit(offlineChannel, n_ticks, v_adc);
	  raw_digits.push_back(raw_digit);
                
	}
    }

  return true;
}

bool PDSPTPCRawDecoder::_processFELIX(art::Event &evt, RawDigits& raw_digits)
{

    // TODO Use LOG_DEBUG
    LOG_INFO("_processFELIX")
      << "-------------------- FELIX RawDecoder -------------------";

  unsigned int n_felix_frags = 0;  

  if (_expect_felix_container_fragments) {
    art::Handle<artdaq::Fragments> cont_frags;
    evt.getByLabel(_felix_input_label, _felix_input_container_instance, cont_frags);  // TODO -- un-hardwire this label

    try { cont_frags->size(); }
    catch(std::exception e) {
      LOG_DEBUG("_processFELIX") << "Container TPC/FELIX data not found " 
                                       << "Run: " << evt.run()
		                       << ", SubRun: " << evt.subRun()
				       << ", Event: " << evt.event();
      return false;
    }
    //Check that the data is valid
    if(!cont_frags.isValid()){
      LOG_ERROR("_processFELIX") << "Container TPC/FELIX fragments Not Valid " 
                                       << "Run: " << evt.run()
		                       << ", SubRun: " << evt.subRun()
				       << ", Event: " << evt.event();
      return false;
    }
    
    for (auto const& cont : *cont_frags)
    {
      artdaq::ContainerFragment cont_frag(cont);
      for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
	{
          if (_process_FELIX_AUX(*cont_frag[ii], raw_digits)) ++n_felix_frags;
	}
    }
  }
  else
  {
    art::Handle<artdaq::Fragments> frags;
    evt.getByLabel(_felix_input_label, _felix_input_noncontainer_instance, frags);
    try { frags->size(); }
    catch(std::exception e) {
	LOG_DEBUG("_process_FELIX_AUX") << "TPC/FELIX fragment data not found " 
					 << "Run: " << evt.run()
					 << ", SubRun: " << evt.subRun()
					 << ", Event: " << evt.event();
	return false;

    }

    //Check that the data is valid
    if(!frags.isValid()){
	LOG_ERROR("_process_FELIX_AUX") << "TPC/FELIX fragments Not Valid " 
					 << "Run: " << evt.run()
					 << ", SubRun: " << evt.subRun()
					 << ", Event: " << evt.event();
	return false;
    }

    for(auto const& frag: *frags)
    {
      if (_process_FELIX_AUX(frag, raw_digits)) ++n_felix_frags;
    }
  }

  LOG_INFO("_processFELIX")
      << " Processed " << n_felix_frags
      << " FELIX Fragments, total size of raw digits is now "
      << raw_digits.size()
      << " RawDigits.";

  return true;
}

bool PDSPTPCRawDecoder::_process_FELIX_AUX(const artdaq::Fragment& frag, RawDigits& raw_digits)
{
  // FIXME: Remove hard-coded fragment type
  //if((unsigned)frag.type() != 2) return false;

  LOG_INFO("_process_FELIX_AUX")
      << "   SequenceID = " << frag.sequenceID()
      << "   fragmentID = " << frag.fragmentID()
      << "   fragmentType = " << (unsigned)frag.type()
      << "   Timestamp =  " << frag.timestamp();
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;
  //Load overlay class.
  dune::FelixFragment felix(frag);
  //Get detector element number
  uint8_t crate = felix.crate_no(0);
  uint8_t slot = felix.slot_no(0);
  uint8_t fiber = felix.fiber_no(0); // two numbers? 
  const unsigned n_frames = felix.total_frames(); // One frame contains 20 ticks.
  //std::cout<<" Nframes = "<<n_frames<<std::endl;
  //_h_nframes->Fill(n_frames);
  const unsigned n_channels = dune::FelixFrame::num_ch_per_frame;// 256

  // check optimization of this -- size not reserved

  raw::RawDigit::ADCvector_t v_adc;
  //v_adc.reserve(n_frames*n_channels);
  // Fill the adc vector.  

  //typedef std::tuple<uint8_t, uint8_t, uint8_t, unsigned> WireInfo_tuple; // unused
  for(unsigned ch = 0; ch < n_channels; ++ch) {
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
    int offlineChannel = -1;
    offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crate, slot, fiber, ch); // FIXME
    // Push to raw_digits.
    raw::RawDigit raw_digit(offlineChannel, n_frames, v_adc);
    raw_digits.push_back(raw_digit);
  }
  return true;
}

DEFINE_ART_MODULE(PDSPTPCRawDecoder)
