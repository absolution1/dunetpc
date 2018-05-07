////////////////////////////////////////////////////////////////////////////
// file: RawData311InputDriver.cc
//
// brief Source to convert raw binary file from 3x1x1 to root file useful to LArSoft
// Adapted from LarRawInputDriverShortBo.cc
//
// Created: April 20th 2017, Last Modified:
// Author: Kevin Fusshoeller, kevin.fusshoeller@cern.ch
////////////////////////////////////////////////////////////////////////////

#include "RawData311InputDriver.h"

#include "larcore/Geometry/Geometry.h"

#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "EventDecoder.h"
#include "dlardaq.h"

#include <iostream>
#include <ios>

// ---------------------------------------------------------------------------------------
// 311 DAQ interface

namespace lris
{
  void SplitAdc(const std::vector<dlardaq::adc16_t> *adc, size_t channel, uint32_t num_samples,
		 std::vector<short> &adclist)
  {
    for(uint32_t i = 0; i < num_samples; i++)
    {
      adclist.emplace_back(adc->at(channel*num_samples + i));
    }
  }// SplitAdc


  //-------------------------------
  size_t Get311Chan(size_t LAr_chan){
    size_t crate = LAr_chan / 320;
    size_t Chan311;

    LAr_chan = 8*(LAr_chan/8+1)-LAr_chan%8 -1;

    if(crate == 0)
      {
	LAr_chan = 32*(LAr_chan/32+1)-LAr_chan%32 -1;
        size_t card = 4 - ((LAr_chan / 32) % 5);
        if(LAr_chan > 159)
          {
            size_t shift = 31 - (LAr_chan % 32);
            Chan311 = (2*card)*32 + shift;
          }
        else
          {
            size_t shift = 31 - (LAr_chan % 32);
            Chan311 = (2*card + 1)*32 + shift;
          }
      }
    else
      {
        size_t new_LAr_chan = LAr_chan - crate*320;
        size_t card = ((new_LAr_chan / 32) % 5);
        if(new_LAr_chan > 159)
          {
            size_t shift = new_LAr_chan % 32;
            Chan311 = (2*card)*32 + shift;
          }
        else
          {
            size_t shift = new_LAr_chan % 32;
            Chan311 = (2*card + 1)*32 + shift;
          }
        Chan311 = Chan311 + crate*320;
      } // end of if/else statementi

    return Chan311;
  } // Get311Chan


  // ----------------------------------------------------------------------
  //
  // ----------------------------------------------------------------------


  void ReadPedestalFile(std::string PedestalFileName, std::vector< std::pair<double, double> > &PedMap){
  //initialize the channel-ped value map
    std::ifstream file;
    file.open(PedestalFileName);
    if( !file.is_open() )
    {
      throw art::Exception( art::errors::FileReadError )
		<< "failed to open input file " << PedestalFileName << "\n";
    }

    while(!file.eof())
    {
      size_t ch, cryo, crate, rawch;
      double mean, rms;
      file >> rawch >> cryo >> crate >> ch >> mean >> rms;
      PedMap.emplace_back(mean, rms);
    }

    file.close();
    return;
  }//Read Pedestal File()


  // ---------------------------------------------------------------------
  //
  // ---------------------------------------------------------------------


  void RawData311InputDriver::process_Event311(std::vector<raw::RawDigit>& digitList,
			   dlardaq::evheader_t &event_head,
			   uint16_t evt_num)
  {
    // one digit for every wire on each plane
    digitList.clear();
    digitList.resize(nchannels);
    // Get the data.
    std::vector<dlardaq::adc16_t> ADCvec311;
    DataDecode.GetEvent(evt_num, event_head, ADCvec311);
    std::vector<dlardaq::adc16_t> *ADCvec311Pointer = &ADCvec311;
    // fill the wires
    std::vector<short> adclist;

    for(size_t LAr_chan = 0; LAr_chan < (size_t)nchannels; LAr_chan++)
    {
      adclist.clear();
      size_t Chan311 = Get311Chan(LAr_chan);
      SplitAdc(ADCvec311Pointer, Chan311, nsamples, adclist);
      short unsigned int nTickReadout = nsamples;
      raw::ChannelID_t channel = LAr_chan;
      raw::Compress_t comp = raw::kNone;
      raw::RawDigit rd(channel, nTickReadout, adclist, comp);

      double pedval = RawData311InputDriver::GetPedMean(Chan311, &fPedMap);
      //std::cout << "Pedval: " << pedval << "\n";
      double pedrms = RawData311InputDriver::GetPedRMS(Chan311, &fPedMap);
      rd.SetPedestal(pedval, pedrms);

      digitList[LAr_chan] = rd;
    }
  }// process_Event311


  //------------------------------------------------------------------
  // class c'tor/d'tor
  RawData311InputDriver::RawData311InputDriver(fhicl::ParameterSet const &p,
					   art::ProductRegistryHelper &helper,
					   art::SourceHelper const &pm)
    :
    fSourceHelper(pm),
    fCurrentSubRunID(),
    fEventCounter(0),
    DataDecode(nchannels, nsamples)
  {
    fPedestalFile = p.get<std::string>("PedestalFile");
    helper.reconstitutes<std::vector<raw::RawDigit>, art::InEvent>("daq");
  }


  // Close File.
  void RawData311InputDriver::closeCurrentFile()
  {
    mf::LogInfo(__FUNCTION__)<<"File boundary: processed " <<fEventCounter <<" events out of " <<fNEvents <<"\n";
    DataDecode.m_file.close();
  }


  // Read File.
  void RawData311InputDriver::readFile(std::string const &name,
				     art::FileBlock* &fb)
  {
    // Read in the pedestal file
    ReadPedestalFile(fPedestalFile, fPedMap);

    filename = name;
    // Fill and return a new Fileblock
    fb = new art::FileBlock(art::FileFormatVersion(1, "311 RawInput 2017"), name);

    // Prepare the EventDecoder
    //short nchannels_Vplane = 320;
    //short nchannels_Zplane = 960;
    //short nchannels = nchannels_Vplane + nchannels_Zplane;
    //uint32_t nsamples = 1667;
    //DataDecode(nchannels, nsamples);

    DataDecode.m_file.open(name.c_str(), std::ios_base::in | std::ios_base::binary);
    if( !DataDecode.m_file.is_open() )
    {
      throw art::Exception( art::errors::FileReadError )
	<< "failed to open input file " << name << "\n";
    }

    // Read in the file header.
    std::vector<dlardaq::BYTE> head_buf;
    head_buf.resize(dlardaq::RunHeadSz);

    DataDecode.m_file.read(&head_buf[0], head_buf.size());
    dlardaq::decode_runhead(&head_buf[0], file_head);

    // Define start of event data.
    std::streampos data_start = DataDecode.m_file.tellg();

    // Read in the file footer.
    std::vector<dlardaq::BYTE> foot_buf;
    foot_buf.resize(dlardaq::FileFootSz);

    DataDecode.m_file.seekg(-dlardaq::FileFootSz, DataDecode.m_file.end);
    DataDecode.m_file.read(&foot_buf[0], foot_buf.size());
    dlardaq::decode_filefoot(&foot_buf[0], file_foot);
    DataDecode.m_file.seekg(data_start);

    fNEvents = file_foot.num_events;
    if(fNEvents > 0 && fNEvents < 1000) //There should be 335 events at most in one file.
    {
      mf::LogInfo("")<<"Opened file " <<name <<" with " << fNEvents <<" events." <<"\n";
    } else
    {
      throw art::Exception( art::errors::FileReadError )
	<<"File " <<name <<" seems to have too many events: " <<fNEvents <<"\n";
    }
  }


  // Read event.
  bool RawData311InputDriver::readNext(art::RunPrincipal* const &/*inR*/,
				     art::SubRunPrincipal* const &/*inSR*/,
				     art::RunPrincipal* &outR,
				     art::SubRunPrincipal* &outSR,
				     art::EventPrincipal* &outE)
  {
    if(fEventCounter == fNEvents)
    {
      mf::LogInfo(__FUNCTION__)<<"All the files have been read in. Checking end of file..." <<"\n";
      std::streampos current_position = DataDecode.m_file.tellg();
      DataDecode.m_file.seekg(0, std::ios::end);
      std::streampos file_length = DataDecode.m_file.tellg();
      if( ((uint8_t)file_length - (uint8_t)current_position) > (uint8_t)100 )
      {
	throw art::Exception( art::errors::FileReadError )
	  <<"Processed " <<fEventCounter <<" events out of " <<fNEvents <<" but there are still "
	  <<(file_length - current_position) <<" bits left." <<"\n";
      }
      mf::LogInfo(__FUNCTION__)<<"Completed reading file and closing output file." <<"\n";
      return false; //Tells readNext that all events have been read.
    }

    mf::LogInfo(__FUNCTION__)<<"Reading event " << fEventCounter << " from " << fNEvents <<"\n";

    // Create empty result, then fill it from current file
    dlardaq::evheader_t event_head;
    std::unique_ptr< std::vector<raw::RawDigit> > tpc_raw_digits( new std::vector<raw::RawDigit> );
    process_Event311(*tpc_raw_digits, event_head, fEventCounter++);



    // Prepare some stuff
    art::RunNumber_t rn;
    rn = file_head.run_num;

    std::uint64_t tthi = event_head.trig_info.ts.tv_sec;
    std::uint64_t thilo = (tthi << 32) + event_head.trig_info.ts.tv_nsec;
    art::Timestamp tstamp(thilo);

    size_t pos = filename.find("-");
		size_t end = filename.find(".");
    std::string str_subrn = filename.substr(pos+1, end-1);
    int int_subrn = std::stoi(str_subrn);
    art::SubRunNumber_t subrn = int_subrn;
    art::SubRunID newID(rn, subrn);
    if( fCurrentSubRunID.runID() != newID.runID() )
    {
      outR = fSourceHelper.makeRunPrincipal(rn, tstamp);
    }
    if( fCurrentSubRunID != newID )
    {
      outSR = fSourceHelper.makeSubRunPrincipal(rn, subrn, tstamp);
      fCurrentSubRunID = newID;
    }

    outE = fSourceHelper.makeEventPrincipal(fCurrentSubRunID.run(), fCurrentSubRunID.subRun(),
					    fEventCounter - 1, tstamp);

    // Put products in the event.
    art::put_product_in_principal(std::move(tpc_raw_digits), *outE, "daq");

    return true;
   }


} // namespace lris
