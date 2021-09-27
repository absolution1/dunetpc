////////////////////////////////////////////////////////////////////////
// Class:       IcebergFELIXBufferDecoderMarch2021
// Plugin Type: producer (art v2_10_03)
// File:        IcebergFELIXBufferDecoderMarch2021_module.cc
//
// Generated at Fri Mar  2 15:36:20 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_02_00.  
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art_root_io/TFileService.h"

#include <memory>
#include <cmath>

// ROOT includes
#include "TMath.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Services/ChannelMap/IcebergChannelMapService.h"

// larsoft includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/raw.h"

// DUNE includes
#include "dune/DuneObj/RDStatus.h"

#include <stdio.h>

class IcebergFELIXBufferDecoderMarch2021 : public art::EDProducer {

public:
  explicit IcebergFELIXBufferDecoderMarch2021(fhicl::ParameterSet const & p);
  IcebergFELIXBufferDecoderMarch2021(IcebergFELIXBufferDecoderMarch2021 const &) = delete;
  IcebergFELIXBufferDecoderMarch2021(IcebergFELIXBufferDecoderMarch2021 &&) = delete;
  IcebergFELIXBufferDecoderMarch2021 & operator = (IcebergFELIXBufferDecoderMarch2021 const &) = delete;
  IcebergFELIXBufferDecoderMarch2021 & operator = (IcebergFELIXBufferDecoderMarch2021 &&) = delete;
  void produce(art::Event & e) override;

private:
  typedef std::vector<raw::RawDigit> RawDigits;
  typedef std::vector<raw::RDTimeStamp> RDTimeStamps;
  typedef art::Assns<raw::RawDigit,raw::RDTimeStamp> RDTsAssocs;
  typedef art::PtrMaker<raw::RawDigit> RDPmkr;
  typedef art::PtrMaker<raw::RDTimeStamp> TSPmkr;
  typedef std::vector<raw::RDStatus> RDStatuses;

  // open files

  std::vector<FILE*> fInputFilePointers;

  // configuration parameters

  std::vector<std::string>   fInputFiles; 
  size_t                     fNSamples;
  std::string                fOutputLabel;
  bool                       fCompressHuffman;
  ULong64_t                  fDesiredStartTimestamp;
  bool                       fFirstRead;

  void computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma);
  // Converts 14 bit packed channel data (56 uint32 words from the WIB) to byte-aligned 16 bit arrays (128 uint16 values)
  void unpack14(const uint32_t *packed, uint16_t *unpacked);
};


IcebergFELIXBufferDecoderMarch2021::IcebergFELIXBufferDecoderMarch2021(fhicl::ParameterSet const & p)
  : EDProducer(p)
{
  fInputFiles = p.get<std::vector<std::string>>("InputFiles");
  fNSamples = p.get<size_t>("NSamples",2000);
  fOutputLabel = p.get<std::string>("OutputDataLabel","daq");
  fCompressHuffman = p.get<bool>("CompressHuffman",false);
  fDesiredStartTimestamp = p.get<ULong64_t>("StartTimestamp",0);

  produces<RawDigits>( fOutputLabel ); //the strings in <> are the typedefs defined above
  produces<RDTimeStamps>( fOutputLabel );
  produces<RDTsAssocs>( fOutputLabel );
  produces<RDStatuses>( fOutputLabel );

  fInputFilePointers.clear();
  for (size_t ifile=0; ifile<fInputFiles.size(); ++ifile)
    {
      fInputFilePointers.push_back(fopen(fInputFiles.at(ifile).data(),"r"));
    }
  fFirstRead = true;
}

void IcebergFELIXBufferDecoderMarch2021::produce(art::Event &e)
{

  art::ServiceHandle<dune::IcebergChannelMapService> channelMap;
  RawDigits raw_digits;
  RDTimeStamps rd_timestamps;
  RDTsAssocs rd_ts_assocs;

  RDPmkr rdpm(e,fOutputLabel);
  TSPmkr tspm(e,fOutputLabel);

  bool discard_data = false;

  uint32_t framebuf[117];
  uint16_t databuf[128];
  uint64_t timestampstart=0;
  uint64_t timestamp=0;

  size_t nfiles = fInputFiles.size();

  // search for the first timestamp on the first read
  // don't do on subsequent reads so we don't always read in a frame just to see what
  // the timestamp is.

  if (fFirstRead)
    {
      for (size_t ifile=0; ifile < nfiles; ++ ifile)
        {
          do
            {
              fread(framebuf, sizeof(uint32_t), 117, fInputFilePointers.at(ifile));
              if (feof(fInputFilePointers.at(ifile)))
                {
                  // don't handle this too gracefully at the moment
                  throw cet::exception("IcebergFELXIBufferDecoderMarch2021") <<
                    "Attempt to read off the end of file " << fInputFiles.at(ifile);
                }
              timestamp= framebuf[3];
              timestamp <<= 32;
              timestamp += framebuf[2];
            }
          while (timestamp+47 < fDesiredStartTimestamp);  
          // criterion so the next timestamp (+32) will be the one we want
        }
      fFirstRead = false;
    }

  // align the readin

  std::vector<std::vector<uint32_t>> fbcache(nfiles);
  std::vector<uint64_t> tscache;
  uint64_t latest_timestamp=0;

  // read one frame in from each file to see what the latest timestamp is

  for (size_t ifile=0; ifile<nfiles; ++ifile)
    {
      fbcache.at(ifile).resize(117);
      fread(fbcache.at(ifile).data(), sizeof(uint32_t), 117, fInputFilePointers.at(ifile));
      if (feof(fInputFilePointers.at(ifile)))
        {
          // don't handle this too gracefully at the moment
          throw cet::exception("IcebergFELXIBufferDecoderMarch2021") <<
            "Attempt to read off the end of file " << fInputFiles.at(ifile);
        }
      timestamp= fbcache.at(ifile)[3];
      timestamp <<= 32;
      timestamp += fbcache.at(ifile)[2];
      tscache.push_back(timestamp);
      if (timestamp > latest_timestamp)
        {
          latest_timestamp = timestamp;
        }
    }

  // read enough frames from the other files so that we align the frames to +- 16 ticks

  for (size_t ifile=0; ifile<nfiles; ++ifile)
    {
      while (tscache.at(ifile) + 16 < latest_timestamp)
        {
          fread(fbcache.at(ifile).data(), sizeof(uint32_t), 117, fInputFilePointers.at(ifile));
          if (feof(fInputFilePointers.at(ifile)))
            {
              // don't handle this too gracefully at the moment
              throw cet::exception("IcebergFELXIBufferDecoderMarch2021") <<
                "Attempt to read off the end of file " << fInputFiles.at(ifile);
            }
          timestamp= fbcache.at(ifile)[3];
          timestamp <<= 32;
          timestamp += fbcache.at(ifile)[2];
          tscache.at(ifile) = timestamp;
        }
    }  

  for (size_t ifile=0; ifile < nfiles; ++ ifile)
    {
      int slot = 0;
      int fiber = 0;

      std::vector<raw::RawDigit::ADCvector_t> adcvv(256);
      for (size_t itick=0; itick<fNSamples; ++itick)
        {
          if (itick == 0)  // already read in the first frame
            {
              for (size_t i=0; i<117; ++i)
                {
                  framebuf[i] = fbcache.at(ifile).at(i);
                }
            }
          else
            {
              fread(framebuf, sizeof(uint32_t), 117, fInputFilePointers.at(ifile));
              if (feof(fInputFilePointers.at(ifile)))
                {
                  // don't handle this too gracefully at the moment
                  throw cet::exception("IcebergFELXIBufferDecoderMarch2021") <<
                    "Attempt to read off the end of file " << fInputFiles.at(ifile);
                }
            }
          int curslot = (framebuf[0] & 0x7000) >> 12;   // assume these are all the same
          int curfiber = (framebuf[0] & 0x8000) >> 15;
          if (itick>0)
            {
              if (curslot != slot)
                {
                  throw cet::exception("IcebergFELXIBufferDecoderMarch2021") <<
                    "Slot mismatch in file: " << curslot << " " << slot;
                }
              if (curfiber != fiber)
                {
                  throw cet::exception("IcebergFELXIBufferDecoderMarch2021") <<
                    "Fiber mismatch in file: " << curfiber << " " << fiber;
                }
            }
          else
            {
              slot = curslot;
              fiber = curfiber;
            }

          uint64_t timestamp= framebuf[3];
          timestamp <<= 32;
          timestamp += framebuf[2];
          //std::cout << std::dec << "   Slot: " << slot << " Fiber: " << fiber 
          //        << " Timestamp: " << std::dec << timestamp <<std::endl;
          if (itick == 0) 
            {
              timestampstart = timestamp;
            }

          // do the data-rearrangement transpose

          unpack14(&(framebuf[4]),databuf);
          for (size_t ichan=0; ichan<128; ++ichan)
            {
              adcvv.at(ichan).push_back(databuf[ichan]);
            }
          unpack14(&(framebuf[4+56]),databuf);
          for (size_t ichan=0; ichan<128; ++ichan)
            {
              adcvv.at(ichan+128).push_back(databuf[ichan]);
            }
        }

      for (size_t ichan=0; ichan<256; ++ichan)
        {
          float median=0;
          float sigma=0;
          computeMedianSigma(adcvv.at(ichan),median,sigma);

          // handle 256 channels on two fibers -- use the channel map that assumes 128 chans per fiber (=FEMB)
    
          unsigned int fiberloc = 0;
          if (fiber == 0) 
            {
              fiberloc = 1;
            }
          else if (fiber == 1)
            {
              fiberloc = 3;
            }

          unsigned int chloc = ichan;
          if (chloc > 127)
            {
              chloc -= 128;
              fiberloc++;
            }
          //unsigned int crateloc = crate;  


          // inverted ordering on back side, Run 2c (=Run 3)
          // note Shekhar's FEMB number is fiber-1, and WIB is slot+1

          auto slotloc2 = slot;
          auto fiberloc2 = fiberloc;

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

          size_t uncompressed_nticks = fNSamples;  
          raw::Compress_t cflag=raw::kNone;
          if (fCompressHuffman)
            {
              cflag = raw::kHuffman;
              raw::Compress(adcvv.at(ichan),cflag);
            }

          raw::RawDigit raw_digit(offlineChannel, uncompressed_nticks, adcvv.at(ichan), cflag);
          raw_digit.SetPedestal(median,sigma);
          raw_digits.push_back(raw_digit);  

          raw::RDTimeStamp rdtimestamp(timestampstart,offlineChannel);
          rd_timestamps.push_back(rdtimestamp);

          //associate the raw digit and the timestamp data products
          auto const rawdigitptr = rdpm(raw_digits.size()-1);
          auto const rdtimestampptr = tspm(rd_timestamps.size()-1);
          rd_ts_assocs.addSingle(rawdigitptr,rdtimestampptr);            

        }
    }

  if (discard_data)
    {
      RawDigits empty_raw_digits;
      RDTimeStamps empty_rd_timestamps;
      RDTsAssocs empty_rd_ts_assocs;
      RDStatuses statuses;
      statuses.emplace_back(true,false,1);
      e.put(std::make_unique<decltype(empty_raw_digits)>(std::move(empty_raw_digits)),fOutputLabel);
      e.put(std::make_unique<decltype(empty_rd_timestamps)>(std::move(empty_rd_timestamps)),fOutputLabel);
      e.put(std::make_unique<decltype(empty_rd_ts_assocs)>(std::move(empty_rd_ts_assocs)),fOutputLabel);
      e.put(std::make_unique<decltype(statuses)>(std::move(statuses)),fOutputLabel);
    }
  else
    {
      RDStatuses statuses;
      unsigned int statword=0;
      statuses.emplace_back(false,false,statword);
      e.put(std::make_unique<decltype(raw_digits)>(std::move(raw_digits)),fOutputLabel);
      e.put(std::make_unique<decltype(rd_timestamps)>(std::move(rd_timestamps)),fOutputLabel);
      e.put(std::make_unique<decltype(rd_ts_assocs)>(std::move(rd_ts_assocs)),fOutputLabel);
      e.put(std::make_unique<decltype(statuses)>(std::move(statuses)),fOutputLabel);
    }
}


// compute median and sigma.  Sigma is half the distance between the upper and lower bounds of the
// 68% region where 34% is above the median and 34% is below ("centered" on the median).

void IcebergFELIXBufferDecoderMarch2021::computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma)
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

void IcebergFELIXBufferDecoderMarch2021::unpack14(const uint32_t *packed, uint16_t *unpacked) {
  for (size_t i = 0; i < 128; i++) { // i == n'th U,V,X value
    const size_t low_bit = i*14;
    const size_t low_word = low_bit / 32;
    const size_t high_bit = (i+1)*14-1;
    const size_t high_word = high_bit / 32;
    //std::cout << "low_word, high_word: " << low_word << " " << high_word << std::endl;
    //glog.log("word %li :: low %li (%li[%li]) high %li (%li[%li])\n",i,low_bit,low_word,low_bit%32,high_bit,high_word,high_bit%32);
    if (low_word == high_word) { //all the bits are in the same word
      unpacked[i] = (packed[low_word] >> (low_bit%32)) & 0x3FFF;
    } else { //some of the bits are in the next word
      size_t high_off = high_word*32-low_bit;
      //glog.log("pre_mask 0x%X post_mask 0x%X\n", (0x3FFF >> (14-high_off)), ((0x3FFF << high_off) & 0x3FFF) );
      unpacked[i] = (packed[low_word] >> (low_bit%32)) & (0x3FFF >> (14-high_off));
      unpacked[i] |= (packed[high_word] << high_off) & ((0x3FFF << high_off) & 0x3FFF);
    }
  }
}
DEFINE_ART_MODULE(IcebergFELIXBufferDecoderMarch2021)
