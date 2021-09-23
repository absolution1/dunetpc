// IcebergDataInterfaceFELIXBufferMarch2021_tool.cc

#include "IcebergDataInterfaceFELIXBufferMarch2021.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <set>
#include "lardataobj/RawData/raw.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Services/ChannelMap/IcebergChannelMapService.h"

IcebergDataInterfaceFELIXBufferMarch2021::IcebergDataInterfaceFELIXBufferMarch2021(fhicl::ParameterSet const& p)
{
  fInputFiles = p.get<std::vector<std::string>>("InputFiles");
  fNSamples = p.get<size_t>("NSamples",2000);
  fCompressHuffman = p.get<bool>("CompressHuffman",false);
  fDesiredStartTimestamp = p.get<ULong64_t>("StartTimestamp",0);

  fInputFilePointers.clear();
  for (size_t ifile=0; ifile<fInputFiles.size(); ++ifile)
    {
      fInputFilePointers.push_back(fopen(fInputFiles.at(ifile).data(),"r"));
    }
  fFirstRead = true;
}

// return all data for the event.  inputLabel is ignored for this as the input files are
// separately specified and data in them are unlabeled.

int IcebergDataInterfaceFELIXBufferMarch2021::retrieveData(art::Event &e, 
                                                           std::string inputLabel, 
                                                           std::vector<raw::RawDigit> &raw_digits, 
                                                           std::vector<raw::RDTimeStamp> &rd_timestamps,
                                                           std::vector<raw::RDStatus> &rdstatuses)
{

  art::ServiceHandle<dune::IcebergChannelMapService> channelMap;

  uint32_t framebuf[117];
  uint16_t databuf[128];
  uint64_t timestampstart=0;
  uint64_t timestamp=0;

  raw_digits.clear();
  rd_timestamps.clear();
  rdstatuses.clear();

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
              if (fInputFilePointers.at(ifile) == 0) return 0;
              fread(framebuf, sizeof(uint32_t), 117, fInputFilePointers.at(ifile));
              if (feof(fInputFilePointers.at(ifile)))
                {
                  // close all the input files and return nothing if we hit eof here
                  for (size_t jfile = 0; jfile < nfiles; ++jfile)
                    {
                      fclose(fInputFilePointers.at(jfile));
                      fInputFilePointers.at(jfile) = 0;
                    }
                  return 0;
                }
              timestamp= framebuf[3];
              timestamp <<= 32;
              timestamp += framebuf[2];
            }
          while (timestamp+16 < fDesiredStartTimestamp);
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
      if (fInputFilePointers.at(ifile) == 0) return 0;
      fread(fbcache.at(ifile).data(), sizeof(uint32_t), 117, fInputFilePointers.at(ifile));
      if (feof(fInputFilePointers.at(ifile)))
        {
          //  close all the input files and return nothing if we hit eof here
          for (size_t jfile = 0; jfile < nfiles; ++jfile)
            {
              fclose(fInputFilePointers.at(jfile));
              fInputFilePointers.at(jfile) = 0;
            }
          return 0;
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
          if (fInputFilePointers.at(ifile) == 0) return 0;
          fread(fbcache.at(ifile).data(), sizeof(uint32_t), 117, fInputFilePointers.at(ifile));
          if (feof(fInputFilePointers.at(ifile)))
            {
              // close all the input files and return nohting if we hit eof here
              for (size_t jfile = 0; jfile < nfiles; ++jfile)
                {
                  fclose(fInputFilePointers.at(jfile));
                  fInputFilePointers.at(jfile) = 0;
                }
              return 0;
            }
          timestamp= fbcache.at(ifile)[3];
          timestamp <<= 32;
          timestamp += fbcache.at(ifile)[2];
          tscache.at(ifile) = timestamp;
        }
    }  

  // actually read in the data now -- we already have the first tick read.

  for (size_t ifile=0; ifile < nfiles; ++ ifile)
    {
      if (fInputFilePointers.at(ifile) == 0) break;
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
                  // close all the input files and return what we have.
                  for (size_t jfile = 0; jfile < nfiles; ++jfile)
                    {
                      fclose(fInputFilePointers.at(jfile));
                      fInputFilePointers.at(jfile) = 0;
                    }
                  break;
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

          timestamp= framebuf[3];
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

          size_t uncompressed_nticks = adcvv.at(0).size();  
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

        }
    }
  // default all good status
  unsigned int statword = 0;
  rdstatuses.emplace_back(false,false,statword);
  //std::cout << "decoder felix tool nrawdigits: " << raw_digits.size() << std::endl;
  return 0;
}

// get data for specified APAs.  There's just one APA in ICEBERG so get all the data

int IcebergDataInterfaceFELIXBufferMarch2021::retrieveDataForSpecifiedAPAs(art::Event &evt, 
                                                                           std::vector<raw::RawDigit> &raw_digits, 
                                                                           std::vector<raw::RDTimeStamp> &rd_timestamps,
                                                                           std::vector<raw::RDStatus> &rdstatuses, 
                                                                           std::vector<int> &apalist)
{
  // ignore the APA list and also the input label
  return retrieveData(evt," ",raw_digits,rd_timestamps,rdstatuses);
}

// get data for a specific label, but only return those raw digits that correspond to APA's on the list
// again, just get all the data for this event

int IcebergDataInterfaceFELIXBufferMarch2021::retrieveDataAPAListWithLabels(art::Event &evt, 
                                                                            std::string inputLabel, 
                                                                            std::vector<raw::RawDigit> &raw_digits, 
                                                                            std::vector<raw::RDTimeStamp> &rd_timestamps,
                                                                            std::vector<raw::RDStatus> &rdstatuses, 
                                                                            std::vector<int> &apalist)
{
  return retrieveData(evt,inputLabel,raw_digits,rd_timestamps,rdstatuses);
}




// compute median and sigma.  

void IcebergDataInterfaceFELIXBufferMarch2021::computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, float &median, float &sigma)
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

void IcebergDataInterfaceFELIXBufferMarch2021::unpack14(const uint32_t *packed, uint16_t *unpacked) {
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

DEFINE_ART_CLASS_TOOL(IcebergDataInterfaceFELIXBufferMarch2021)
