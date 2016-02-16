// -*- mode: c++; c-basic-offset: 2; -*-
// art includes

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib/search_path.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include <fstream>

// lbnecode/daqinput35t includes

#include "SSPReformatterAlgs.h"

// larsoft includes

#include "larcore/Geometry/Geometry.h"

// lbne-raw-data includes
#include "lbne-raw-data/Overlays/SSPFragment.hh"
#include "lbne-raw-data/Overlays/anlTypes.hh"



//////////////////////////////////////////////////////////////////////////////////////////////////

DAQToOffline::SSPReformatterAlgs::SSPReformatterAlgs(fhicl::ParameterSet const & pset) :
                                                     NOvAClockFrequency(pset.get<double>("NOvAClockFrequency")), // in MHz
                                                     m1(pset.get<int>("SSPm1")),
                                                     i1(pset.get<int>("SSPi1")),
                                                     i2(pset.get<int>("SSPi2")),
                                                     SPESize(pset.get<float>("SPESize"))

{
  BuildOpDetChannelMap(pset.get<std::string>("OpDetChannelMapFile"));
}


//////////////////////////////////////////////////////////////////////////////////////////////////


void DAQToOffline::SSPReformatterAlgs::SSPFragmentToWaveformsAndHits(artdaq::Fragments const& rawFragments,
                                                                     std::vector<raw::OpDetWaveform> &opDetWaveformVector,
                                                                     std::vector<recob::OpHit> &opHitVector)
{
  unsigned int numFragments = rawFragments.size();
  
  for (size_t idx = 0; idx < numFragments; ++idx) {
    const auto& frag(rawFragments[idx]);
    lbne::SSPFragment sspf(frag);
    
    unsigned int nTriggers = CheckAndGetNTriggers(frag, sspf);
    
    const unsigned int* dataPointer = sspf.dataBegin();
    
    
    for (unsigned int triggersProcessed = 0;
         (nTriggers==0 || triggersProcessed < nTriggers) && dataPointer < sspf.dataEnd();
         ++triggersProcessed) {
      
      try {
        auto daqHeader = GetHeaderAndAdvance(dataPointer);
        if ( GetWaveformLength(daqHeader) > 0) {
          // We have a waveform
          opDetWaveformVector.emplace_back( ConstructWaveformAndAdvance(daqHeader, dataPointer) );
          
        }
        else {
          // We only have a header
          opHitVector.emplace_back( ConstructOpHit(daqHeader) );
        }
        
      }
      catch (cet::exception e) {
        continue;
      }
    } // End of loop over triggers
  } // End of loop over fragments (rawFragments)

  
}


//////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<raw::OpDetWaveform> DAQToOffline::SSPReformatterAlgs::SSPFragmentToOpDetWaveform(artdaq::Fragments const& rawFragments)
{
  std::vector<raw::OpDetWaveform> opDetWaveformVector;

  unsigned int numFragments = rawFragments.size();

  for (size_t idx = 0; idx < numFragments; ++idx) {
    const auto& frag(rawFragments[idx]);
    lbne::SSPFragment sspf(frag);

    unsigned int nTriggers = CheckAndGetNTriggers(frag, sspf);
      
    const unsigned int* dataPointer = sspf.dataBegin();

        
    for (unsigned int triggersProcessed = 0;
         (nTriggers==0 || triggersProcessed < nTriggers) && dataPointer < sspf.dataEnd();
         ++triggersProcessed) {

      try {
        auto daqHeader = GetHeaderAndAdvance(dataPointer);
        opDetWaveformVector.emplace_back( ConstructWaveformAndAdvance(daqHeader, dataPointer) );
      }
      catch (cet::exception e) {
        continue;
      }
      
    } // End of loop over triggers
  } // End of loop over fragments (rawFragments)

  return opDetWaveformVector;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<recob::OpHit> DAQToOffline::SSPReformatterAlgs::SSPHeaderToOpHit(artdaq::Fragments const& rawFragments)
{
  std::vector<recob::OpHit> opHitVector;
  
  unsigned int numFragments = rawFragments.size();
  
  
  for (size_t idx = 0; idx < numFragments; ++idx) {
    const auto& frag(rawFragments[idx]);
    lbne::SSPFragment sspf(frag);
    
    unsigned int nTriggers = CheckAndGetNTriggers(frag, sspf);
    const unsigned int* dataPointer = sspf.dataBegin();
    
    
    for (unsigned int triggersProcessed = 0;
         (nTriggers==0 || triggersProcessed < nTriggers) && dataPointer < sspf.dataEnd();
         ++triggersProcessed) {
      
      
      try {
        auto daqHeader = GetHeaderAndAdvance(dataPointer);
        opHitVector.emplace_back( ConstructOpHit(daqHeader) );
        
        // Advance the dataPointer to the next header
        unsigned int nADC = GetWaveformLength(daqHeader);
        dataPointer+=nADC/2;
      }
      catch (cet::exception e) {
        continue;
      }
    } // End of loop over triggers
  } // End of loop over fragments (rawFragments)
  
  return opHitVector;
}



//////////////////////////////////////////////////////////////////////////////////////////////////

void DAQToOffline::SSPReformatterAlgs::PrintHeaderInfo(const SSPDAQ::EventHeader *daqHeader)
{
  auto FirstSample = GetGlobalFirstSample(daqHeader);
  auto InternalSample = GetInternalFirstSample(daqHeader);
  auto TimeStamp = ((double)FirstSample)/NOvAClockFrequency;
  auto InternalTimeStamp = ((double)InternalSample)/NOvAClockFrequency;
  
  mf::LogDebug("DAQToOffline")
  << "Header:                             " << daqHeader->header   << "\n"
  << "Length:                             " << daqHeader->length   << "\n"
  << "Trigger type:                       " << ((daqHeader->group1 & 0xFF00) >> 8) << "\n"
  << "Status flags:                       " << ((daqHeader->group1 & 0x00F0) >> 4) << "\n"
  << "Header type:                        " << ((daqHeader->group1 & 0x000F) >> 0) << "\n"
  << "Trigger ID:                         " << daqHeader->triggerID << "\n"
  << "Module ID:                          " << ((daqHeader->group2 & 0xFFF0) >> 4) << "\n"
  << "Channel ID:                         " << ((daqHeader->group2 & 0x000F) >> 0) << "\n"
  << "OpChannel:                          " << GetOpChannel(daqHeader) << "\n"
  << "External (NOvA) timestamp:          " << FirstSample << " ticks" << "\n"
  << "                                    " << TimeStamp << " microseconds" << "\n"
  // these first_ ouptuts are a little ill-defined anyway since they are not
  // reset by run but rather with the job
  //<< "Since first sample this run:        " << FirstSample-first_FirstSample << " ticks" << "\n"
  //<< "                                    " << TimeStamp-first_TimeStamp << " microseconds" << "\n"
  << "Peak sum:                           " << GetPeakSum(daqHeader) << "\n"
  << "Peak time:                          " << GetPeakTime(daqHeader) << "\n"
  
  << "Baseline Sum (Prerise):             " << GetBaselineSum(daqHeader) << "\n"
  << "Integrated sum:                     " << GetIntegratedSum(daqHeader) << "\n"
  << "Baseline:                           " << daqHeader->baseline << "\n"
  << "CFD Timestamp interpolation points: " << daqHeader->cfdPoint[0] << " " << daqHeader->cfdPoint[1]
  << " " << daqHeader->cfdPoint[2] << " " << daqHeader->cfdPoint[3] << "\n"
  << "Internal interpolation point:       " << daqHeader->intTimestamp[0] << "\n"
  << "Internal timestamp:                 " << InternalSample << " ticks\n"
  << "                                    " << InternalTimeStamp << " microseconds" << "\n"
  //<< "Relative internal timestamp:        " << InternalSample-first_InternalSample << " ticks" << "\n"
  //<< "                                    " << InternalTimeStamp-first_InternalTimeStamp << " microseconds"
  << ""  ;
  
}



//////////////////////////////////////////////////////////////////////////////////////////////////


unsigned int DAQToOffline::SSPReformatterAlgs::CheckAndGetNTriggers(const artdaq::Fragment& frag, const lbne::SSPFragment sspf)
{
    mf::LogDebug("DAQToOffline") << "\n"
                                 << "SSP fragment "     << frag.fragmentID() 
                                 << " has total size: " << sspf.hdr_event_size()
                                 << " and run number: " << sspf.hdr_run_number()
                                 << " with " << sspf.total_adc_values() << " total ADC values"
                                 << "\n"
                                 << "\n";

    const SSPDAQ::MillisliceHeader* meta=0;
    //get the information from the header
    if(frag.hasMetadata())
    {
        meta = &(frag.metadata<lbne::SSPFragment::Metadata>()->sliceHeader);
            
        mf::LogInfo("DAQToOffline")
          << "===Slice metadata====" << "\n"
          << "  Start time         " << meta->startTime << "\n"
          << "  End time           " << meta->endTime << "\n"
          << "  Packet length      " << meta->length << "\n"
          << "  Number of triggers " << meta->nTriggers << "\n"
          << "=====================";
    }
    else
    {
        mf::LogWarning("DAQToOffline") << "SSP fragment has no metadata associated with it.";
    }

    // No metadata, no trigger count
    if (meta == 0) return 0;

    return meta->nTriggers;
}


//////////////////////////////////////////////////////////////////////////////////////////////////


void DAQToOffline::SSPReformatterAlgs::BuildOpDetChannelMap(std::string fChannelMapFile)
{
    theChannelMap.clear();

    int onlineChannel;
    int offlineChannel;
    
    std::string fullname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fChannelMapFile, fullname);
    
    if ( fullname.empty() ) {
      mf::LogWarning("DAQToOffline") << "Input spectrum file "
                                        << fChannelMapFile
                                        << " not found in FW_SEARCH_PATH, using debugging map!\n";
    } 
    else {
      mf::LogVerbatim("DAQToOffline") << "Build Online->Offline channel Map from " << fullname;
      std::ifstream infile(fullname);
      while (infile.good()) {
        infile >> onlineChannel >> offlineChannel;
        theChannelMap[onlineChannel] = offlineChannel;
        mf::LogVerbatim("DAQToOffline") << "   " << onlineChannel << " -> " << offlineChannel;
      }
    }
    
}


//////////////////////////////////////////////////////////////////////////////////////////////////

const SSPDAQ::EventHeader* DAQToOffline::SSPReformatterAlgs::GetHeaderAndAdvance(const unsigned int* &dataPointer)
{
  // Load the event header, advance the pointer
  const SSPDAQ::EventHeader* daqHeader=reinterpret_cast<const SSPDAQ::EventHeader*>(dataPointer);
  dataPointer += sizeof(SSPDAQ::EventHeader)/sizeof(unsigned int);
  return daqHeader;
}



//////////////////////////////////////////////////////////////////////////////////////////////////


uint32_t DAQToOffline::SSPReformatterAlgs::GetPeakSum(const SSPDAQ::EventHeader* daqHeader)
{
  uint32_t peaksum = ((daqHeader->group3 & 0x00FF) >> 16) + daqHeader->peakSumLow;
  if(peaksum & 0x00800000) {
    peaksum |= 0xFF000000;
  }
  return peaksum;
}



//////////////////////////////////////////////////////////////////////////////////////////////////


unsigned short DAQToOffline::SSPReformatterAlgs::GetOpChannel(const SSPDAQ::EventHeader* daqHeader)
{
  unsigned short OpChannel = -1;
  
  // Extract values we need for the data product
  if ( theChannelMap.size() == 0) {
    // No channel map, default to debugging map
    int HardwareChannel = ((daqHeader->group2 & 0x000F) >> 0); // Channel Number
    int SSPNumber       = ((daqHeader->group2 & 0x00F0) >> 4); // Module Number
    OpChannel = 100*SSPNumber + HardwareChannel;
  }
  else if ( theChannelMap.find(daqHeader->group2) != theChannelMap.end() ) {
    OpChannel = theChannelMap[daqHeader->group2];
  }
  else {
    int HardwareChannel = ((daqHeader->group2 & 0x000F) >> 0); // Channel Number
    int SSPNumber       = ((daqHeader->group2 & 0x00F0) >> 4); // Module Number
    mf::LogWarning("DAQToOffline") << "SSP " << SSPNumber << " Channel " << HardwareChannel
                                             << "(" << daqHeader->group2 << ") "
                                             << " not in the map (OK for uninstrumented channels), skipping." << std::endl;
    throw cet::exception( "SSP Channel Invalid" );
  }

  return OpChannel;
}


//////////////////////////////////////////////////////////////////////////////////////////////////


unsigned long DAQToOffline::SSPReformatterAlgs::GetGlobalFirstSample(const SSPDAQ::EventHeader* daqHeader)
{
  return (   ( (unsigned long)daqHeader->timestamp[3] << 48 )
           + ( (unsigned long)daqHeader->timestamp[2] << 32 )
           + ( (unsigned long)daqHeader->timestamp[1] << 16 )
           + ( (unsigned long)daqHeader->timestamp[0] ) );
}


//////////////////////////////////////////////////////////////////////////////////////////////////


unsigned long DAQToOffline::SSPReformatterAlgs::GetInternalFirstSample(const SSPDAQ::EventHeader *daqHeader)
{
  return (   ((uint64_t)((uint64_t)daqHeader->intTimestamp[3] << 32))
           + ((uint64_t)((uint64_t)daqHeader->intTimestamp[2]) << 16)
           + ((uint64_t)((uint64_t)daqHeader->intTimestamp[1])) );
}


//////////////////////////////////////////////////////////////////////////////////////////////////

unsigned long DAQToOffline::SSPReformatterAlgs::GetBaselineSum(const SSPDAQ::EventHeader *daqHeader)
{
  return ((daqHeader->group4 & 0x00FF) << 16) + daqHeader->preriseLow;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////

unsigned long DAQToOffline::SSPReformatterAlgs::GetIntegratedSum(const SSPDAQ::EventHeader *daqHeader)
{
  return ((unsigned int)(daqHeader->intSumHigh) << 8) + (((unsigned int)(daqHeader->group4) & 0xFF00) >> 8);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int DAQToOffline::SSPReformatterAlgs::GetPeakTime(const SSPDAQ::EventHeader *daqHeader)
{
  return (daqHeader->group3 & 0xFF00) >> 8 ;
}


//////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int DAQToOffline::SSPReformatterAlgs::GetWaveformLength(const SSPDAQ::EventHeader *daqHeader)
{
  unsigned int nADC=(daqHeader->length-sizeof(SSPDAQ::EventHeader)/sizeof(unsigned int))*2;
  return nADC;
}


//////////////////////////////////////////////////////////////////////////////////////////////////

raw::OpDetWaveform DAQToOffline::SSPReformatterAlgs::ConstructWaveformAndAdvance(const SSPDAQ::EventHeader* daqHeader,
                                                                                 const unsigned int* &dataPointer)
{
  
  // Get basic information from the header
  unsigned short     OpChannel   = GetOpChannel(daqHeader);                  ///< Derived Optical channel
  unsigned long      FirstSample = GetGlobalFirstSample(daqHeader);          ///< first sample time in ticks
  double             TimeStamp   = ((double)FirstSample)/NOvAClockFrequency; ///< first sample time in microseconds

  // Print some debugging info
  PrintHeaderInfo(daqHeader);
  
  // Get ADC Count, create pointer to adcs
  unsigned int          nADC       = GetWaveformLength(daqHeader);
  const unsigned short* adcPointer = reinterpret_cast<const unsigned short*>(dataPointer);
  
  // Initialize the waveform
  raw::OpDetWaveform Waveform(TimeStamp, OpChannel, nADC);
  
  // Build up the waveform
  for(size_t idata = 0; idata < nADC; idata++) {
    const unsigned short* adc = adcPointer + idata;
    Waveform.push_back(*adc);
  }
  
  // Advance the dataPointer to the next header
  dataPointer+=nADC/2;
  
  return Waveform;
}


//////////////////////////////////////////////////////////////////////////////////////////////////

recob::OpHit DAQToOffline::SSPReformatterAlgs::ConstructOpHit(const SSPDAQ::EventHeader* daqHeader)
{
  // Get basic information from the header
  unsigned short     OpChannel   = GetOpChannel(daqHeader);                  ///< Derived Optical channel
  unsigned long      FirstSample = GetGlobalFirstSample(daqHeader);;         ///< first sample time in ticks
  double             TimeStamp   = ((double)FirstSample)/NOvAClockFrequency; ///< first sample time in microseconds
  
  auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
  
  double peakTime = ((double)GetPeakTime(daqHeader)) * ts->OpticalClock().TickPeriod(); // microseconds
  double width = ((double)i1) * ts->OpticalClock().TickPeriod(); // microseconds
  
  double pedestal = ((double)GetBaselineSum(daqHeader)) / ((double)i1);
  double area = ((double)GetIntegratedSum(daqHeader))  - pedestal * i2;
  double peak = ((double)GetPeakSum(daqHeader)) / ((double)m1) - pedestal;
  
  
  recob::OpHit ophit(OpChannel,
                     TimeStamp+peakTime,  // Relative Time
                     TimeStamp+peakTime,  // Absolute time
                     0,          // Frame, not used by DUNE
                     width,
                     area,
                     peak,
                     area / SPESize, // PE
                     0.);
  
  return ophit;
}

