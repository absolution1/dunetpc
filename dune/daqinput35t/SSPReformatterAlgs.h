// -*- mode: c++; c-basic-offset: 2; -*-
#ifndef SSPReformatterAlgs_h
#define SSPReformatterAlgs_h

//artdaq includes
#include "artdaq-core/Data/Fragment.hh"


//larsoft includes
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "fhiclcpp/ParameterSet.h"

// lbne-raw-data includes
#include "lbne-raw-data/Overlays/SSPFragment.hh"
#include "lbne-raw-data/Overlays/anlTypes.hh"

//lbnecode/daqinput35t includes
#include "utilities/UnpackFragment.h"
#include <vector>
#include <map>

// Unpack the given artdaq::Fragment objects, and create a vector of raw::OpDetWaveform objects. The
// Fragments are expected to be carrying Optical Detector data; this is not
// checked.

namespace DAQToOffline {

  class SSPReformatterAlgs
  {
    
  public:
    SSPReformatterAlgs(fhicl::ParameterSet const & pset);
    
    /**
     * Construct OpDetWaveform objects from triggers with full waveforms and
     * construct OpHit objects from triggers which only have headers.
     */
    void SSPFragmentToWaveformsAndHits(artdaq::Fragments const& rawFragments,
                                       std::vector<raw::OpDetWaveform> &opDetWaveformVector,
                                       std::vector<recob::OpHit> &opHitVector);
    
    /// Construct a waveform from each trigger
    std::vector<raw::OpDetWaveform> SSPFragmentToOpDetWaveform(artdaq::Fragments const& raw);
    std::vector<recob::OpHit> SSPHeaderToOpHit(artdaq::Fragments const& raw);
    
    /// Print out header information
    void PrintHeaderInfo(const SSPDAQ::EventHeader *daqHeader);

    /// Load the milislice
    unsigned int CheckAndGetNTriggers(const artdaq::Fragment& frag, const lbne::SSPFragment sspf);
    
    /// Load the header
    const SSPDAQ::EventHeader* GetHeaderAndAdvance(const unsigned int* &dataPointer);

    // Extract data from the header
    uint32_t GetPeakSum(const SSPDAQ::EventHeader* daqHeader);
    unsigned short GetOpChannel(const SSPDAQ::EventHeader* daqHeader);
    unsigned long GetGlobalFirstSample(const SSPDAQ::EventHeader* daqHeader);
    unsigned long GetInternalFirstSample(const SSPDAQ::EventHeader *daqHeader);
    unsigned long GetBaselineSum(const SSPDAQ::EventHeader *daqHeader);
    unsigned long GetIntegratedSum(const SSPDAQ::EventHeader *daqHeader);
    unsigned int  GetPeakTime(const SSPDAQ::EventHeader *daqHeader);
    unsigned int  GetWaveformLength(const SSPDAQ::EventHeader *daqHeader);

    /// Return the NOvAClockFrequency
    double ClockFrequency() { return NOvAClockFrequency; }

  private:
    /// Construct a waveform from the adc vector, advance the data pointer when done
    raw::OpDetWaveform ConstructWaveformAndAdvance(const SSPDAQ::EventHeader* daqHeader,
                                                   const unsigned int* &dataPointer);
    /// Construct an OpHit object from the daqHeader
    recob::OpHit       ConstructOpHit(const SSPDAQ::EventHeader* daqHeader);

    
    double NOvAClockFrequency;
    std::map<int,int> theChannelMap;
    
    int m1;
    int i1;
    int i2;
    double SPESize;
    

    
    void BuildOpDetChannelMap(std::string fChannelMapFile);
    
  };
  
}
#endif
