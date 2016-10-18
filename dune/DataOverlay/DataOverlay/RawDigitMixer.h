/**
 * \file RawDigitMixer.h
 *
 * \ingroup DataOverlay
 * 
 * \brief Mixer function for putting together two raw digit collections
 *
 * @author wketchum
 */

/** \addtogroup DataOverlay

    @{*/
#ifndef OVERLAY_DATAOVERLAY_RAWDIGITMIXER_H
#define OVERLAY_DATAOVERLAY_RAWDIGITMIXER_H

#include <vector>
#include <string>
#include <unordered_map>

#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "RawDigitAdder_35t.h"


/**
   \class RawDigitMixer
   Add two raw digit collections together.
   
*/

namespace mix {
  class RawDigitMixer;
}

class mix::RawDigitMixer{

public:

  /// Default constructor
  RawDigitMixer(bool p=false):
  _printWarnings(p){};

  void DeclareData(std::vector<raw::RawDigit> const& dataVector);
  void Mix(std::vector<raw::RawDigit> const& mcVector,
	   std::unordered_map<raw::ChannelID_t,float> const& map);

  void FillRawDigitOutput(std::vector<raw::RawDigit> & output);

  //void SetSaturationPoint(short x)
  //{ fRDAdderAlg.SetSaturationPoint(x); }
  
  /// Default destructor
  virtual ~RawDigitMixer(){};
  
  void SetTickOffsetForData(size_t o) { _offset=o; }

  void SetStuckBitRetentionMethod(bool b) { _stuckRetention = b; }

 private:

  bool _stuckRetention;
  bool _printWarnings;
  size_t _offset;

  //this is just to make storing this info easier, for RawDigit construction later
  struct RD_Info{
    std::vector<short> waveform;
    raw::ChannelID_t   channel;
    float              ped;
    float              sigma;
  };
  std::vector< RD_Info > fOutputWaveforms;

  std::unordered_map<raw::ChannelID_t,size_t> fChannelIndexMap;

  RawDigitAdder_35t fRDAdderAlg;

};

#endif
/** @} */ // end of doxygen group 

