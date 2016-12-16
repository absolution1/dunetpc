/**
 * \file RawDigitAdder_35t.h
 *
 * \ingroup DataOverlay
 * 
 * \brief Defintion for a class to add two vectors together,
 *        and give an "added" waveform.
 *
 *
 * @author wketchum
 * @author mthiesse
 */

/** \addtogroup DataOverlay

    @{*/
#ifndef OVERLAY_DATAOVERLAY_RAWDIGITADDER_35T_H
#define OVERLAY_DATAOVERLAY_RAWDIGITADDER_35T_H

#include <vector>
#include <string>
#include "RawDigitAdder.h"

/**
   \class RawDigitAdder_35t
   Add two vectors together. 
   Take stuck bits into account -- i.e. do not add when stuck
   Allows for a scale factor to be applied to inputs: defaults to 1.
   
*/
namespace mix {
  class RawDigitAdder_35t;
}

class mix::RawDigitAdder_35t : public mix::RawDigitAdder {

public:

  RawDigitAdder_35t(bool t=true);

  
  void SetScaleFirstInput(float f1)  { SetScaleInput(f1,_scale1); }
  void SetScaleSecondInput(float f2) { SetScaleInput(f2,_scale2); }

  void SetStuckBitRetentionMethod(bool s) { _forceStuckBitRetention = s; }

  void SetScaleInputs(float f1, float f2)
  { SetScaleFirstInput(f1); SetScaleSecondInput(f2); }
  
  std::string Name() { return "RawDigitAdder_35t"; }
  
 private:

  bool _forceStuckBitRetention;
  float _scale1,_scale2;
  void SetScaleInput(float f, float& _scale);
  void AddRawDigit( short const&, short const&, short&);

};

#endif
/** @} */ // end of doxygen group 

