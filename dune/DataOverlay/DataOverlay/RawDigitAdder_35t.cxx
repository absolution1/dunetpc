#ifndef OVERLAY_DATAOVERLAY_RAWDIGITADDER_35T_CXX
#define OVERLAY_DATAOVERLAY_RAWDIGITADDER_35T_CXX

#include "RawDigitAdder_35t.h"
#include <limits>
#include <stdexcept>
#include <cmath>
#include "dune/DetSim/Utility/AdcCodeHelper.h"

mix::RawDigitAdder_35t::RawDigitAdder_35t(bool t):
  RawDigitAdder(t),
  _forceStuckBitRetention(false),
  _scale1(1),
  _scale2(1)
{}

void mix::RawDigitAdder_35t::SetScaleInput(float f, float& _scale)
{
  if(f<0){
    if(_throw)
      throw std::runtime_error("Error in RawDigitAdder_35t::SetScaleInput : scale < 0");
    return;
  }
  _scale = f;
}

void mix::RawDigitAdder_35t::AddRawDigit(short const& d1, short const& d2, short& d_out)
{
  // d1 ==> MC
  // d2 ==> Real Data
  AdcCodeHelper ach;
  if (ach.hasStickyBits(d2) && _forceStuckBitRetention) d_out = (short)(std::round((float)d2 * _scale2)); 
  else d_out = (short)(std::round((float)d1 * _scale1)) + (short)(std::round((float)d2 * _scale2));
  FixOverflow(d_out);
}

#endif
