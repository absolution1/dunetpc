#ifndef OVERLAY_DATAOVERLAY_RAWDIGITADDER_35T_CXX
#define OVERLAY_DATAOVERLAY_RAWDIGITADDER_35T_CXX

#include "RawDigitAdder_35t.h"
#include <limits>
#include <stdexcept>
#include <cmath>

mix::RawDigitAdder_35t::RawDigitAdder_35t(bool t):
  RawDigitAdder(t),
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
  d_out = (short)(std::round((float)d1 * _scale1)) + (short)(std::round((float)d2 * _scale2));
  FixOverflow(d_out);
}

#endif
