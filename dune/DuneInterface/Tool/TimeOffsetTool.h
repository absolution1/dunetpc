// TimeOffsetTool.h

#ifndef TimeOffsetTool_H
#define TimeOffsetTool_H

// Interface for a tool providing access to a timing offset,
// e.g. for the ADC samples from a TPC.

#include<string>

class TimeOffsetTool {

public:

  using Index = unsigned int;
  using LongIndex = unsigned long;
  using Name = std::string;

  // Data on which the offset may depend.
  class Data {
  public:
    Index run =0;
    Index subrun =0;
    Index event =0;
    Index channel =0;
    LongIndex triggerClock =0;
  };

  // Class describing a time offset.
  //  value - Integral part of the offset (typically sufficient?)
  //    rem - value remainder if offset in non-integral
  //   unit - Unit for offset, e.g. tick, us, ns
  class Offset {
  public:
    int status =0;
    long value =0;
    double rem =0.0;
    Name unit;
    bool isValid() const { return status == 0; }
    const Offset& setStatus(int val) { status = val; return *this; }
  };

  virtual ~TimeOffsetTool() =default;

  virtual Offset offset(const Data& dat) const =0;

};

#endif
