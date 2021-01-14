// AdcSimulator.h

#ifndef AdcSimulator_H
#define AdcSimulator_H

// Interface for a tool that simulates the response of an ADC,
// i.e. converts an input voltage to an ADC count.
//
// The method count(...) returns the ADC count (i.e. bin) for a given
// input voltage. The caller also has the option of providing the 
// channel number and tick to allow the repective possibilities of
// channel and time dependence.
//
// Longer time dependence, e.g. event-to-event or day-to-day is not
// (yet) directly supported. We assume that concrete subclasses are
// otherwise made aware of the current event or time if there is to
// be such dependence.

class AdcSimulator {

public:

  using Channel = unsigned int;
  using Count = unsigned short;
  using Tick = unsigned int;

  virtual ~AdcSimulator() =default;

  virtual Count count(double vin, Channel chan =0, Tick tick =0) const =0;

};

#endif
