// TH1Manipulator.h
//
// David Adams
// November 2015
//
// Root macros to add top and/or right axis to the current pad.
//
// If histogram is specified, the axis attributes (tick size and
// # divisions) are taken from the histogram.
//
// After modifying plot range, repeat call to fix top/right axis.

#ifndef TH1Manipulator_H
#define TH1Manipulator_H

class TH1;

class TH1Manipulator {

public:

  // Draw top x-axis and right y-axis.
  // If ph != 0, axis attributes are taken from the histogram.
  static int addaxis(TH1* ph =0);

  // Draw top x-axis with range taken from a histogram.
  // If ph != 0, axis attributes are taken from the histogram.
  static int addaxistop(TH1* ph =0);

  // Draw top x-axis with specified attributes.
  static int addaxistop(double ticksize, int ndiv);

  // Draw right y-axis with range taken from a histogram.
  // If ph != 0, axis attributes are taken from the histogram.
  static int addaxisright(TH1* ph =0);

  // Draw right y-axis with specified attributes.
  static int addaxisright(double ticksize, int ndiv);

  // Fix the BG color for 2D histos to be the same as the lowest color.
  // Otherwise underflows look like zeros.
  static int fixFrameFillColor();

};

#endif
