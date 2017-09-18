// TPadManipulator.h
//
// David Adams
// November 2015
//
// Root macros to manipulate a root drawing pad.
//

#ifndef TPadManipulator_H
#define TPadManipulator_H

#include <vector>

class TVirtualPad;
class TH1;
class TLine;

class TPadManipulator {

public:

  // Ctor from a pad. If null the current pad is used.
  TPadManipulator(TVirtualPad* ppad =nullptr);

  // Dtor.
  // This removes the lines.
  ~TPadManipulator();

  // Return the axis coordinates in pad units (pixels).
  double xminPad() const { return m_xminPad; }
  double xmaxPad() const { return m_xmaxPad; }
  double yminPad() const { return m_yminPad; }
  double ymaxPad() const { return m_ymaxPad; }

  // Return the axis coordinates in drawn units.
  double xmin() const { return m_xmin; }
  double xmax() const { return m_xmax; }
  double ymin() const { return m_ymin; }
  double ymax() const { return m_ymax; }

  // Return the pad.
  TVirtualPad* pad() const { return m_ppad; }

  // Return the first histogram for this pad.
  TH1* hist() const { return m_ph; }

  // Return the lines associated with this pad.
  const std::vector<TLine*>& lines() const { return m_lines; }
  // Update the coordinates and histogram for this pad.
  int update();

  // Draw top and right axis using the first histogram on the pad.
  int addAxis();

  // Draw top x-axis with range taken from the first histogram on the pad.
  int addAxisTop();

  // Draw top x-axis with specified attributes.
  int addAxisTop(double ticksize, int ndiv);

  // Draw right y-axis with range taken from the first histogram on the pad.
  // If ph != 0, axis attributes are taken from the histogram.
  int addAxisRight();

  // Draw right y-axis with specified attributes.
  int addAxisRight(double ticksize, int ndiv);

  // Add vertical modulus lines.
  // I.e at x = xoff, xoff+/-xmod, xoff+/-2*xmod, ...
  int addVerticalModLines(double xmod, double xoff =0.0);

  // Fix the BG color for 2D histos to be the same as the lowest color.
  // Otherwise underflows have the color of zeros.
  int fixFrameFillColor();

private:

  TVirtualPad* m_ppad;
  double m_xminPad;
  double m_xmaxPad;
  double m_yminPad;
  double m_ymaxPad;
  double m_xmin;
  double m_xmax;
  double m_ymin;
  double m_ymax;
  TH1* m_ph;
  std::vector<TLine*> m_lines;

};

#endif
