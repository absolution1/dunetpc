// TPadManipulator.h
//
// David Adams
// November 2015
//
// Root macros to manipulate a root drawing pad.
//

#ifndef TPadManipulator_H
#define TPadManipulator_H

class TVirtualPad;
class TH1;

class TPadManipulator {

public:

  // Ctor from a pad. If null the current pad is used.
  TPadManipulator(TVirtualPad* ppad =nullptr);

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

  // Return the first histogram for this pad.
  TH1* hist() const { return m_ph; }

  // Update the coordinates and histogram for this pad.
  int update();

  // Draw top and right axis using the first histogram on the pad.
  int addaxis();

  // Draw top x-axis with range taken from the first histogram on the pad.
  int addaxistop();

  // Draw top x-axis with specified attributes.
  int addaxistop(double ticksize, int ndiv);

  // Draw right y-axis with range taken from the first histogram on the pad.
  // If ph != 0, axis attributes are taken from the histogram.
  int addaxisright();

  // Draw right y-axis with specified attributes.
  int addaxisright(double ticksize, int ndiv);

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

};

#endif
