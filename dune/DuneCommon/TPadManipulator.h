// TPadManipulator.h
//
// David Adams
// Octobe 2017
//
// Root macro to manipulate a root drawing pad.
//
// The first histogram or graph in the pad is cloned and the pad is updated
// using the clone. Caller then has the options to
//  change X and Y range
//  add top and bottom axis
//  display underflow and overflow (histograms only)
//  add vertical lines at regular intervals
//  display any of the functions (e.g. the fit) associated with the histogram or graph
//
// It is also posssible to and multiple subpads and draw into those instead of the
// top-level pad. There is a manipulator for each subpad.

#ifndef TPadManipulator_H
#define TPadManipulator_H

#include <vector>
#include <memory>
#include "TLine.h"

class TVirtualPad;
class TH1;
class TGraph;
class TAxis;

class TPadManipulator {

public:

  using Index = unsigned int;
  using TLinePtr = std::shared_ptr<TLine>;

  // Ctor from a pad. If null the current pad is used.
  TPadManipulator(TVirtualPad* ppad =nullptr);

  // Ctors that creates a new canvas.
  // The canvas size is wx x wy.
  // If nPadX > 0, calls split(nPadX, nPadY) if nPadY > 0
  //                  or split(nPadX, nPadX) if nPadY == 0
  TPadManipulator(Index wx, Index wy, Index nPadX =0, Index nPadY =0);

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

  // Return the top-level pad.
  // Use pad()->SetFillColor(...), etc. to change the appearance of the pad.
  // and pad()->SetFrameFillColor(...), etc. to change the appearance of the frame in the pad.
  TVirtualPad* pad() const { return m_ppad; }

  // Return the number of subpads.
  unsigned int npad() const { return m_subMans.size(); }

  // Return the manipulator for a subpad.
  // If there are no subpads and ipad=0,  the top manipulator is returned.
  TPadManipulator* man(unsigned int ipad =0) {
    if ( npad() == 0 && ipad == 0 ) return this;
    return ipad<npad() ? &m_subMans[ipad] : nullptr;
  }

  // Return the histogram or graph for this pad.
  TH1* hist() const { return m_ph.get(); }
  TGraph* graph() const { return m_pg.get(); }

  // Return the vertical mod lines associated with this pad.
  const std::vector<TLinePtr>& verticalModLines() const { return m_vmlLines; }

  // Add a subpad covering (x1, y1) to (x2, y2) in NDC units, i.e. 0 < x,y < 1.
  int addPad(double x1, double y1, double x2, double y2, int icol =-1);

  // Tile the top pad with nx x ny  or nx x nx subpads.
  int split(Index nx, Index ny);
  int split(Index nx);

  // Add an object (histogram or graph) to the pad.
  // For now there can only be one object/pad. If replace is true, the old object is
  // removed, otherwise this methd fails.
  int add(unsigned int ipad, TObject* pobj, std::string sopt, bool replace =false);
  int add(TObject* pobj, std::string sopt, bool replace =false);

  // Remove histograms and graphs from top and subpads.
  int clear();

  // Update the coordinates and histogram for this pad.
  int update();

  // Get the axes of the histogram or graph.
  TAxis* getXaxis();
  TAxis* getYaxis();

  // Set the histogram range
  int setRangeX(double x1, double x2);
  int setRangeY(double y1, double y2);
  int setRanges(double x1, double x2, double y1, double y2);

  // Add top and right axis using attributes of the first histogram on the pad.
  int addAxis();

  // Add top x-axis with attributes taken from the first histogram on the pad.
  int addAxisTop();

  // Add top x-axis with specified attributes.
  int addAxisTop(double ticksize, int ndiv);

  // Draw the top axis.
  int drawAxisTop();

  // Add right y-axis with attributes from the first histogram on the pad.
  int addAxisRight();

  // Add right y-axis with specified attributes.
  int addAxisRight(double ticksize, int ndiv);

  // Draw right y-axis.
  int drawAxisRight();

  // Add underflow bin to the plot (false to remove it).
  int showUnderflow(bool show =true);
  int showOverflow(bool show =true);

  // Return the under/overflow histogram.
  TH1* flowHistogram() { return m_flowHist.get(); }

  // Add vertical modulus lines.
  // I.e at x = xoff, xoff+/-xmod, xoff+/-2*xmod, ...
  // The lines are draw from the bottom to lenfrac*height
  int addVerticalModLines(double xmod, double xoff =0.0, double lenfrac =1.0);

  // Add histogram function ifun to the pad.
  int addHistFun(unsigned int ifun =0);

  // Draw the hist functions.
  int drawHistFuns();

  // Draw the current mod lines.
  int drawVerticalModLines();

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
  std::shared_ptr<TH1> m_ph;
  std::shared_ptr<TGraph> m_pg;
  std::shared_ptr<TH1> m_flowHist;
  bool m_showUnderflow;
  bool m_showOverflow;
  bool m_top;
  double m_topTicksize;
  double m_topNdiv;
  bool m_right;
  double m_rightTicksize;
  double m_rightNdiv;
  std::vector<unsigned int> m_histFuns;
  double m_vmlXmod;
  double m_vmlXoff;
  double m_vmlXLength;
  std::vector<std::shared_ptr<TLine>> m_vmlLines;
  std::vector<TPadManipulator> m_subMans;

};

#endif
