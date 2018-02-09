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
#include "TLatex.h"
#include "TLine.h"

class TVirtualPad;
class TCanvas;
class TH1;
class TGraph;
class TAxis;
class TLegend;
class TFrame;

class TPadManipulator {

public:

  struct Bounds {
    double x1;
    double y1;
    double x2;
    double y2;
    Bounds() : x1(0.0), y1(0.0), x2(0.0), y2(0.0) { }
    Bounds(double ax1, double ay1, double ax2, double ay2)
    : x1(ax1), y1(ay1), x2(ax2), y2(ay2) { }
  };

  struct Bounds3 : public Bounds {
    double z1;
    double z2;
    Bounds3() : z1(0.0), z2(0.0) { }
    Bounds3(double ax1, double ay1, double ax2, double ay2, double az1, double az2)
    : Bounds(ax1, ay1, ax2, ay2), z1(az1), z2(az2) { }
  };

  using Index = unsigned int;
  using TLinePtr = std::shared_ptr<TLine>;
  using TObjPtr = std::shared_ptr<TObject>;
  using TObjVector = std::vector<TObjPtr>;
  using BoundsVector = std::vector<Bounds>;

  // Default ctor.
  // Creates an empty top-level object.
  TPadManipulator();

  // Ctor from a pad. Use ppad = gPad to grab the current pad.
  // Copies the objects from the pad.
  TPadManipulator(TVirtualPad* ppad);

  // Ctors that creates a new canvas.
  // The canvas size is wx x wy.
  // If nPadX > 0, calls split(nPadX, nPadY) if nPadY > 0
  //                  or split(nPadX, nPadX) if nPadY == 0
  TPadManipulator(Index wx, Index wy, Index nPadX =0, Index nPadY =0);

  // Copy and assignment.
  // Objects from the RHS are cloned so that mods in one pad don't affect the other.
  // These can be used to copy an existing pad into a subpad of another.
  // In assignment, the original objects belonging to this pad are deleted.
  TPadManipulator(const TPadManipulator& rhs);
  TPadManipulator& operator=(const TPadManipulator& rhs);

  // Dtor.
  // This removes the lines.
  ~TPadManipulator();

  // Return the axis coordinates in pad units (pixels).
  //double xminPad() const { return m_xminPad; }
  //double xmaxPad() const { return m_xmaxPad; }
  //double yminPad() const { return m_yminPad; }
  //double ymaxPad() const { return m_ymaxPad; }

  // Return the axis limits in drawn units.
  // These will be the user-specified range if valid or, for X and Y, the drawn pad.
  double xmin() const;
  double xmax() const;
  double ymin() const;
  double ymax() const;
  double zmin() const;
  double zmax() const;

  // Return the pad.
  // Use pad()->SetFillColor(...), etc. to change the appearance of the pad.
  // and pad()->SetFrameFillColor(...), etc. to change the appearance of the frame in the pad.
  TVirtualPad* pad() const { return m_ppad; }
  bool havePad() const { return pad() != nullptr; }

  // Return the canvas holding this pad. It may be held by an ancestor.
  // If doDraw is true, a draw is done if the canvas does not yet exist.
  TCanvas* canvas() const;
  TCanvas* canvas(bool doDraw =false);
  bool haveCanvas() const { return canvas() != nullptr; }

  // Save the canvas holding this pad to the specified file.
  int print(std::string fname);

  // Return the top-level manipulator, i.e. the ancestor that holds (or would hold)
  // the canvas for this pad.
  TPadManipulator* progenitor();

  // Return the number of subpads.
  unsigned int npad() const { return m_subMans.size(); }

  // Return the manipulator for a subpad.
  // If there are no subpads and ipad=0,  the top manipulator is returned.
  TPadManipulator* man(unsigned int ipad =0) {
    if ( npad() == 0 && ipad == 0 ) return this;
    return ipad<npad() ? &m_subMans[ipad] : nullptr;
  }

  // Return the parent pad.
  TPadManipulator* parent() { return m_parent; }
  bool haveParent() const { return m_parent != nullptr; }

  // Return the primary histogram or graph for this pad.
  TH1* hist() const { return m_ph.get(); }
  bool haveHist() const { return hist() != nullptr; }
  TGraph* graph() const { return m_pg.get(); }
  bool haveGraph() const { return graph() != nullptr; }
  TObject* object() const;
  bool haveHistOrGraph() const { return object() != nullptr; }
  std::string drawOpt() const { return m_dopt; }

  // Return the overlaid objects and options.
  const TObjVector& objects() const { return m_objs; }
  TH1* getHist(unsigned int iobj);
  const std::vector<std::string>& objOpts() const { return m_opts; }

  // Return a histogram (primary or secondary) by name.
  // This is the drawn hist and may be used to changes its visible properties
  // or reference those in a legend.
  TH1* getHist(std::string hnam);

  // Return information about the canvas holding this pad.
  int canvasPixelsX() const;
  int canvasPixelsY() const;

  // Return information about the full drawing area for this pad.
  double padPixelsX() const;
  double padPixelsY() const;

  // Return info about the frame that holds the histogram or graph.
  // The frame does not exist until draw() is called..
  // The frame histogram has the drawn axis.
  TFrame* frame() const;
  bool haveFrame() const { return frame() != nullptr; }
  int framePixelsX() const;
  int framePixelsY() const;
  TH1* frameHist() const;
  bool haveFrameHist() const { return frameHist() != nullptr; }

  // Return the vertical mod lines associated with this pad.
  const std::vector<TLinePtr>& verticalModLines() const { return m_vmlLines; }

  // Add a subpad covering (x1, y1) to (x2, y2) in NDC units, i.e. 0 < x,y < 1.
  int addPad(double x1, double y1, double x2, double y2, int icol =-1);

  // Tile the top pad with nx x ny  or nx x nx subpads.
  int split(Index nx, Index ny);
  int split(Index nx);

  // Add an object (histogram or graph) to the pad.
  // The first object must be a histogram or graph.
  int add(unsigned int ipad, TObject* pobj, std::string sopt ="", bool replace =false);
  int add(TObject* pobj, std::string sopt ="", bool replace =false);

  // Add a legend.
  // This is added to the list of objects.
  TLegend* addLegend(double x1, double y1, double x2, double y2);

  // Set and get the title associated with this pad.
  // The initial value for this is taken from the primary object.
  int setTitle(std::string sttl);
  std::string getTitle() const { return m_title.GetTitle(); }

  // Remove histograms and graphs from top and subpads.
  int clear();

  // Update the coordinates and histogram for this pad.
  int update();

  // Get the axes of the histogram or graph.
  TAxis* getXaxis() const;
  TAxis* getYaxis() const;

  // Set drawing attributes.
  int setCanvasSize(int wx, int wy);
  int setFillColor(int col) { m_fillColor = col; return update(); }
  int setFrameFillColor(int col) { m_frameFillColor = col; return update(); }
  int setGridX(bool flag =true) { m_gridX = flag; return update(); }
  int setGridY(bool flag =true) { m_gridY = flag; return update(); }
  int setGrid(bool flag =true) { m_gridX = flag; m_gridY=flag; return update(); }
  int setLogX(bool flag =true) { m_logX = flag; return update(); }
  int setLogY(bool flag =true) { m_logY = flag; return update(); }
  int setLogZ(bool flag =true) { m_logZ = flag; return update(); }

  // Set the tick lengths. This is a fraction of the axis length.
  // If the X length is zero, those ticks are not drawn.
  // If the Y length is zero, those ticks on both axis have the same pixel lengths and
  // that length is detemined by the shorter axis.
  int setTickLength(double len) { m_tickLengthX = len; m_tickLengthY = 0.0; return update(); }
  int setTickLengthX(double len) { m_tickLengthX = len; return update(); }
  int setTickLengthY(double len) { m_tickLengthY = len; return update(); }

  // Set the displayed ranges.
  // If x1 >= x2 (default), then the range is that of the primary object.
  int setRangeX(double x1, double x2);
  int setRangeY(double y1, double y2);
  int setRangeZ(double y1, double y2);
  int setRangeXY(double x1, double x2, double y1, double y2);

  // Set the displayed ranges for logarithmic axes.
  // If x1 >= x2 or either <=0, then the standard range is used.
  int setLogRangeX(double x1, double x2);
  int setLogRangeY(double y1, double y2);
  int setLogRangeZ(double y1, double y2);

  // Add or remove top and right axis.
  int addAxis(bool flag =true);

  // Add or remove top x-axis.
  int addAxisTop(bool flag =true);

  // Add or remove right y-axis.
  int addAxisRight(bool flag =true);

  // Add underflow bin to the plot (false to remove it).
  int showUnderflow(bool show =true);
  int showOverflow(bool show =true);

  // Return the under/overflow histogram.
  TH1* flowHistogram() { return m_flowHist.get(); }

  // Remove all lines.
  int clearLines();

  // Add a vertical line at x=xoff or horizonatl at y=yoff.
  // The lines are draw with style isty from the low edge to lenfrac*width
  int addVerticalLine(double xoff =0.0, double lenfrac =1.0, int isty =1);
  int addHorizontalLine(double yoff =0.0, double lenfrac =1.0, int isty =1);

  // Add vertical modulus lines.
  // I.e at x = xoff, xoff+/-xmod, xoff+/-2*xmod, ...
  // The lines are draw with style isty from the low edge to lenfrac*width
  int addVerticalModLines(double xmod, double xoff =0.0, double lenfrac =1.0, int isty =3);
  int addHorizontalModLines(double ymod, double yoff =0.0, double lenfrac =1.0, int isty =3);

  // Add histogram function ifun to the pad.
  int addHistFun(unsigned int ifun =0);

  // Fix the BG color for 2D histos to be the same as the lowest color.
  // Otherwise underflows have the color of zeros.
  int fixFrameFillColor();

  // Draw the canvas.
  // Until this is called, nothing will apear on the screen.
  // If this is a subpad, then the parent is drawn.
  int draw();

  // Delete the canvas holding this pad.
  int erase();

  // Redraw a canvas.
  // The current canvas is deleted and a new one created.
  // This is useful if the batch mode has changed or to bring a window to the top.
  int redraw();

  // Draw the top axis.
  int drawAxisTop();

  // Draw right y-axis.
  int drawAxisRight();

  // Draw the hist functions.
  int drawHistFuns();

  // Draw the current lines.
  int drawLines();

private:

  TPadManipulator* m_parent;
  TVirtualPad* m_ppad;
  int m_canvasWidth;
  int m_canvasHeight;
  std::shared_ptr<TH1> m_ph;
  std::shared_ptr<TGraph> m_pg;
  std::string m_dopt;
  TObjVector m_objs;
  std::vector<std::string> m_opts;
  Bounds3 m_setBounds;
  Bounds3 m_setBoundsLog;
  int m_fillColor;
  int m_frameFillColor;
  bool m_gridX;
  bool m_gridY;
  bool m_logX;
  bool m_logY;
  bool m_logZ;
  double m_tickLengthX;
  double m_tickLengthY;
  std::shared_ptr<TH1> m_flowHist;
  bool m_showUnderflow;
  bool m_showOverflow;
  bool m_top;
  bool m_right;
  TLatex m_title;
  std::vector<unsigned int> m_histFuns;
  std::vector<double> m_hmlXmod;
  std::vector<double> m_hmlXoff;
  std::vector<int> m_hmlXStyle;
  std::vector<double> m_hmlXLength;
  std::vector<double> m_vmlXmod;
  std::vector<double> m_vmlXoff;
  std::vector<int> m_vmlXStyle;
  std::vector<double> m_vmlXLength;
  std::vector<std::shared_ptr<TLine>> m_vmlLines;
  BoundsVector m_subBounds;
  std::vector<TPadManipulator> m_subMans;

};

#endif
