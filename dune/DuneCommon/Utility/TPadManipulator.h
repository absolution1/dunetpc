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
#include "TLegend.h"
#include "TExec.h"

class TVirtualPad;
class TCanvas;
class TH1;
class TGraph;
class TAxis;
class TLegend;
class TFrame;
class TBuffer;
class TDirectory;

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
  using LineVector = std::vector<TLine*>;
  using TObjVector = std::vector<TObject*>;
  using BoundsVector = std::vector<Bounds>;
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using HistPtr = std::shared_ptr<TH1>;

  // Get object with name onam from Root directory rdir.
  static TPadManipulator* get(Name onam ="tpad", TDirectory* tdir =nullptr);

  // Read object with name onam from the Root file fnam.
  static TPadManipulator* read(Name fnam, Name onam ="tpad");

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

  // Save this object with name onam to root directory tdir.
  // If tdir is null, the current Root directory is used.
  int put(Name onam ="tpad", TDirectory* tdir =nullptr) const;

  // Save this object with name onam to root file fnam.
  // If tdir is null, the current Root directory is used.
  // The root file is opened in update mode and then closed.
  int write(Name fnam, Name onam ="tpad") const;

  // Create an image file. Name suffix should be a known format: png, pdf, ....
  // If the suffix is .root or .tpad, then the object is written to that file.
  int printOnce(std::string fname);

  // Print one or more image files using StringManipulator pattern split.
  // E.g. myfile.{png,pdf,tpad} --> myfile.png, myfile.pdf and myfile.tpad.
  int print(std::string fname, std::string spat ="{,}");

  // Return the top-level manipulator, i.e. the ancestor that holds (or would hold)
  // the canvas for this pad.
  TPadManipulator* progenitor();

  // Return the number of subpads.
  unsigned int npad() const;

  // Return the manipulator for a subpad.
  // If there are no subpads and ipad=0,  the top manipulator is returned.
  TPadManipulator* man(unsigned int ipad =0);

  // Return the parent pad.
  TPadManipulator* parent() { return m_parent; }
  const TPadManipulator* parent() const { return m_parent; }
  bool haveParent() const { return m_parent != nullptr; }

  // Return the primary histogram or graph for this pad.
  TH1* hist() const { return m_ph; }
  bool haveHist() const { return hist() != nullptr; }
  TGraph* graph() const { return m_pg; }
  bool haveGraph() const { return graph() != nullptr; }
  TObject* object() const;
  bool haveHistOrGraph() const { return object() != nullptr; }
  std::string drawOpt() const { return m_dopt; }

  // Return the overlaid objects and options.
  const TObjVector& objects() const { return m_objs; }
  TObject* object(Index iobj) const { return iobj<objects().size() ? objects()[iobj] : nullptr; }
  TH1* getHist(unsigned int iobj);
  const std::vector<std::string>& objOpts() const { return m_opts; }

  // Return a histogram (primary or secondary) by name.
  // This is the drawn hist and may be used to changes its visible properties
  // or reference those in a legend.
  TH1* getHist(std::string hnam);

  // Return the last object added to the pad.
  TObject* lastObject() const;

  // Return information about the canvas holding this pad.
  int canvasPixelsX() const;
  int canvasPixelsY() const;

  // Return information about the full drawing area for this pad.
  // Pad must be drawn (or printed) to set these values for subpads.
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

  // Delete the line objects.
  void clearLineObjects();

  // Return the lines associated with this pad.
  const LineVector& verticalModLines() const { return m_lines; }

  // Add a subpad covering (x1, y1) to (x2, y2) in NDC units, i.e. 0 < x,y < 1.
  int addPad(double x1, double y1, double x2, double y2, int icol =-1);

  // Tile the top pad with nx x ny  or nx x nx subpads.
  int split(Index nx, Index ny);
  int split(Index nx);

  // Add an object (histogram or graph) to the pad.
  // The first object must be a histogram or graph.
  int add(unsigned int ipad, TObject* pobj, std::string sopt ="", bool replace =false);
  int add(TObject* pobj, std::string sopt ="", bool replace =false);

  // Set and fetch the TExec command run when the pad is drawn.
  void setTExec(Name com) { m_exec.SetTitle(com.c_str()); }
  const TExec& getTExec() const { return m_exec; }

  // Set the palette for this pad.
  // Added to TExec.
  // Negative for no platte command.
  void setPalette(int ipal);

  // Set margins. Negative value uses default margin.
  // The default is set here to be reasonable for a wide range of aspect ratio.
  // The margin is the fraction of the pad between the frame and pad edges.
  // The corresponding axis title offset is also changed so that the title
  // does not move relative to the pad, i.e. remains close to the pad edge.
  // This behavior is not (yet) implemented for the z-axis (right) or title (top).
  void setMarginLeft(double xmar) { m_marginLeft = xmar; }
  void setMarginRight(double xmar) { m_marginRight = xmar; }
  void setMarginBottom(double xmar) { m_marginBottom = xmar; }
  void setMarginTop(double xmar) { m_marginTop = xmar; }

  // Center axis labels.
  void centerAxisTitles(bool center =true) { m_axisTitleOpt = center; }

  // Add a legend.
  // This is added to the list of objects.
  TLegend* addLegend(double x1, double y1, double x2, double y2);

  // Return the legend.
  TLegend* getLegend() const { return dynamic_cast<TLegend*>(object(m_iobjLegend)); }

  // Set and get the title associated with this pad.
  // The initial value for this is taken from the primary object.
  // The title is positioned in the middle of the top margin and its
  // size is limited to 90% of the margin height.
  int setTitle(std::string sttl, float height =-1.0);
  std::string getTitle() const { return m_title.GetTitle(); }

  // Set and get the label associated with this pad.
  // The label object can be used to modify the label.
  // This is written in the lower left corner.
  int setLabel(std::string slab);
  std::string getLabel() const { return m_label.GetTitle(); }
  TLatex& getLabelObject() { return m_label; }

  // Remove histograms and graphs from top and subpads.
  int clear();

  // Update the coordinates and histogram for this pad.
  int update();

  // Get the axes of the histogram or graph.
  TAxis* getXaxis() const;
  TAxis* getYaxis() const;
  TAxis* getZaxis() const;

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

  // Set the tick/label division guides.
  int setNdivisionsX(int ndiv) { m_ndivX = ndiv; return 0; };
  int setNdivisionsY(int ndiv) { m_ndivY = ndiv; return 0; };

  // Set the axis label sizes. Default is zero.
  // Unless set explcitely, Sub-pads will use a scaled version of this
  // so all labels are the same size.
  // The z-axis label size is the same as that for y.
  int setLabelSizeX(double siz) { m_labSizeX = siz; return 0; }
  int setLabelSizeY(double siz) { m_labSizeY = siz; return 0; }
  int setTitleSize(double siz) { m_ttlSize = siz; return 0; }

  // Fetch the label size for this pad.
  // Returns the set value if nonzero, otherwise
  // if parent, returns the parent value scaled by ratio of pad heights
  // otherwise returns zero.
  double getLabelSizeX() const;
  double getLabelSizeY() const;
  double getTitleSize() const;

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

  // Set the time offset in second for axes using time format.
  // 0 (default) is UTC (GMT).
  int setTimeOffset(double toff);

  // Set the time format for an axis.
  // If blank (default), time format is not used.
  // The format is similar to that of strftime.
  // See TAxis::SetTimeFormat for more information.
  // Append "%Fyyyy-mm-dd hh:mm:ss" to set the time offset.
  int setTimeFormatX(std::string sfmt);
  int setTimeFormatY(std::string sfmt);

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
  TH1* flowHistogram() { return m_flowHist; }

  // Show overflow points for graphs.
  // Any point off scale in any of the indicated directions is drawn at
  // that boundary with the indicated marker and color.
  // Direction: B=bottom, T=top, L=left, R=right, "" to show nothing.
  int showGraphOverflow(std::string sopt ="BTLR", int imrk =38, int icol =1);

  // Return the under/overflow graph.
  TGraph* flowGraph() { return m_flowGraph; }

  // Remove all lines.
  int clearLines();

  // Add a vertical line at x=xoff or horizontal at y=yoff.
  // The lines are draw with style isty from the low edge to lenfrac*width
  int addVerticalLine(double xoff =0.0, double lenfrac =1.0, int isty =1);
  int addHorizontalLine(double yoff =0.0, double lenfrac =1.0, int isty =1);

  // Add a line with slop dydx and y = y0 at x = 0.
  // The lines are draw with style isty between frame boundaries.
  int addSlopedLine(double slop, double yoff =0.0, int isty =1);

  // Add vertical modulus lines.
  // I.e at x = xoff, xoff+/-xmod, xoff+/-2*xmod, ...
  // The lines are draw with style isty from the low edge to lenfrac*width
  int addVerticalModLines(double xmod, double xoff =0.0, double lenfrac =1.0, int isty =3);
  int addHorizontalModLines(double ymod, double yoff =0.0, double lenfrac =1.0, int isty =3);

  // Set text bin labels. 
  // Only used if primary object is a histogram (2D for y).
  int setBinLabelsX(const NameVector& labs);
  int setBinLabelsY(const NameVector& labs);

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

  // Custom streamer.
  void Streamer(TBuffer& buf);

private:

  // Set the parent for all children.
  // If recurse is true, the operation is also performed on their children.
  void setParents(bool recurse);

private:

  TPadManipulator* m_parent;  //! ==> Do not stream
  TVirtualPad* m_ppad;  //! ==> Do not stream.
  int m_canvasWidth;
  int m_canvasHeight;
  TExec m_exec;
  double m_marginLeft;
  double m_marginRight;
  double m_marginBottom;
  double m_marginTop;
  TH1* m_ph;
  TGraph* m_pg;
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
  int m_ndivX;
  int m_ndivY;
  double m_labSizeX;
  double m_labSizeY;
  double m_ttlSize;
  TH1* m_flowHist;        //! ==> Do not stream.
  bool m_showUnderflow;
  bool m_showOverflow;
  TGraph* m_flowGraph;    //! ==> Do not stream.
  std::string m_gflowOpt;
  int m_gflowMrk;
  int m_gflowCol;
  bool m_top;
  bool m_right;
  TLatex m_title;
  TLatex m_label;
  std::vector<unsigned int> m_histFuns;
  std::vector<double> m_hmlXmod;
  std::vector<double> m_hmlXoff;
  std::vector<int> m_hmlXStyle;
  std::vector<double> m_hmlXLength;
  std::vector<double> m_vmlXmod;
  std::vector<double> m_vmlXoff;
  std::vector<int> m_vmlXStyle;
  std::vector<double> m_vmlXLength;
  std::vector<double> m_slSlop;
  std::vector<double> m_slYoff;
  std::vector<int> m_slStyl;
  LineVector m_lines;  //! ==> Do not stream.
  NameVector m_binLabelsX;
  NameVector m_binLabelsY;
  double m_timeOffset =0.0;
  std::string m_timeFormatX;
  std::string m_timeFormatY;
  BoundsVector m_subBounds;
  std::vector<TPadManipulator> m_subMans;
  Index m_iobjLegend;
  int m_axisTitleOpt;

};

#endif
