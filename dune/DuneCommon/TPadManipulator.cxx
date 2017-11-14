// TPadManipulator.cxx

#include "TPadManipulator.h"
#include <iostream>
#include <sstream>
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TList.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLine.h"
#include "TF1.h"

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using Index = unsigned int;

bool dbg = 0;

//**********************************************************************

TPadManipulator::TPadManipulator()
: m_parent(nullptr), m_ppad(nullptr),
  m_canvasWidth(0), m_canvasHeight(0),
  m_fillColor(0), m_frameFillColor(0),
  m_gridX(false), m_gridY(false),
  m_logX(false), m_logY(false),
  m_showUnderflow(false), m_showOverflow(false),
  m_top(false), m_right(false) {
  const string myname = "TPadManipulator::ctor: ";
  if ( dbg ) cout << myname << this << endl;
}

//**********************************************************************

TPadManipulator::TPadManipulator(TVirtualPad* ppad)
: TPadManipulator() {
  const string myname = "TPadManipulator::ctor(pad): ";
  cout << myname << "Not implemented--will copy objects from pad." << endl;
  if ( dbg ) cout << myname << this << endl;
}

//**********************************************************************

TPadManipulator::TPadManipulator(Index wx, Index wy, Index nPadX, Index nPadY)
: TPadManipulator() {
  m_canvasWidth = wx;
  m_canvasHeight = wy;
  Index nx = nPadX;
  Index ny = nPadY > 0 ? nPadY : nx;
  if ( nx ) split(nx, ny);
}

//**********************************************************************

TPadManipulator::TPadManipulator(const TPadManipulator& rhs)
: TPadManipulator() {
  *this = rhs;
}

//**********************************************************************

TPadManipulator& TPadManipulator::operator=(const TPadManipulator& rhs) {
  if ( this == &rhs ) return *this;
  bool saveAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  if ( m_ppad != nullptr ) {
    m_ppad->Clear();
  }
  m_canvasWidth = rhs.m_canvasWidth;
  m_canvasHeight = rhs.m_canvasHeight;
  // Clone drawn objects.
  m_ph.reset();
  m_pg.reset();
  if ( rhs.m_ph ) {
    m_ph.reset((TH1*)rhs.m_ph->Clone());
  } else if ( rhs.m_pg ) {
    m_pg.reset((TGraph*)rhs.m_pg->Clone());
  }
  m_dopt = rhs.m_dopt;
  m_objs.clear();
  for ( TObjPtr pobj : rhs.m_objs ) {
    m_objs.emplace_back(pobj ? pobj->Clone() : nullptr);
  }
  m_opts = rhs.m_opts;
  TH1::AddDirectory(saveAddDirectory);
  // Copy decoration instructions.
  m_fillColor = rhs.m_fillColor;
  m_frameFillColor = rhs.m_frameFillColor;
  m_gridX = rhs.m_gridX;
  m_gridY = rhs.m_gridY;
  m_logX = rhs.m_logX;
  m_logY = rhs.m_logY;
  m_showUnderflow = rhs.m_showUnderflow;
  m_showOverflow = rhs.m_showOverflow;
  m_top = rhs.m_top;
  m_right = rhs.m_right;
  m_histFuns = rhs.m_histFuns;
  m_hmlXmod = rhs.m_hmlXmod;
  m_hmlXoff = rhs.m_hmlXoff;
  m_hmlXStyle = rhs.m_hmlXStyle;
  m_hmlXLength = rhs.m_hmlXLength;
  m_vmlXmod = rhs.m_vmlXmod;
  m_vmlXoff = rhs.m_vmlXoff;
  m_vmlXStyle = rhs.m_vmlXStyle;
  m_vmlXLength = rhs.m_vmlXLength;
  m_vmlLines.clear();
  m_subMans.clear();
  for ( const TPadManipulator& man : rhs.m_subMans ) {
    m_subMans.emplace_back(man);
  }
  update();
  return *this;
}

//**********************************************************************

TPadManipulator::~TPadManipulator() {
  //cout << "Destructed @" << this << endl;
}

//**********************************************************************

double TPadManipulator::xmin() const {
  if ( m_ppad != nullptr ) return m_ppad->GetUxmin();
  TAxis* pax = getXaxis();
  if ( pax == nullptr ) return 0.0;
  return pax->GetXmin();
}

//**********************************************************************

double TPadManipulator::xmax() const {
  if ( m_ppad != nullptr ) return m_ppad->GetUxmax();
  TAxis* pax = getXaxis();
  if ( pax == nullptr ) return 0.0;
  return pax->GetXmax();
}

//**********************************************************************

double TPadManipulator::ymin() const {
  if ( m_ppad != nullptr ) return m_ppad->GetUymin();
  TAxis* pax = getYaxis();
  if ( pax == nullptr ) return 0.0;
  return pax->GetXmin();
}

//**********************************************************************

double TPadManipulator::ymax() const {
  if ( m_ppad != nullptr ) return m_ppad->GetUymax();
  TAxis* pax = getYaxis();
  if ( pax == nullptr ) return 0.0;
  return pax->GetXmax();
}

//**********************************************************************

int TPadManipulator::addPad(double x1, double y1, double x2, double y2, int icol) {
  if ( x2 <= x1 ) return 1;
  if ( y2 <= y1 ) return 2;
  if ( x1 < 0.0 ) return 3;
  if ( x2 > 1.0 ) return 4;
  if ( y1 < 0.0 ) return 5;
  if ( y2 > 1.0 ) return 6;
  m_subBounds.emplace_back(x1, y1, x2, y2);
  m_subMans.emplace_back();
  m_subMans.back().m_parent = this;
  return 0;
}

//**********************************************************************

int TPadManipulator::split(Index nx, Index ny) {
  if ( nx < 1 ) return 1;
  if ( ny < 1 ) return 2;
  double dx = 1.0/nx;
  double dy = 1.0/ny;
  double y2 = 1.0;
  for ( Index iy=0; iy<ny; ++iy ) {
    double y1 = y2 - dy;
    double x1 = 0.0;
    for ( Index ix=0; ix<nx; ++ix ) {
      double x2 = x1 + dx;
      int rstat = addPad(x1, y1, x2, y2);
      if ( rstat ) return 100 + rstat;
      x1 = x2;
    }
    y2 = y1;
  }
  return 0;
}

//**********************************************************************

int TPadManipulator::split(Index nx) {
  return split(nx, nx);
}

//**********************************************************************

int TPadManipulator::add(Index ipad, TObject* pobj, string sopt, bool replace) {
  const string myname = "TPadManipulator::add: ";
  if ( dbg ) cout << myname << this << endl;
  if ( pobj == nullptr ) return 101;
  TPadManipulator* pman = man(ipad);
  if ( pman == nullptr ) return 102;
  bool haveHistOrGraph = m_ph != nullptr || m_pg != nullptr;
  TH1* ph = dynamic_cast<TH1*>(pobj);
  TGraph* pg = dynamic_cast<TGraph*>(pobj);
  bool isNotHistOrGraph = ph==nullptr && pg==nullptr;
  // If we already have the primary histogram or graph, or this is not one of those,
  // this is an overlaid object.
  if ( (!replace && haveHistOrGraph) || isNotHistOrGraph ) {
    TObject* pobjc = pobj->Clone();
    TH1* phc = dynamic_cast<TH1*>(pobjc);
    if ( phc != nullptr ) phc->SetDirectory(nullptr);
    m_objs.emplace_back(pobjc);
    m_opts.push_back(sopt);
  // Otherwise, the passed object becomes the primary object and must be
  // a histogram or graph.
  } else {
    if ( ph == nullptr && pg == nullptr ) return 103;
    if ( pman->hist() != nullptr ) {
      if ( replace ) pman->m_ph.reset();
      else return 104;
    }
    if ( pman->graph() != nullptr ) {
      if ( replace ) pman->m_pg.reset();
      else return 105;
    }
    //if ( m_ppad == nullptr ) return 106;
    if ( ph != nullptr ) {
      TH1* phc = (TH1*) ph->Clone();
      phc->SetDirectory(nullptr);
      pman->m_ph.reset(phc);
    } else {
      TGraph* pgc = (TGraph*) pg->Clone();
      pman->m_pg.reset(pgc);
    }
    pman->m_dopt = sopt;
  }
  return update();
} 

//**********************************************************************

int TPadManipulator::add(TObject* pobj, string sopt, bool replace) {
  return add(0, pobj, sopt, replace);
}

//**********************************************************************

int TPadManipulator::clear() {
  m_ph.reset();
  m_pg.reset();
  if ( pad() != nullptr && npad() == 0 ) pad()->Clear();
  for ( TPadManipulator& man : m_subMans ) man.clear();
  return 0;
}

//**********************************************************************

int TPadManipulator::update() {
  const string myname = "TPadManipulator::update: ";
  if ( dbg ) cout << myname << this << endl;
  int rstat = 0;
  // No action if cavas/pad is not yet created.
  if ( m_ppad == nullptr ) return 0;
  // Update the subpads.
  // If needed, create the TPad for each subpad.
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  if ( m_subMans.size() ) {
    for ( Index isub=0; isub<m_subMans.size(); ++isub ) {
      TPadManipulator& man = m_subMans[isub];
      if ( man.m_ppad == nullptr ) {
        if ( m_subBounds.size() <= isub ) {
          cout << myname << "ERROR: Too few bounds." << endl;
          break;
        }
        const Bounds& bnd = m_subBounds[isub];
        ostringstream ssnam;
        ssnam << m_ppad->GetName() << "_sub" << isub;
        string spnam = ssnam.str();
        string spttl = spnam;
        TPad* ppad = new TPad(spnam.c_str(), spttl.c_str(), bnd.x1, bnd.y1, bnd.x2, bnd.y2);
        ppad->Draw();
        man.m_ppad = ppad;
      }
      int srstat = man.update();
      if ( srstat ) {
        cout << myname << "WARNING: Error " << srstat << " updating pad " << isub << endl;
        rstat += srstat;
      }
    }
  }
  // Set the TPad attributes.
  m_ppad->SetFillColor(m_fillColor);
  m_ppad->SetFrameFillColor(m_frameFillColor);
  m_ppad->SetGridx(m_gridX);
  m_ppad->SetGridy(m_gridY);
  m_ppad->SetLogx(m_logX);
  m_ppad->SetLogy(m_logY);
  // Make sure the axis range are up to date before fetching them.
  if ( m_ph != nullptr ) m_ph->Draw("AXIS");
  else if ( m_pg != nullptr ) m_pg->Draw("A");
  //gPad->Update();
  string sopt = "-US";
  const TList* prims = m_ppad->GetListOfPrimitives();
  bool noHist = m_ph == nullptr;
  bool noGraph = m_pg == nullptr;
  bool setMargins = noHist && noGraph;
  setMargins = true;  // Do this always so it is also called when primary is added
  if ( noHist && noGraph ) {
    for ( TObjLink* plnk=prims->FirstLink(); plnk; plnk=plnk->Next() ) {
      TObject* pobj = plnk->GetObject();
      string sopt = plnk->GetOption();
      TH1* ph = dynamic_cast<TH1*>(pobj);
      if ( ph != nullptr )  {
        m_ph.reset(dynamic_cast<TH1*>(ph->Clone("hmanip0")));
        m_ph->SetDirectory(nullptr);
        noHist = false;
        m_dopt = sopt;
        break;
      }
      TGraph* pg = dynamic_cast<TGraph*>(pobj);
      if ( pg != nullptr )  {
        m_pg.reset(dynamic_cast<TGraph*>(pg->Clone("gmanip0")));
        noGraph = false;
        m_dopt = sopt;
        break;
      }
    }
  }
  if ( noHist && noGraph ) {
    //cout << myname << "Pad does not have a histogram or graph." << endl;
    gPad = pPadSave;
    return 0;
  }
  // Set margins the first time the histogram is found.
  // After this, user can override with pad()->SetRightMargin(...), ...
  bool isTH = m_ph != nullptr;
  bool isTH2 = dynamic_cast<TH2*>(m_ph.get()) != nullptr;
  bool isTH1 = isTH && !isTH2;
  if ( setMargins ) {
    if ( isTH2 ) m_ppad->SetRightMargin(0.10);
    else m_ppad->SetRightMargin(0.03);
    m_ppad->SetLeftMargin(0.12);
    m_ppad->SetTopMargin(0.07);
    getXaxis()->SetTitleOffset(1.20);
    if ( isTH ) m_ph->SetTitleOffset(1.20);
  }
  int nbin = isTH1 ? m_ph->GetNbinsX() : 0;
  int flowcol = kAzure - 9;
  // Redraw everything.
  m_ppad->Clear();
  if ( isTH ) m_ph->Draw(m_dopt.c_str());
  else m_pg->Draw(m_dopt.c_str());
  if ( (m_showUnderflow || m_showOverflow) && nbin > 0 ) {
    if ( m_flowHist == nullptr ) {
      m_flowHist.reset((TH1*) m_ph->Clone("hmaniptmp"));
      m_flowHist->SetDirectory(nullptr);
      m_flowHist->SetStats(0);
      m_flowHist->SetLineColor(flowcol);
      m_flowHist->SetFillColor(flowcol);
      m_flowHist->SetLineWidth(1);
    }
    m_flowHist->Reset();
    if ( m_showUnderflow ) m_flowHist->SetBinContent(1, m_ph->GetBinContent(0));
    if ( m_showOverflow ) m_flowHist->SetBinContent(nbin, m_ph->GetBinContent(nbin+1));
    m_flowHist->DrawCopy("same");
    m_ph->Draw("same");
    m_ph->Draw("axis same");
  }
  drawHistFuns();
  for ( Index iobj=0; iobj<m_objs.size(); ++iobj ) {
    TObject* pobj = m_objs[iobj].get();
    string sopt = m_opts[iobj];
    if ( pobj != nullptr ) pobj->Draw(sopt.c_str());
  }
  m_ppad->Update();   // Need an update here to get correct results for xmin, ymin, xmax, ymax
  drawLines();
  if ( m_top ) drawAxisTop();
  if ( m_right ) drawAxisRight();
  pad()->RedrawAxis();  // In case they are covered
  gPad = pPadSave;
  return rstat;
}

//**********************************************************************

TAxis* TPadManipulator::getXaxis() const {
  if ( m_ph != nullptr ) return m_ph->GetXaxis();
  if ( m_pg != nullptr ) return m_pg->GetXaxis();
  return nullptr;
}

//**********************************************************************

TAxis* TPadManipulator::getYaxis() const {
  if ( m_ph != nullptr ) return m_ph->GetYaxis();
  if ( m_pg != nullptr ) return m_pg->GetYaxis();
  return nullptr;
}

//**********************************************************************

int TPadManipulator::setRangeX(double x1, double x2) {
  TAxis* pax = getXaxis();
  if ( pax == nullptr ) return 1;
  if ( x2 <= x1 ) return 2;
  if ( graph() != nullptr ) pax->SetLimits(x1, x2);
  else pax->SetRangeUser(x1, x2);
  return update();
}

//**********************************************************************

int TPadManipulator::setRangeY(double x1, double x2) {
  TAxis* pax = getYaxis();
  if ( pax == nullptr ) return 1;
  if ( x2 <= x1 ) return 2;
  pax->SetRangeUser(x1, x2);
  return update();
}

//**********************************************************************

int TPadManipulator::setRanges(double x1, double x2, double y1, double y2) {
  TAxis* pax = getXaxis();
  TAxis* pay = getYaxis();
  if ( pax == nullptr || pay == nullptr ) return 1;
  if ( x2 <= x1 ) return 2;
  if ( y2 <= y1 ) return 2;
  if ( graph() != nullptr ) pax->SetLimits(x1, x2);
  else pax->SetRangeUser(x1, x2);
  pay->SetRangeUser(y1, y2);
  return update();
}

//**********************************************************************

int TPadManipulator::addAxis(bool flag) {
  return addAxisTop(flag) + addAxisRight(flag);
}

//**********************************************************************

int TPadManipulator::addAxisTop(bool flag) {
  m_top = flag;
  return update();
}

//**********************************************************************

int TPadManipulator::addAxisRight(bool flag) {
  m_right = flag;
  return update();
}

//**********************************************************************

int TPadManipulator::addHistFun(unsigned int ifun) {
  m_histFuns.push_back(ifun);
  return update();
}

//**********************************************************************

int TPadManipulator::showUnderflow(bool show) {
  if ( show == m_showUnderflow ) return 0;
  m_showUnderflow = show;
  return update();
}

//**********************************************************************

int TPadManipulator::showOverflow(bool show) {
  if ( show == m_showOverflow ) return 0;
  m_showOverflow = show;
  return update();
}

//**********************************************************************

int TPadManipulator::clearLines() {
  m_vmlXmod.clear();
  m_vmlXoff.clear();
  m_vmlXStyle.clear();
  m_vmlXLength.clear();
  m_hmlXmod.clear();
  m_hmlXoff.clear();
  m_hmlXStyle.clear();
  m_hmlXLength.clear();
  return update();
}

//**********************************************************************

int TPadManipulator::addVerticalLine(double xoff, double lenfrac, int isty) {
  m_vmlXmod.push_back(0.0);
  m_vmlXoff.push_back(xoff);
  m_vmlXStyle.push_back(isty);
  m_vmlXLength.push_back(lenfrac);
  return update();
}
//**********************************************************************

int TPadManipulator::addHorizontalLine(double xoff, double lenfrac, int isty) {
  m_hmlXmod.push_back(0.0);
  m_hmlXoff.push_back(xoff);
  m_hmlXStyle.push_back(isty);
  m_hmlXLength.push_back(lenfrac);
  drawLines();
  return 0;
}
//**********************************************************************

int TPadManipulator::addVerticalModLines(double xmod, double xoff, double lenfrac, int isty) {
  m_vmlXmod.push_back(xmod);
  m_vmlXoff.push_back(xoff);
  m_vmlXStyle.push_back(isty);
  m_vmlXLength.push_back(lenfrac);
  drawLines();
  return 0;
}

//**********************************************************************

int TPadManipulator::addHorizontalModLines(double xmod, double xoff, double lenfrac, int isty) {
  m_hmlXmod.push_back(xmod);
  m_hmlXoff.push_back(xoff);
  m_hmlXStyle.push_back(isty);
  m_hmlXLength.push_back(lenfrac);
  drawLines();
  return 0;
}

//**********************************************************************

int TPadManipulator::fixFrameFillColor() {
  if ( gStyle == nullptr ) return 2;
  m_ppad->SetFrameFillColor(gStyle->GetColorPalette(0));
  return 0;
}

//**********************************************************************

int TPadManipulator::draw() {
  if ( m_parent != nullptr ) return m_parent->draw();
  if ( m_ppad == nullptr ) {
    TCanvas* pcan = new TCanvas;
    if ( m_canvasWidth > 0 && m_canvasHeight > 0 ) {
      string snam = pcan->GetName();
      string sttl = pcan->GetTitle();
      delete pcan;
      pcan = new TCanvas(snam.c_str(), sttl.c_str(), m_canvasWidth, m_canvasHeight);
      //pcan->SetCanvasSize(m_canvasWidth, m_canvasHeight);
      //pcan->SetWindowSize(m_canvasWidth, m_canvasHeight);
    }
    m_ppad = pcan;
  }
  return update();
}

//**********************************************************************

int TPadManipulator::drawAxisTop() {
  if ( ! m_top ) return 0;
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  double xminPad, yminPad, xmaxPad, ymaxPad;
  m_ppad->GetRangeAxis(xminPad, yminPad, xmaxPad, ymaxPad);
  string sopt = "-US";
  double x1 = xmin();
  double x2 = xmax();
  if ( m_ppad->GetLogx() ) {
    xminPad = pow(10.0, xminPad);
    xmaxPad = pow(10.0, xmaxPad);
    x1 = pow(10.0, x1);
    x2 = pow(10.0, x2);
    sopt += "G";
  }
  if ( m_ppad->GetLogy() ) {
    yminPad = pow(10.0, yminPad);
    ymaxPad = pow(10.0, ymaxPad);
  }
  TAxis* paxold = getXaxis();
  if ( paxold == 0 ) return 2;
  double ticksize = paxold->GetTickLength();
  int ndiv = paxold->GetNdivisions();
  TGaxis* paxnew = new TGaxis(xminPad, ymaxPad, xmaxPad, ymaxPad,
                              x1, x2, ndiv, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickLength(ticksize);
  string name = "TopAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = m_ppad->GetListOfPrimitives();
  for ( int iobj=pobjs->GetEntries()-1; iobj>=0; --iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold == nullptr ) continue;
    if ( paxold->GetName() == name ) {
      pobjs->RemoveAt(iobj);
      break;
    }
  }
  paxnew->Draw();
  gPad = pPadSave;
  return 0;
}

//**********************************************************************

int TPadManipulator::drawAxisRight() {
  if ( ! m_right ) return 0;
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  double xminPad, yminPad, xmaxPad, ymaxPad;
  m_ppad->GetRangeAxis(xminPad, yminPad, xmaxPad, ymaxPad);
  double y1 = ymin();
  double y2 = ymax();
  string sopt = "+US";
  if ( m_ppad->GetLogx() ) {
    xminPad = pow(10.0, xminPad);
    xmaxPad = pow(10.0, xmaxPad);
  }
  if ( m_ppad->GetLogy() ) {
    yminPad = pow(10.0, yminPad);
    ymaxPad = pow(10.0, ymaxPad);
    y1 = pow(10.0, y1);
    y2 = pow(10.0, y2);
    sopt += "G";
  }
  TAxis* paxold = getYaxis();
  if ( paxold == 0 ) return 2;
  double ticksize = paxold->GetTickLength();
  int ndiv = paxold->GetNdivisions();
  TGaxis* paxnew = new TGaxis(xmaxPad, yminPad, xmaxPad, ymaxPad,
                              y1, y2, ndiv, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickLength(ticksize);
  string name = "RightAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = m_ppad->GetListOfPrimitives();
  for ( int iobj=pobjs->GetEntries()-1; iobj>=0; --iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold == nullptr ) continue;
    if ( paxold->GetName() == name ) {
      pobjs->RemoveAt(iobj);
      break;
    }
  }
  paxnew->Draw("");
  gPad = pPadSave;
  return 0;
}

//**********************************************************************

int TPadManipulator::drawHistFuns() {
  if ( m_ppad == nullptr ) return 1;
  const TList* pfuns = nullptr;
  if      ( m_ph != nullptr ) pfuns = m_ph->GetListOfFunctions();
  else if ( m_pg != nullptr ) pfuns = m_pg->GetListOfFunctions();
  else    return 1;
  const TList& funs = *pfuns;
  unsigned int nfun = funs.GetEntries();
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  for ( unsigned int ifun : m_histFuns ) {
    if ( ifun >= nfun ) continue;
    TF1* pfun = dynamic_cast<TF1*>(funs.At(ifun));
    pfun->Draw("same");
  }
  gPad = pPadSave;
  return 0;
}

//**********************************************************************

int TPadManipulator::drawLines() {
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  m_vmlLines.clear();
  for ( Index iset=0; iset<m_vmlXmod.size(); ++iset ) {
    double xmod = m_vmlXmod[iset];
    double xoff = m_vmlXoff[iset];
    double xlen = m_vmlXLength[iset];
    int isty = m_vmlXStyle[iset];
    if ( xmod < 0.0 ) continue;
    double x = xoff;
    if ( xmod > 0.0 ) {
      while ( x >= xmin() ) x -= xmod;
      while ( x < xmin() ) x += xmod;
    }
    double ytop = ymax() + (xlen - 1.0)*(ymax() - ymin());
    while ( x <= xmax() ) {
      std::shared_ptr<TLine> pline(new TLine(x, ymin(), x, ytop));
      pline->SetLineStyle(isty);
      m_vmlLines.push_back(pline);
      pline->Draw();
      if ( xmod == 0.0 ) break;
      x += xmod;
    }
  }
  for ( Index iset=0; iset<m_hmlXmod.size(); ++iset ) {
    double ymod = m_hmlXmod[iset];
    double yoff = m_hmlXoff[iset];
    double ylen = m_hmlXLength[iset];
    int isty = m_hmlXStyle[iset];
    if ( ymod < 0.0 ) continue;
    double y = yoff;
    if ( ymod > 0.0 ) {
      while ( y >= ymin() ) y -= ymod;
      while ( y < ymin() ) y += ymod;
    }
    double xtop = xmax() + (ylen - 1.0)*(xmax() - xmin());
    while ( y <= ymax() ) {
      std::shared_ptr<TLine> pline(new TLine(xmin(), y, xtop, y));
      pline->SetLineStyle(isty);
      m_vmlLines.push_back(pline);
      pline->Draw();
      if ( ymod == 0.0 ) break;
      y += ymod;
    }
  }
  gPad = pPadSave;
  return 0;
}

//**********************************************************************
