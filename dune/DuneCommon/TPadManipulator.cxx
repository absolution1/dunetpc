// TPadManipulator.cxx

#include "TPadManipulator.h"
#include <iostream>
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
using Index = unsigned int;

//**********************************************************************

TPadManipulator::TPadManipulator(TVirtualPad* ppad)
: m_ppad(ppad == nullptr ? gPad : ppad),
  m_showUnderflow(false), m_showOverflow(false),
  m_top(false), m_right(false),
  m_vmlXmod(0.0), m_vmlXoff(0.0) {
  update();
  //cout << "Constructed @" << this << endl;
}

//**********************************************************************

TPadManipulator::TPadManipulator(Index wx, Index wy, Index nPadX, Index nPadY)
: m_ppad(new TCanvas("mycan", "mycan", wx, wy)),
  m_showUnderflow(false), m_showOverflow(false),
  m_top(false), m_right(false),
  m_vmlXmod(0.0), m_vmlXoff(0.0) {
  Index nx = nPadX;
  Index ny = nPadY > 0 ? nPadY : nx;
  if ( nx ) split(nx, ny);
}

//**********************************************************************

TPadManipulator::~TPadManipulator() {
  //cout << "Destructed @" << this << endl;
}

//**********************************************************************

int TPadManipulator::addPad(double x1, double y1, double x2, double y2, int icol) {
  if ( x2 <= x1 ) return 1;
  if ( y2 <= y1 ) return 2;
  if ( x1 < 0.0 ) return 3;
  if ( x2 > 1.0 ) return 4;
  if ( y1 < 0.0 ) return 5;
  if ( y2 > 1.0 ) return 6;
  string snam = "pad";
  string sttl = "pad";
  TPad* ppad = new TPad(snam.c_str(), sttl.c_str(), x1, y1, x2, y2, icol);
  m_ppad->cd();
  ppad->Draw();
  m_subMans.emplace_back(ppad);
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
  const string myname = "TPadManipulator::update: ";
  TH1* ph = dynamic_cast<TH1*>(pobj);
  TGraph* pg = dynamic_cast<TGraph*>(pobj);
  if ( ph == nullptr && pg == nullptr ) return 1;
  TPadManipulator* pman = man(ipad);
  if ( pman == nullptr ) return 2;
  if ( pman->hist() != nullptr ) {
    if ( replace ) pman->m_ph.reset();
    else return 3;
  }
  if ( pman->graph() != nullptr ) {
    if ( replace ) pman->m_pg.reset();
    else return 4;
  }
  if ( m_ppad == nullptr ) m_ppad = new TCanvas;
  pman->pad()->cd();
  pobj->Draw(sopt.c_str());
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
  if ( m_subMans.size() ) {
    int rstat = 0;
    for ( TPadManipulator& man : m_subMans ) rstat += man.update();
    return rstat;
  }
  if ( m_ppad == nullptr ) m_ppad = gPad;
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = m_ppad == gPad ? nullptr : gPad;
  if ( pPadSave != nullptr ) m_ppad->cd();
  gPad->Update();
  gPad->GetRangeAxis(m_xminPad, m_yminPad, m_xmaxPad, m_ymaxPad);
  m_xmin = gPad->GetUxmin();
  m_xmax = gPad->GetUxmax();
  m_ymin = gPad->GetUymin();
  m_ymax = gPad->GetUymax();
  string sopt = "-US";
  if ( gPad->GetLogx() ) {
    m_xminPad = pow(10.0, m_xminPad);
    m_xmaxPad = pow(10.0, m_xmaxPad);
    m_xmin = pow(10.0, m_xmin);
    m_xmax = pow(10.0, m_xmax);
  }
  if ( gPad->GetLogy() ) {
    m_yminPad = pow(10.0, m_yminPad);
    m_ymaxPad = pow(10.0, m_ymaxPad);
    m_ymin = pow(10.0, m_ymin);
    m_ymax = pow(10.0, m_ymax);
  }
  const TList* prims = gPad->GetListOfPrimitives();
  bool noHist = m_ph == nullptr;
  bool noGraph = m_pg == nullptr;
  bool firstDraw = noHist && noGraph;
  if ( noHist && noGraph ) {
    for ( int iprm=0; iprm<prims->GetEntries(); ++iprm ) {
      TObject* pobj = prims->At(iprm);
      TH1* ph = dynamic_cast<TH1*>(pobj);
      if ( ph != nullptr )  {
        m_ph.reset(dynamic_cast<TH1*>(ph->Clone("hmanip0")));
        m_ph->SetDirectory(nullptr);
        noHist = false;
        break;
      }
      TGraph* pg = dynamic_cast<TGraph*>(pobj);
      if ( pg != nullptr )  {
        m_pg.reset(dynamic_cast<TGraph*>(pg->Clone("gmanip0")));
        noGraph = false;
        break;
      }
    }
  }
  if ( noHist && noGraph ) {
    //cout << myname << "Pad does not have a histogram or graph." << endl;
    return 0;
  }
  // Set margins the first time the histogram is found.
  // After this, user can override with pad()->SetRightMargin(...), ...
  bool isTH = m_ph != nullptr;
  bool isTH2 = dynamic_cast<TH2*>(m_ph.get()) != nullptr;
  bool isTH1 = isTH && !isTH2;
  if ( firstDraw ) {
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
  gPad->Clear();
  if ( isTH ) m_ph->Draw();
  else m_pg->Draw("AP");
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
  if ( m_top ) drawAxisTop();
  if ( m_right ) drawAxisRight();
  if ( m_vmlXmod >= 0.0 ) drawVerticalModLines();
  pad()->RedrawAxis();  // In case they are covered
  if ( pPadSave != nullptr ) pPadSave->cd();
  return 0;
}

//**********************************************************************

TAxis* TPadManipulator::getXaxis() {
  if ( m_ph != nullptr ) return m_ph->GetXaxis();
  if ( m_pg != nullptr ) return m_pg->GetXaxis();
  return nullptr;
}

//**********************************************************************

TAxis* TPadManipulator::getYaxis() {
  if ( m_ph != nullptr ) return m_ph->GetYaxis();
  if ( m_pg != nullptr ) return m_pg->GetYaxis();
  return nullptr;
}

//**********************************************************************

int TPadManipulator::setRangeX(double x1, double x2) {
  TAxis* pax = getXaxis();
  if ( pax == nullptr ) return 1;
  if ( x2 <= x1 ) return 2;
  pax->SetRangeUser(x1, x2);
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
  pax->SetRangeUser(x1, x2);
  pay->SetRangeUser(y1, y2);
  return update();
}

//**********************************************************************

int TPadManipulator::addAxis() {
  return addAxisTop() + addAxisRight();
}

//**********************************************************************

int TPadManipulator::addAxisTop() {
  double ticksize = 0;
  int ndiv = 0;
  TAxis* paxold = getXaxis();
  if ( paxold == nullptr ) return 2;
  ticksize = paxold->GetTickLength();
  ndiv = paxold->GetNdivisions();
  return addAxisTop(ticksize, ndiv);
}

//**********************************************************************

int TPadManipulator::addAxisTop(double ticksize, int ndiv) {
  m_top = true;
  m_topTicksize = ticksize;
  m_topNdiv = ndiv;
  return drawAxisTop();
}

//**********************************************************************

int TPadManipulator::drawAxisTop() {
  if ( ! m_top ) return 0;
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = m_ppad == gPad ? nullptr : gPad;
  if ( pPadSave != nullptr ) m_ppad->cd();
  string sopt = "-US";
  if ( gPad->GetLogx() ) sopt += "G";
  double ticksize = m_topTicksize;
  int ndiv = m_topNdiv;
  TGaxis* paxnew = new TGaxis(m_xminPad, m_ymaxPad, m_xmaxPad, m_ymaxPad,
                              m_xmin, m_xmax, ndiv, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickLength(ticksize);
  string name = "TopAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = gPad->GetListOfPrimitives();
  for ( int iobj=pobjs->GetEntries()-1; iobj>=0; --iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold == nullptr ) continue;
    if ( paxold->GetName() == name ) {
      pobjs->RemoveAt(iobj);
      break;
    }
  }
  paxnew->Draw();
  if ( pPadSave != nullptr ) pPadSave->cd();
  return 0;
}

//**********************************************************************

int TPadManipulator::addAxisRight() {
  double ticksize = 0;
  int ndiv = 0;
  TAxis* paxold = getYaxis();
  if ( paxold == nullptr ) return 2;
  ticksize = paxold->GetTickLength();
  ndiv = paxold->GetNdivisions();
  return addAxisRight(ticksize, ndiv);
}

//**********************************************************************

int TPadManipulator::addAxisRight(double ticksize, int ndiv) {
  m_right = true;
  m_rightTicksize = ticksize;
  m_rightNdiv = ndiv;
  return drawAxisRight();
}

//**********************************************************************

int TPadManipulator::drawAxisRight() {
  if ( ! m_right ) return 0;
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = m_ppad == gPad ? nullptr : gPad;
  if ( pPadSave != nullptr ) m_ppad->cd();
  if ( gPad == 0 ) return 2;
  gPad->Update();
  string sopt = "+US";
  if ( gPad->GetLogy() ) sopt += "G";
  double ticksize = m_rightTicksize;
  int ndiv = m_rightNdiv;
  TGaxis* paxnew = new TGaxis(m_xmaxPad, m_yminPad, m_xmaxPad, m_ymaxPad,
                              m_ymin, m_ymax, ndiv, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickLength(ticksize);
  string name = "RightAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = gPad->GetListOfPrimitives();
  for ( int iobj=pobjs->GetEntries()-1; iobj>=0; --iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold == nullptr ) continue;
    if ( paxold->GetName() == name ) {
      pobjs->RemoveAt(iobj);
      break;
    }
  }
  paxnew->Draw("");
  if ( pPadSave != nullptr ) pPadSave->cd();
  return 0;
}

//**********************************************************************

int TPadManipulator::addHistFun(unsigned int ifun) {
  m_histFuns.push_back(ifun);
  drawHistFuns();
  return 0;
}

//**********************************************************************

int TPadManipulator::drawHistFuns() {
  const TList* pfuns = nullptr;
  if      ( m_ph != nullptr ) pfuns = m_ph->GetListOfFunctions();
  else if ( m_pg != nullptr ) pfuns = m_pg->GetListOfFunctions();
  else    return 1;
  const TList& funs = *pfuns;
  unsigned int nfun = funs.GetEntries();
  for ( unsigned int ifun : m_histFuns ) {
    if ( ifun >= nfun ) continue;
    TF1* pfun = dynamic_cast<TF1*>(funs.At(ifun));
    pfun->Draw("same");
  }
  return 0;
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

int TPadManipulator::addVerticalModLines(double xmod, double xoff, double lenfrac) {
  m_vmlXmod = xmod;
  m_vmlXoff = xoff;
  m_vmlXLength = lenfrac;
  drawVerticalModLines();
  return 0;
}

//**********************************************************************

int TPadManipulator::drawVerticalModLines() {
  double xmod = m_vmlXmod;
  double xoff = m_vmlXoff;
  if ( xmod <= 0.0 ) return 1;
  double x = xoff;
  while ( x >= xmin() ) x -= xmod;
  while ( x < xmin() ) x += xmod;
  m_vmlLines.clear();
  double ytop = ymax() + (m_vmlXLength - 1.0)*(ymax() - ymin());
  while ( x <= xmax() ) {
    std::shared_ptr<TLine> pline(new TLine(x, ymin(), x, ytop));
    pline->SetLineStyle(3);
    m_vmlLines.push_back(pline);
    pline->Draw();
    x += xmod;
  }
  return 0;
}

//**********************************************************************

int TPadManipulator::fixFrameFillColor() {
  if ( gStyle == nullptr ) return 2;
  m_ppad->SetFrameFillColor(gStyle->GetColorPalette(0));
  return 0;
}

//**********************************************************************
