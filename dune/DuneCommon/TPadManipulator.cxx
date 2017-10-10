// TPadManipulator.cxx

#include "TPadManipulator.h"
#include "TPad.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TList.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLine.h"
#include "TF1.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

TPadManipulator::TPadManipulator(TVirtualPad* ppad)
: m_ppad(ppad),
  m_top(false), m_right(false),
  m_vmlXmod(0.0), m_vmlXoff(0.0) {
  if ( m_ppad == 0 ) m_ppad = gPad;
  update();
  cout << "Constructed @" << this << endl;
}

//**********************************************************************

TPadManipulator::~TPadManipulator() {
  cout << "Deleting " << m_vmlLines.size() << " lines @" << this << endl;
}

//**********************************************************************

int TPadManipulator::update() {
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
  m_ph = nullptr;
  for ( int iprm=0; iprm<prims->GetEntries(); ++iprm ) {
    TObject* pobj = prims->At(iprm);
    m_ph = dynamic_cast<TH1*>(pobj);
    if ( m_ph != nullptr ) break;
  }
  if ( m_top ) drawAxisTop();
  if ( m_right ) drawAxisRight();
  if ( m_vmlXmod >= 0.0 ) drawVerticalModLines();
  if ( pPadSave != nullptr ) pPadSave->cd();
  return 0;
}

//**********************************************************************

int TPadManipulator::setRangeX(double x1, double x2) {
  TH1* ph = hist();
  if ( ph == nullptr ) return 1;
  if ( x2 <= x1 ) return 2;
  ph->GetXaxis()->SetRangeUser(x1, x2);
  ph->Draw();
  return update();
}

//**********************************************************************

int TPadManipulator::setRangeY(double y1, double y2) {
  TH1* ph = hist();
  if ( ph == nullptr ) return 1;
  if ( y2 <= y1 ) return 2;
  ph->GetYaxis()->SetRangeUser(y1, y2);
  ph->Draw();
  return update();
}

//**********************************************************************

int TPadManipulator::setRanges(double x1, double x2, double y1, double y2) {
  TH1* ph = hist();
  if ( ph == nullptr ) return 1;
  if ( x2 <= x1 ) return 2;
  if ( y2 <= y1 ) return 2;
  ph->GetXaxis()->SetRangeUser(x1, x2);
  ph->GetYaxis()->SetRangeUser(y1, y2);
  ph->Draw();
  return update();
}

//**********************************************************************

int TPadManipulator::addAxis() {
  return addAxisTop() + addAxisRight();
}

//**********************************************************************

int TPadManipulator::addAxisTop() {
  if ( m_ph == nullptr ) return 1;
  double ticksize = 0;
  int ndiv = 0;
  TAxis* paxold = m_ph->GetXaxis();
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
  if ( m_ph == nullptr ) return 1;
  double ticksize = 0;
  int ndiv = 0;
  TAxis* paxold = m_ph->GetYaxis();
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
  if ( m_ph == nullptr ) return 0;
  const TList& funs =  *m_ph->GetListOfFunctions();
  unsigned int nfun = funs.GetEntries();
  for ( unsigned int ifun : m_histFuns ) {
    if ( ifun >= nfun ) continue;
    TF1* pfun = dynamic_cast<TF1*>(funs.At(ifun));
    pfun->Draw("same");
  }
  return 0;
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
