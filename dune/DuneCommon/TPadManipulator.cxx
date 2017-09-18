// TPadManipulator.cxx

#include "TPadManipulator.h"
#include "TPad.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TList.h"
#include "TStyle.h"
#include "TH1.h"

using std::string;

//**********************************************************************

TPadManipulator::TPadManipulator(TVirtualPad* ppad) : m_ppad(ppad) {
  if ( m_ppad == 0 ) m_ppad = gPad;
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
  auto prims = gPad->GetListOfPrimitives();
  for ( int iprm=0; iprm<prims->GetEntries(); ++iprm ) {
    TObject* pobj = prims->At(iprm);
    m_ph = dynamic_cast<TH1*>(pobj);
    if ( m_ph != nullptr ) break;
  }
  if ( pPadSave != nullptr ) pPadSave->cd();
  return 0;
}

//**********************************************************************

int TPadManipulator::addaxis() {
  return addaxistop() + addaxisright();
}

//**********************************************************************

int TPadManipulator::addaxistop() {
  if ( m_ph == nullptr ) return 1;
  double ticksize = 0;
  int ndiv = 0;
  TAxis* paxold = m_ph->GetXaxis();
  if ( paxold == nullptr ) return 2;
  ticksize = paxold->GetTickLength();
  ndiv = paxold->GetNdivisions();
  return addaxistop(ticksize, ndiv);
}

//**********************************************************************

int TPadManipulator::addaxistop(double ticksize, int ndiv) {
  TVirtualPad* pPadSave = m_ppad == gPad ? nullptr : gPad;
  if ( pPadSave != nullptr ) m_ppad->cd();
  string sopt = "-US";
  if ( gPad->GetLogx() ) sopt += "G";
  TGaxis* paxnew = new TGaxis(m_xminPad, m_ymaxPad, m_xmaxPad, m_ymaxPad,
                              m_xmin, m_xmax, 510, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickSize(ticksize);
  if ( ndiv > 0 ) paxnew->SetNdivisions(ndiv);
  string name = "TopAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = gPad->GetListOfPrimitives();
  for ( int iobj=0; iobj<pobjs->GetEntries(); ++iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold != nullptr ) continue;
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

int TPadManipulator::addaxisright() {
  if ( m_ph == nullptr ) return 1;
  double ticksize = 0;
  int ndiv = 0;
  TAxis* paxold = m_ph->GetYaxis();
  if ( paxold == nullptr ) return 2;
  ticksize = paxold->GetTickLength();
  ndiv = paxold->GetNdivisions();
  return addaxisright(ticksize, ndiv);
}

//**********************************************************************

int TPadManipulator::addaxisright(double ticksize, int ndiv) {
  TVirtualPad* pPadSave = m_ppad == gPad ? nullptr : gPad;
  if ( pPadSave != nullptr ) m_ppad->cd();
  if ( gPad == 0 ) return 1;
  gPad->Update();
  string sopt = "+US";
  if ( gPad->GetLogy() ) sopt += "G";
  TGaxis* paxnew = new TGaxis(m_xmaxPad, m_yminPad, m_xmaxPad, m_ymaxPad,
                              m_ymin, m_ymax, 510, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickLength(ticksize);
  if ( ndiv > 0 ) paxnew->SetNdivisions(ndiv);
  string name = "RightAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = gPad->GetListOfPrimitives();
  for ( int iobj=0; iobj<pobjs->GetEntries(); ++iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold != nullptr ) continue;
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

int TPadManipulator::fixFrameFillColor() {
  if ( gStyle == nullptr ) return 2;
  m_ppad->SetFrameFillColor(gStyle->GetColorPalette(0));
  return 0;
}

//**********************************************************************
