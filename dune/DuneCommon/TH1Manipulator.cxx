// TH1Manipulator.cxx

#include "TH1Manipulator.h"
#include "TH1.h"
#include "TPad.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TList.h"
#include "TStyle.h"

using std::string;

//**********************************************************************

int TH1Manipulator::addaxis(TH1* ph) {
  addaxistop(ph);
  return addaxisright(ph);
}

//**********************************************************************

int TH1Manipulator::addaxistop(TH1* ph) {
  double ticksize = 0;
  int ndiv = 0;
  if ( ph != 0 ) {
    TAxis* paxold = ph->GetXaxis();
    if ( paxold != 0 ) {
      ticksize = paxold->GetTickLength();
      ndiv = paxold->GetNdivisions();
    }
  }
  return addaxistop(ticksize, ndiv);
}

//**********************************************************************

int TH1Manipulator::addaxistop(double ticksize, int ndiv) {
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  if ( gPad == 0 ) return 3;
  gPad->Update();
  gPad->GetRangeAxis(xmin, ymin, xmax, ymax);
  double xax1 = gPad->GetUxmin();
  double xax2 = gPad->GetUxmax();
  double yax = gPad->GetUymax();
  string sopt = "-US";
  if ( gPad->GetLogx() ) {
    xax1 = pow(10.0, xax1);
    xax2 = pow(10.0, xax2);
    xmin = pow(10.0, xmin);
    xmax = pow(10.0, xmax);
    sopt += "G";
  }
  if ( gPad->GetLogy() ) {
    yax = pow(10.0, yax);
  }
  TGaxis* paxnew = new TGaxis(xax1, yax, xax2, yax,
                              xmin, xmax, 510, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickSize(ticksize);
  if ( ndiv > 0 ) paxnew->SetNdivisions(ndiv);
  string name = "TopAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = gPad->GetListOfPrimitives();
  for ( int iobj=0; iobj<pobjs->GetEntries(); ++iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold != 0 ) {
      
      if ( paxold->GetName() == name ) {
        pobjs->RemoveAt(iobj);
        break;
      }
    }
  }
  paxnew->Draw();
  return 0;
}

//**********************************************************************

int TH1Manipulator::addaxisright(TH1* ph) {
  double ticksize = 0;
  int ndiv = 0;
  if ( ph != 0 ) {
    TAxis* paxold = ph->GetYaxis();
    if ( paxold != 0 ) {
      ticksize = paxold->GetTickLength();
      ndiv = paxold->GetNdivisions();
    }
  }
  return addaxisright(ticksize, ndiv);
}

//**********************************************************************

int TH1Manipulator::addaxisright(double ticksize, int ndiv) {
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  if ( gPad == 0 ) return 3;
  gPad->Update();
  gPad->GetRangeAxis(xmin, ymin, xmax, ymax);
  double xax = gPad->GetUxmax();
  double yax1 = gPad->GetUymin();
  double yax2 = gPad->GetUymax();
  string sopt = "+US";
  if ( gPad->GetLogx() ) {
    xax = pow(10.0, xax);
  }
  if ( gPad->GetLogy() ) {
    yax1 = pow(10.0, yax1);
    yax2 = pow(10.0, yax2);
    ymin = pow(10.0, ymin);
    ymax = pow(10.0, ymax);
    sopt += "G";
  }
  if ( gPad->GetLogy() ) sopt += "G";
  TGaxis* paxnew = new TGaxis(xax, yax1, xax, yax2,
                              ymin, ymax, 510, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickLength(ticksize);
  if ( ndiv > 0 ) paxnew->SetNdivisions(ndiv);
  string name = "RightAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = gPad->GetListOfPrimitives();
  for ( int iobj=0; iobj<pobjs->GetEntries(); ++iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold != 0 ) {
      if ( paxold->GetName() == name ) {
        pobjs->RemoveAt(iobj);
        break;
      }
    }
  }
  paxnew->Draw("");
  return 0;
}

//**********************************************************************

int TH1Manipulator::fixFrameFillColor() {
  if ( gPad == nullptr ) return 1;
  if ( gStyle == nullptr ) return 2;
  gPad->SetFrameFillColor(gStyle->GetColorPalette(0));
  return 0;
}

//**********************************************************************
//**********************************************************************
