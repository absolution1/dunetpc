// shiftHistFit.cxx

#include "dune/DuneCommon/shiftHistFit.h"
#include "quietHistFit.h"
#include "TFitResult.h"

#include "TF1.h"
#include "TH1F.h"

//**********************************************************************

int shiftHistFit(TH1* ph, TF1* pf, std::string fopt, int ipar, double xshift) {
  if ( ipar < 0 || ipar >= pf->GetNpar() ) return 101;
  // Shift the histogram
  TAxis* pax = ph->GetXaxis();
  double xh1Min = pax->GetXmin();
  double xh1Max = pax->GetXmax();
  double xh2Min = xh1Min - xshift;
  double xh2Max = xh1Max - xshift;
  pax->SetLimits(xh2Min, xh2Max);
  // Shift the function.
  double xf1Min = pf->GetXmin();
  double xf1Max = pf->GetXmax();
  double xf1Par = pf->GetParameter(ipar);
  double xf1Llo = 0.0;
  double xf1Lhi = 0.0;
  pf->GetParLimits(ipar, xf1Llo, xf1Lhi);
  bool haveLimits = xf1Llo || xf1Lhi;
  double xf2Min = xf1Min - xshift;
  double xf2Max = xf1Max - xshift;
  double xf2Par = xf1Par - xshift;
  double xf2Llo = xf1Llo - xshift;
  double xf2Lhi = xf1Lhi - xshift;
  pf->SetRange(xf2Min, xf2Max);
  pf->SetParameter(ipar, xf2Par);
  if ( haveLimits ) pf->SetParLimits(ipar, xf2Llo, xf2Lhi);
  // Perform fit.
  int fstat = quietHistFit(ph, pf, fopt);
  // Shift function back.
  xf2Par = pf->GetParameter(ipar);
  xf1Min = xf2Min + xshift;
  xf1Max = xf2Max + xshift;
  xf1Par = xf2Par + xshift;
  pf->SetRange(xf1Min, xf1Max);
  pf->SetParameter(ipar, xf1Par);
  if ( haveLimits ) pf->SetParLimits(ipar, xf1Llo, xf1Lhi);
  // Shift histogram back.
  pax->SetLimits(xh1Min, xh1Max);
  return fstat;
}

//**********************************************************************

int shiftHistFit(TH1* ph, TF1* pf, std::string fopt, std::string spar, double xshift) {
  int ipar = pf->GetParNumber(spar.c_str());
  return shiftHistFit(ph, pf, fopt, ipar, xshift);
}

//**********************************************************************
