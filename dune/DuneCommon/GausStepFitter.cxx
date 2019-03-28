// GausStepFitter.cxx

#include "GausStepFitter.h"
#include "dune/DuneCommon/gausTF1.h"
#include "dune/DuneCommon/quietHistFit.h"
#include "TF1.h"
#include "TH1.h"
#include <iostream>

using std::cout;
using std::endl;

//**********************************************************************

GausStepFitter::GausStepFitter(double pos, double sigma, double height, Name fnam, Name fopt)
: m_pos(pos), m_sigma(sigma), m_height(height), m_fnam(fnam), m_fopt(fopt),
  m_LogLevel(0) { }

//**********************************************************************

int GausStepFitter::fit(TH1* ph) const {
  Name myname = "GausStepFitter::fit: ";
  double sigma = m_sigma;
  double height = m_height;
  double sigfac = 2.0;
  for ( int ifit=0; ifit<5; ++ifit ) {
    TF1* pffix = gausTF1(height, m_pos, sigma, m_fnam);
    double sigmax = 1.1*sigfac*sigma;
    double sigmin = 0.9*sigma/sigfac;
    if ( m_LogLevel >= 1 ) cout << myname << "  Doing constrained fit " << ifit
                                << " with pos=" << m_pos << ", sigma=" << sigma
                                << " (" << sigmin << ", " << sigmax << ")" << endl;
    pffix->SetParameter(2, sigma);
    pffix->SetParLimits(2, sigmin, sigmax);
    if ( m_height == 0.0 ) {
      pffix->SetParLimits(0, 0.0, 1.e10);   // Don't let height go negative.
    } else if ( height > 0.0 ) {
      pffix->SetParLimits(0, 0.1*height, 2.0*height);   // Don't let height go negative.
    } else if ( height < 0.0 ) {
      pffix->SetParLimits(0, 2.0*height, 0.1*height);   // Don't let height go positive.
    }
    Name fopt = m_fopt;
    if ( m_LogLevel >= 2 ) fopt += "Q";
    int fstat = quietHistFit(ph, pffix, fopt.c_str());
    double signew = pffix->GetParameter(2);
    bool atHiLimit = signew > 0.999*sigmax;
    bool atLoLimit = signew < 1.001*sigmin;
    if ( m_LogLevel >= 1 ) cout << myname << "  status " << fstat << ", fit sigma=" << signew  << endl;
    if ( fstat == 0 && !atHiLimit && !atLoLimit ) {
      if ( m_LogLevel >= 1 ) cout << myname << "Fit suceeded." << endl;
      return 0;
    }
    // Otherwise we may get a crash when we try to view saved copy of histo:
    ph->GetListOfFunctions()->Clear();
    delete pffix;
    if ( atLoLimit ) sigma /= sigfac;
    else             sigma *= sigfac;
  }
  if ( m_LogLevel >= 1 ) cout << myname << "Fit failed." << endl;
  return 1;
}

//**********************************************************************
